"""Fit the application grid refinements and EURUSD marginal transformations."""

from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
from numba import njit
from scipy.fft import rfft
from scipy.optimize import minimize
from scipy.stats import rankdata
from cse import Bank, features
from check_numerics import interpolated

ROOT = Path(__file__).resolve().parents[1]
LOWER = np.array([.005, .01])
SPAN = np.array([.995, 3.])-LOWER


def read_returns(path, series):
    data = pd.read_csv(path, na_values=["."])
    date = next(x for x in data if x.lower() in ["date", "observation_date"])
    values = pd.to_numeric(data[series], errors="coerce")
    data = pd.DataFrame({"date":pd.to_datetime(data[date]), "price":values}).dropna()
    data = data.sort_values("date")
    data = data[(data.date >= "2016-12-30") & (data.date <= "2025-12-31")]
    if data.date.duplicated().any() or (data.price <= 0).any():
        raise ValueError("Prices must be positive and dates unique.")
    if data.empty or data.date.iloc[0] != pd.Timestamp("2016-12-30") or data.date.iloc[-1] != pd.Timestamp("2025-12-31"):
        raise ValueError("The price file must cover 30 December 2016 through 31 December 2025.")
    returns = 100*np.diff(np.log(data.price.to_numpy()))
    dates = data.date.iloc[1:].to_numpy()
    selected = dates >= np.datetime64("2017-01-03")
    answer = returns[selected]
    expected = {"DEXUSEU":2246,"SP500":2262}[series]
    if len(answer) != expected:
        raise ValueError(f"Expected {expected} returns for {series}; found {len(answer)}.")
    return answer


def minimize_criterion(fun, extra=None):
    starts = np.array([(p,v) for v in [.1,.5,.9] for p in [.1,.5,.9]])
    if extra is not None:
        starts = np.vstack([starts, np.clip((np.asarray(extra)-LOWER)/SPAN,0,1)])
    values = np.array([fun(LOWER+z*SPAN) for z in starts])
    scale = max(1e-9, np.median(values[np.isfinite(values)]))
    def objective(z):
        return fun(LOWER+z*SPAN)/scale
    def gradient(z):
        answer = np.empty(2)
        for j in range(2):
            a, b = z.copy(), z.copy()
            a[j] = max(0,a[j]-2e-5)
            b[j] = min(1,b[j]+2e-5)
            answer[j] = (objective(b)-objective(a))/(b[j]-a[j])
        return answer
    fits = [minimize(objective, starts[i], jac=gradient, method="L-BFGS-B",
                     bounds=[(0,1),(0,1)], options={"maxiter":300,"ftol":1e7*np.finfo(float).eps,"gtol":1e-7})
            for i in np.argsort(values,kind="stable")[:4]]
    fit = min(fits,key=lambda f:f.fun)
    return {"theta":LOWER+fit.x*SPAN,"objective":float(fit.fun*scale),"converged":bool(fit.success)}


def subset_bank(bank,tau):
    indices = np.searchsorted(bank.tau,np.round(tau,12))
    if not np.allclose(bank.tau[indices],tau,atol=1e-12,rtol=0):
        raise ValueError("Missing coefficient levels.")
    obj = Bank.__new__(Bank)
    obj.__dict__.update(bank.__dict__)
    obj.values = bank.values[:,indices]
    obj.derivative = bank.derivative[:,indices]
    obj.tau = tau
    return obj


def accelerated_objective(y,spec,bank):
    n, q, length = len(y), len(spec.tau), len(spec.lam)
    kk = np.floor(n*spec.lam/(2*np.pi)+1e-12).astype(int)
    scores = rankdata(y,method="max")/n
    d = rfft((scores[:,None] <= spec.tau).astype(float),axis=0)[1:kk.max()+1]
    dr, di = d.real, d.imag
    g0 = np.minimum.outer(spec.tau,spec.tau)-np.outer(spec.tau,spec.tau)
    cumulative = np.zeros((q,q))
    previous, constant = 0, 0.
    for k, lam in zip(kk,spec.lam):
        cumulative += dr[previous:k].T@dr[previous:k]+di[previous:k].T@di[previous:k]
        target = cumulative/n**2-lam*g0/(2*np.pi)
        constant += np.sum(target*target)/(length*q*q)
        previous = k
    def objective(theta):
        b = interpolated(bank,theta[1])
        r = theta[0]**np.arange(1,b.shape[1]+1)
        envelope = np.max(np.abs(b),axis=0)**2*r/(np.pi*(1-r))
        remaining = np.r_[np.cumsum(envelope[::-1])[::-1][1:],0.]
        keep = int(np.flatnonzero(remaining <= 1e-13)[0])+1
        b, r = b[:,:keep],r[:keep]
        angles = np.arctan(np.sin(spec.lam[:,None])*r/(1-np.cos(spec.lam[:,None])*r))/np.pi
        z = (dr@b)**2+(di@b)**2
        cumulative_z = np.vstack([np.zeros(keep),np.cumsum(z,axis=0)])
        h = cumulative_z[kk]/n**2-np.outer(spec.lam/(2*np.pi),np.sum(b*(g0@b),axis=0))
        gram = (b.T@b)**2
        return constant+(np.sum((angles@gram)*angles)-2*np.sum(angles*h))/(length*q*q)
    return objective


def relative_target_difference(theta,reference,spec,bank):
    b0,b1 = interpolated(bank,reference[1]),interpolated(bank,theta[1])
    r0 = reference[0]**np.arange(1,b0.shape[1]+1)
    r1 = theta[0]**np.arange(1,b1.shape[1]+1)
    def angles(r):
        return np.arctan(np.sin(spec.lam[:,None])*r/(1-np.cos(spec.lam[:,None])*r))/np.pi
    a0,a1 = angles(r0),angles(r1)
    n0 = np.sum((a0@((b0.T@b0)**2))*a0)
    n1 = np.sum((a1@((b1.T@b1)**2))*a1)
    cross = np.sum((a0@((b0.T@b1)**2))*a1)
    return np.sqrt(max(0,n0+n1-2*cross)/n0)


@njit(cache=True)
def qml_value(theta,x):
    phi,v = theta
    measurement = np.pi**2/2
    innovation = v*(1-phi**2)
    a,b,P,srr,srd,sdd,ld = 0.,0.,v,0.,0.,0.,0.
    for xt in x:
        variance = P+measurement
        residual = xt-a
        derivative = 1+b
        srr += residual*residual/variance
        srd += residual*derivative/variance
        sdd += derivative*derivative/variance
        ld += np.log(variance)
        gain = P/variance
        a = phi*(a+gain*residual)
        b = phi*(b-gain*derivative)
        P = phi**2*(P-gain*P)+innovation
    return (ld+srr-srd*srd/sdd)/len(x)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir",type=Path,required=True)
    parser.add_argument("--output",type=Path,default=ROOT/"results"/"sensitivity")
    parser.add_argument("--qml-only",action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    eurusd = read_returns(args.data_dir/"DEXUSEU.csv","DEXUSEU")
    x = (eurusd-eurusd.mean())/eurusd.std(ddof=1)
    transformations = {"identity":x,"asinh":np.arcsinh(x),"cubic":x+.2*x**3,
                       "skew":np.sinh(np.arcsinh(x)+.8)}
    margins = []
    for name,y in transformations.items():
        centered = y-y.mean()
        if np.any(centered == 0):
            raise ValueError("QML needs nonzero centered observations.")
        fit = minimize_criterion(lambda theta:qml_value(theta,2*np.log(np.abs(centered))))
        margins.append({"margin":name,"phi":fit["theta"][0],"v":fit["theta"][1],
                        "objective":fit["objective"],"converged":fit["converged"]})
    pd.DataFrame(margins).to_csv(args.output/"EURUSD_margin_qml.csv",index=False)
    print(pd.DataFrame(margins).to_string(index=False),flush=True)
    if args.qml_only:
        return
    series = {"EURUSD":eurusd,"SP500":read_returns(args.data_dir/"SP500.csv","SP500")}
    grids = [(4,4),(8,8),(16,16),(16,32),(32,32),(16,64)]
    tau = np.unique(np.concatenate([features(len(y),q,f).tau for y in series.values() for q,f in grids]))
    bank = Bank(tau,upper_phi=.995)
    output = []
    for name,y in series.items():
        point = pd.read_csv(ROOT/"results"/"reported"/f"{name}_point_estimate.csv").iloc[0]
        reference = point[["phi","v"]].to_numpy(dtype=float)
        ref_spec = features(len(y),16,32)
        ref_bank = subset_bank(bank,ref_spec.tau)
        fits = []
        for q,f in grids:
            spec = features(len(y),q,f)
            reduced = subset_bank(bank,spec.tau)
            fit = minimize_criterion(accelerated_objective(y,spec,reduced),extra=reference)
            fits.append((spec,fit))
            print(name,len(spec.tau),len(spec.lam),fit["theta"],flush=True)
        fitted_reference = fits[3][1]["theta"]
        for spec,fit in fits:
            difference = relative_target_difference(fit["theta"],fitted_reference,ref_spec,ref_bank)
            output.append({"application":name,"n":len(y),"N_u":len(spec.tau),"N_lambda":len(spec.lam),
                           "phi":fit["theta"][0],"v":fit["theta"][1],"objective":fit["objective"],
                           "converged":fit["converged"],"target_difference_percent":100*difference})
    pd.DataFrame(output).to_csv(args.output/"grid_summary.csv",index=False)
    settings = {"parameter_region":[LOWER.tolist(),(LOWER+SPAN).tolist()],
                "grid_multipliers":grids,"Hermite_tolerance":1e-5,"interpolation_nodes":81,
                "sample_sizes":{k:len(v) for k,v in series.items()},"return_formula":"100 * diff(log(price))"}
    (args.output/"settings.json").write_text(json.dumps(settings,indent=2)+"\n")


if __name__ == "__main__":
    main()
