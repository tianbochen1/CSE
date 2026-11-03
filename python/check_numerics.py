"""Check spectral evaluation on the full application quantile grid."""

from pathlib import Path
import json
import numpy as np
import pandas as pd
from cse import Bank, coefficients, gaussian_rule


def hermite_rule(order, terms):
    x, w = gaussian_rule(order)
    h = np.empty((len(x), terms + 1))
    h[:, 0] = np.sqrt(w)
    h[:, 1] = x * h[:, 0]
    for j in range(1, terms):
        h[:, j + 1] = (x*h[:, j]-np.sqrt(j)*h[:, j-1])/np.sqrt(j+1)
    return x, w, h[:, 1:]


def model_matrices(b, phi, endpoints):
    powers = phi ** np.arange(1, b.shape[1] + 1)
    return np.array([(b*(np.arctan(powers*np.sin(lam)/(1-powers*np.cos(lam)))/np.pi)) @ b.T
                     for lam in endpoints])


def interpolated(bank, v):
    step = bank.s[1]-bank.s[0]
    i = min(len(bank.s)-2, max(0, int((np.sqrt(v)-bank.s[0])/step)))
    a = (np.sqrt(v)-bank.s[i])/step
    return ((2*a**3-3*a*a+1)*bank.values[i]+(a**3-2*a*a+a)*step*bank.derivative[i]
            +(-2*a**3+3*a*a)*bank.values[i+1]+(a**3-a*a)*step*bank.derivative[i+1])


def main():
    root = Path(__file__).resolve().parents[1]
    tau = np.round(.05+.9*(np.arange(128)+.5)/128, 12)
    endpoints = np.pi*np.array([.05,.2,.5,.8,.95])
    phi = .995
    bank = Bank(tau, upper_phi=phi)
    direct_rule = hermite_rule(bank.order, bank.jmax)
    fine_rule = hermite_rule(bank.order+128, bank.jmax)
    fine_terms = int(np.ceil(np.log(4*np.pi*1e-6/(1+4*np.pi*1e-6))/np.log(phi)-1))
    trunc_order = 64*int(np.ceil((fine_terms+24)/64))
    trunc_rule = hermite_rule(trunc_order, fine_terms)
    positions = np.unique(np.rint(np.linspace(0, len(bank.s)-2, 17)).astype(int))
    variances = (bank.s[positions]+.37*(bank.s[1]-bank.s[0]))**2
    records = []
    for v in variances:
        exact = coefficients(tau, v, *direct_rule)
        error = np.max(np.abs(model_matrices(interpolated(bank,v),phi,endpoints)
                              -model_matrices(exact,phi,endpoints)))
        records.append({"comparison":"interpolation", "v":v, "max_error":error})
    for v in [.03,.25,1.,2.,2.9]:
        b = coefficients(tau,v,*direct_rule)
        bf = coefficients(tau,v,*fine_rule)
        error = np.max(np.abs(model_matrices(b,phi,endpoints)-model_matrices(bf,phi,endpoints)))
        records.append({"comparison":"Gaussian quadrature", "v":v, "max_error":error})
        b = coefficients(tau,v,*trunc_rule)
        error = np.max(np.abs(model_matrices(b[:,:bank.jmax],phi,endpoints)-model_matrices(b,phi,endpoints)))
        records.append({"comparison":"Hermite truncation", "v":v, "max_error":error})
    result = pd.DataFrame(records)
    result.to_csv(root/"results/application_numerical_accuracy.csv", index=False)
    settings = {"upper_phi":phi, "quantile_count":len(tau), "endpoints_over_pi":(endpoints/np.pi).tolist(),
                "interpolation_nodes":len(bank.s), "Hermite_orders":[bank.jmax,fine_terms],
                "Gaussian_orders":[bank.order,bank.order+128,trunc_order],
                "maximum_errors":result.groupby("comparison").max_error.max().to_dict()}
    (root/"results/application_numerical_settings.json").write_text(json.dumps(settings,indent=2)+"\n")
    print(json.dumps(settings,indent=2))


if __name__ == "__main__":
    main()
