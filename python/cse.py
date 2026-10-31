"""Copula-spectral estimation for Gaussian stochastic volatility."""

from dataclasses import dataclass
from pathlib import Path
import math
import numpy as np
from numba import njit
from scipy.fft import rfft
from scipy.optimize import brentq, minimize
from scipy.signal import lfilter
from scipy.special import ndtr, ndtri, roots_hermitenorm
from scipy.stats import rankdata


@dataclass
class Features:
    tau: np.ndarray
    lam: np.ndarray
    pairs: np.ndarray
    weights: np.ndarray


def features(n, q_multiplier=1, f_multiplier=1, diagnostic=False):
    if diagnostic:
        tau = np.array([.1, .2, .35, .5, .65, .8, .9])
        count = 32
    else:
        k = int(np.ceil(np.log(n))) * q_multiplier
        tau = np.round(.05 + .9 * (np.arange(k) + .5) / k, 12)
        count = int(np.ceil(4 * np.log(n))) * f_multiplier
    lam = np.pi * (np.arange(count) + .5) / count
    pairs = np.array([(i, j) for j in range(len(tau)) for i in range(j + 1)])
    weights = np.repeat(np.where(pairs[:, 0] == pairs[:, 1], 1., 2.), count)
    return Features(tau, lam, pairs, weights / weights.sum())


def gaussian_rule(order):
    x, w = roots_hermitenorm(order)
    take = w > 0
    return x[take], w[take] / w[take].sum()


def coefficients(tau, v, x, w, hermite):
    scale = np.exp(-np.sqrt(v) * x / 2 + v / 4)
    q = []
    for u in tau:
        if abs(u - .5) < 1e-13:
            q.append(0.)
            continue
        p = min(u, 1 - u)
        hi = max(1., abs(ndtri(p)))
        fun = lambda z: np.dot(w, ndtr(-z * scale)) - p
        while fun(hi) > 0:
            hi *= 2
        z = brentq(fun, 0, hi, xtol=1e-12)
        q.append(-z if u < .5 else z)
    values = ndtr(scale[:, None] * np.array(q)) - tau
    out = (values * np.sqrt(w)[:, None]).T @ hermite
    out[np.abs(tau - .5) < 1e-13] = 0
    return out


class Bank:
    def __init__(self, tau, upper_phi=.995, tolerance=1e-5, nodes=81):
        self.tau = np.sort(np.unique(np.round(tau, 12)))
        self.s = np.linspace(.1, np.sqrt(3), nodes)
        self.upper_phi = upper_phi
        self.tolerance = tolerance
        jmax = max(1, int(np.ceil(np.log(4*np.pi*tolerance/(1+4*np.pi*tolerance)) /
                                  np.log(upper_phi) - 1)))
        self.order = max(512, int(64 * np.ceil((jmax + 24) / 64)))
        x, w = gaussian_rule(self.order)
        h = np.empty((len(x), jmax + 1))
        h[:, 0] = np.sqrt(w)
        h[:, 1] = x * h[:, 0]
        for j in range(1, jmax):
            h[:, j + 1] = (x*h[:, j] - np.sqrt(j)*h[:, j - 1]) / np.sqrt(j + 1)
        values = np.array([coefficients(self.tau, s*s, x, w, h[:, 1:]) for s in self.s])
        step = self.s[1] - self.s[0]
        derivative = np.gradient(values, step, axis=0, edge_order=2)
        envelope = np.max(np.abs(values), axis=(0, 1)) + 8*step/27*np.max(np.abs(derivative), axis=(0, 1))
        r = upper_phi ** np.arange(1, jmax + 1)
        bound = envelope**2 * r / (np.pi*(1-r))
        remainder = np.r_[np.cumsum(bound[::-1])[::-1][1:], 0.]
        keep = int(np.flatnonzero(remainder <= 1e-13)[0]) + 1
        self.values = np.ascontiguousarray(values[:, :, :keep])
        self.derivative = np.ascontiguousarray(derivative[:, :, :keep])
        self.jmax = jmax
        self.removed_bound = float(remainder[keep-1])
        self.truncation_bound = float(upper_phi**(jmax+1)/(4*np.pi*(1-upper_phi**(jmax+1))))
        if not np.isfinite(self.values).all():
            raise FloatingPointError("Nonfinite model coefficients")

    def save(self, path):
        np.savez_compressed(path, **self.__dict__)

    @classmethod
    def load(cls, path):
        obj = cls.__new__(cls)
        for key, value in np.load(path).items():
            setattr(obj, key, value.item() if value.ndim == 0 else value)
        return obj


@njit(cache=True)
def model_jacobian(theta, s, values, deriv, lam, pairs):
    phi, v = theta
    root = np.sqrt(v)
    step = s[1] - s[0]
    cell = min(len(s)-2, max(0, int((root-s[0])/step)))
    a = (root-s[cell])/step
    h0, h1 = 2*a**3-3*a*a+1, -2*a**3+3*a*a
    g0, g1 = a**3-2*a*a+a, a**3-a*a
    dh0, dh1 = (6*a*a-6*a)/step, (-6*a*a+6*a)/step
    dg0, dg1 = 3*a*a-4*a+1, 3*a*a-2*a
    b = h0*values[cell] + g0*step*deriv[cell] + h1*values[cell+1] + g1*step*deriv[cell+1]
    db = (dh0*values[cell] + dg0*deriv[cell] + dh1*values[cell+1] + dg1*deriv[cell+1]) / (2*root)
    length = len(lam)
    out = np.zeros(length*len(pairs))
    jac = np.zeros((len(out), 2))
    keep = b.shape[1]
    remainder = 0.
    for j in range(b.shape[1]-1, -1, -1):
        bb = np.max(np.abs(b[:, j]))
        dd = np.max(np.abs(db[:, j]))
        r = phi**(j+1)
        rp = (j+1)*phi**j
        bound = bb*bb*r/(np.pi*(1-r))
        bound += bb*bb*rp/(np.pi*(1-r)**2) + 2*bb*dd*r/(np.pi*(1-r))
        if remainder + bound > 1e-12:
            break
        remainder += bound
        keep = j
    sn_all, cs_all = np.sin(lam), np.cos(lam)
    for j in range(keep):
        r = phi**(j+1)
        rp = (j+1)*phi**j
        for l in range(length):
            sn, cs = sn_all[l], cs_all[l]
            angle = np.arctan(r*sn/(1-r*cs))/np.pi
            dangle = rp*sn/(np.pi*(1-2*r*cs+r*r))
            for p in range(len(pairs)):
                i, k = pairs[p]
                ix = p*length + l
                prod = b[i, j]*b[k, j]
                out[ix] += angle*prod
                jac[ix, 0] += dangle*prod
                jac[ix, 1] += angle*(db[i, j]*b[k, j]+b[i, j]*db[k, j])
    return out, jac


def empirical(y, spec):
    y = np.atleast_2d(y)
    n = y.shape[1]
    ranks = rankdata(y, method="max", axis=1) / n
    indicators = ranks[:, :, None] <= spec.tau
    d = rfft(indicators.astype(float), axis=1, workers=1)
    kk = np.floor(n*spec.lam/(2*np.pi)+1e-12).astype(int)
    out = np.empty((len(y), len(spec.weights)), dtype=complex)
    L = len(spec.lam)
    for p, (i, j) in enumerate(spec.pairs):
        products = d[:, 1:, i] * d[:, 1:, j].conj()
        cumulative = np.c_[np.zeros(len(y)), np.cumsum(products, axis=1)] / n**2
        g0 = min(spec.tau[i], spec.tau[j]) - spec.tau[i]*spec.tau[j]
        out[:, p*L:(p+1)*L] = cumulative[:, kk] - spec.lam*g0/(2*np.pi)
    return out


def generate(theta, n, rng, count=1):
    phi, v = theta
    noise = rng.normal(size=(count, n)) * np.sqrt(v*(1-phi*phi))
    initial = rng.normal(size=(count, 1))*np.sqrt(v)*phi
    h, _ = lfilter([1.], [1., -phi], noise, axis=1, zi=initial)
    return np.exp(h/2-v/4) * rng.normal(size=(count, n))


class Estimator:
    def __init__(self, spec, bank, upper_phi=.995):
        self.spec = spec
        self.bank = bank
        indices = np.searchsorted(bank.tau, np.round(spec.tau, 12))
        if not np.allclose(bank.tau[indices], spec.tau, atol=1e-12, rtol=0):
            raise ValueError("Missing quantile levels")
        self.b = np.ascontiguousarray(bank.values[:, indices])
        self.db = np.ascontiguousarray(bank.derivative[:, indices])
        self.lower = np.array([.005, .01])
        self.span = np.array([upper_phi, 3.]) - self.lower
        self.starts = np.array([(i, j) for i in [.1, .5, .9] for j in [.1, .5, .9]])
        self.start_models = np.array([self.model(self.lower+x*self.span)[0] for x in self.starts])

    def dense_starts(self, count=25):
        phi = np.r_[np.linspace(.005, .95, count-5), .96, .975, .985, .99, .995]
        variance = np.linspace(.1, np.sqrt(3), count)**2
        points = np.array([(a, b) for a in phi for b in variance])
        self.starts = (points-self.lower)/self.span
        self.start_models = np.array([self.model(t)[0] for t in points])

    def batch_fit(self, targets):
        targets = targets.real
        w = self.spec.weights
        costs = np.sum(self.start_models**2*w, axis=1)[:, None] - 2*(self.start_models*w)@targets.T
        indices = np.argmin(costs, axis=0)
        answers = []
        for target, index in zip(targets, indices):
            scale = max(1e-9, np.mean(target**2))
            def objective(x):
                model, jac = self.model(self.lower+x*self.span)
                residual = model-target
                wr = w*residual
                return np.dot(wr, residual)/scale, 2*(wr@jac)*self.span/scale
            opt = minimize(objective, self.starts[index], method="L-BFGS-B", jac=True,
                           bounds=[(0,1),(0,1)], options={"maxiter":300,"ftol":1e7*np.finfo(float).eps,"gtol":1e-7})
            answers.append({"theta":self.lower+opt.x*self.span,"objective":float(opt.fun*scale),
                            "boundary":bool(np.any((opt.x<=.002)|(opt.x>=.998))),"converged":bool(opt.success)})
        return answers

    def model(self, theta):
        return model_jacobian(np.asarray(theta, dtype=float), self.bank.s, self.b, self.db,
                              self.spec.lam, self.spec.pairs)

    def fit(self, empirical_values, extra=None, nstarts=4):
        target = np.asarray(empirical_values).real
        starts = self.starts
        values = np.sum((self.start_models-target)**2*self.spec.weights, axis=1)
        if extra is not None:
            starts = np.vstack([starts, np.clip((extra-self.lower)/self.span, 0, 1)])
            values = np.r_[values, np.sum((self.model(extra)[0]-target)**2*self.spec.weights)]
        scale = max(1e-9, float(np.median(values)))
        def objective(x):
            model, jac = self.model(self.lower+x*self.span)
            residual = model-target
            weighted = self.spec.weights*residual
            return np.dot(weighted, residual)/scale, 2*(weighted@jac)*self.span/scale
        results = [minimize(objective, starts[i], method="L-BFGS-B", jac=True,
                            bounds=[(0, 1), (0, 1)],
                            options={"maxiter":300, "ftol":1e7*np.finfo(float).eps, "gtol":1e-7})
                   for i in np.argsort(values)[:nstarts]]
        best = min(results, key=lambda x:x.fun)
        theta = self.lower+best.x*self.span
        boundary = bool(np.any((best.x <= .002) | (best.x >= .998)))
        return {"theta":theta, "objective":float(best.fun*scale),
                "boundary":boundary, "converged":bool(best.success)}


def bootstrap(y, estimator, diagnostic_estimator, rng, draws=499):
    observed = empirical(y, estimator.spec)[0]
    fit = estimator.fit(observed)
    theta = fit["theta"]
    simulated = generate(theta, len(y), rng, draws)
    em = empirical(simulated, estimator.spec)
    fits = estimator.batch_fit(em)
    estimates = np.array([z["theta"] for z in fits])
    limits = 2*theta-np.quantile(estimates, [.975, .025], axis=0)
    limits[0] = np.maximum(limits[0], estimator.lower)
    limits[1] = np.minimum(limits[1], estimator.lower+estimator.span)
    percentile = np.quantile(estimates, [.025, .975], axis=0)
    ed = empirical(simulated, diagnostic_estimator.spec)
    md = np.array([diagnostic_estimator.model(t)[0] for t in estimates])
    weights = diagnostic_estimator.spec.weights
    bs_real = len(y)*np.sum(weights*(ed.real-md)**2, axis=1)
    bs_imag = len(y)*np.sum(weights*ed.imag**2, axis=1)
    obs_d = empirical(y, diagnostic_estimator.spec)[0]
    obs_m = diagnostic_estimator.model(theta)[0]
    obs_real = len(y)*np.sum(weights*(obs_d.real-obs_m)**2)
    obs_imag = len(y)*np.sum(weights*obs_d.imag**2)
    p = (1+np.sum(np.c_[bs_real+bs_imag, bs_real, bs_imag] >=
                  np.array([obs_real+obs_imag, obs_real, obs_imag]), axis=0))/(draws+1)
    return {"fit":fit, "limits":limits, "percentile":percentile, "se":estimates.std(axis=0, ddof=1),
            "p":p, "bootstrap_boundary":np.mean([z["boundary"] for z in fits]),
            "bootstrap_failures":sum(not z["converged"] for z in fits)}
