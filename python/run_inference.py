"""Monte Carlo coverage and size for copula-spectral bootstrap inference."""

import os
for name in ["OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"]:
    os.environ[name] = "1"
import argparse
import csv
import json
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import numpy as np
from cse import Bank, Estimator, features, generate, bootstrap

CASES = [(n, p, v) for n in [512, 1024, 2048]
         for p, v in [(.6, 1.), (.9, 1.), (.98, .25), (.98, 1.)]]
CONTEXT = {}


def initialize(path):
    CONTEXT["bank"] = Bank.load(path)


def run_chunk(task):
    cell, first, last, draws, seed = task
    n, phi, variance = CASES[cell]
    if n not in CONTEXT:
        e = Estimator(features(n), CONTEXT["bank"])
        e.dense_starts()
        d = Estimator(features(n, diagnostic=True), CONTEXT["bank"])
        CONTEXT[n] = (e, d)
    e, d = CONTEXT[n]
    output = []
    for replicate in range(first, last):
        rng = np.random.default_rng(np.random.SeedSequence([seed, cell, replicate]))
        truth = np.array([phi, variance])
        y = generate(truth, n, rng)[0]
        ans = bootstrap(y, e, d, rng, draws)
        row = {"cell":cell, "replicate":replicate, "n":n, "phi":phi, "v":variance,
               "bootstrap_draws":draws, "boundary":int(ans["fit"]["boundary"]),
               "converged":int(ans["fit"]["converged"]),
               "bootstrap_boundary":ans["bootstrap_boundary"],
               "bootstrap_failures":ans["bootstrap_failures"]}
        for k, name in enumerate(["phi", "v"]):
            row[name+"_hat"] = ans["fit"]["theta"][k]
            row[name+"_se"] = ans["se"][k]
            for kind in ["limits", "percentile"]:
                lo, hi = ans[kind][:, k]
                prefix = "basic" if kind == "limits" else kind
                row[prefix+"_"+name+"_coverage"] = int(lo <= truth[k] <= hi)
                row[prefix+"_"+name+"_length"] = max(0., hi-lo)
        for k, name in enumerate(["global", "real", "imaginary"]):
            row[name+"_p"] = ans["p"][k]
            row[name+"_reject"] = int(ans["p"][k] <= .05)
        output.append(row)
    return output


def summarize(path):
    import pandas as pd
    data = pd.read_csv(path)
    rows = []
    for cell, z in data.groupby("cell", sort=True):
        row = {"cell":cell, "n":int(z.n.iloc[0]), "phi":z.phi.iloc[0],
               "v":z.v.iloc[0], "replications":len(z), "bootstrap_draws":int(z.bootstrap_draws.iloc[0])}
        for name in ["phi", "v"]:
            row[name+"_bias"] = np.mean(z[name+"_hat"]-z[name])
            row[name+"_sd"] = z[name+"_hat"].std(ddof=1)
            row[name+"_mean_se"] = z[name+"_se"].mean()
            for kind in ["basic", "percentile"]:
                col = kind+"_"+name+"_coverage"
                row[col] = z[col].mean()
                row[col+"_mcse"] = np.sqrt(row[col]*(1-row[col])/len(z))
                row[kind+"_"+name+"_length"] = z[kind+"_"+name+"_length"].mean()
        for name in ["global", "real", "imaginary"]:
            col = name+"_reject"
            row[col] = z[col].mean()
            row[col+"_mcse"] = np.sqrt(row[col]*(1-row[col])/len(z))
        row["boundary"] = z.boundary.mean()
        row["fit_failures"] = int((1-z.converged).sum())
        row["bootstrap_failures"] = int(z.bootstrap_failures.sum())
        rows.append(row)
    result = pd.DataFrame(rows)
    result.to_csv(path.with_name("inference_summary.csv"), index=False)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--replications", type=int, default=1000)
    parser.add_argument("--bootstrap", type=int, default=199)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--seed", type=int, default=20261101)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[1]/"results")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    settings = {"replications":args.replications,"bootstrap":args.bootstrap,"seed":args.seed,"cases":CASES}
    config = args.output/"inference_settings.json"
    if config.exists() and json.loads(config.read_text()) != json.loads(json.dumps(settings)):
        raise ValueError("The output directory contains a different experiment design")
    config.write_text(json.dumps(settings, indent=2)+"\n")
    bank_path = args.output/"inference_coefficients.npz"
    if not bank_path.exists():
        levels = np.unique(np.concatenate([features(n).tau for n in [512,1024,2048]] +
                                           [features(1024, diagnostic=True).tau]))
        Bank(levels).save(bank_path)
    path = args.output/"inference_replicates.csv"
    done = set()
    if path.exists():
        with path.open() as file:
            done = {(int(r["cell"]),int(r["replicate"])) for r in csv.DictReader(file)}
    tasks = []
    for cell in range(len(CASES)):
        for rep in range(0,args.replications,5):
            last = min(rep+5,args.replications)
            status = [(cell,r) in done for r in range(rep,last)]
            if all(status):
                continue
            if any(status):
                tasks.extend((cell,r,r+1,args.bootstrap,args.seed) for r in range(rep,last) if (cell,r) not in done)
            else:
                tasks.append((cell,rep,last,args.bootstrap,args.seed))
    began = time.time()
    finished = len(done)
    with ProcessPoolExecutor(max_workers=args.workers, initializer=initialize, initargs=(str(bank_path),)) as pool:
        futures = [pool.submit(run_chunk,t) for t in tasks]
        for future in as_completed(futures):
            rows = future.result()
            exists = path.exists()
            with path.open("a",newline="") as file:
                writer = csv.DictWriter(file,fieldnames=list(rows[0]))
                if not exists:
                    writer.writeheader()
                writer.writerows(rows)
            finished += len(rows)
            if finished % 100 == 0:
                print(f"Completed {finished}/{len(CASES)*args.replications}; elapsed {time.time()-began:.0f} s",flush=True)
    print(summarize(path).to_string(index=False),flush=True)


if __name__ == "__main__":
    main()
