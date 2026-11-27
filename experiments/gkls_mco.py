"""Парное сравнение двух режимов MCOMethod на ограниченных задачах GKLS 2D."""
import argparse
import csv
import os
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np


def frontier(points):
    if not len(points):
        return np.empty((0, 2))
    points = points[np.lexsort((points[:, 1], points[:, 0]))]
    keep = np.r_[True, points[1:, 1] < np.minimum.accumulate(points[:-1, 1])]
    return points[keep]


def hv(points):
    # После нормировки одна опорная точка (1.1, 1.1) для обоих методов.
    front = frontier(points[np.all(points < 1.1, axis=1)])
    if not len(front):
        return 0.0
    return np.dot(1.1 - front[:, 0], -np.diff(np.r_[1.1, front[:, 1]]))


def igd_plus(reference, points):
    if not len(points):
        return float("inf")
    distance = 0.0
    for start in range(0, len(reference), 64):
        delta = np.maximum(points[None, :, :] - reference[start:start + 64, None, :], 0)
        distance += np.sqrt(np.min(np.sum(delta * delta, axis=2), axis=1)).sum()
    return distance / len(reference)


def first_hit(trials, target):
    def reached(n):
        prefix = trials[:n]
        return hv(prefix[np.isfinite(prefix[:, 0])]) >= target

    if not reached(len(trials)):
        return None
    left, right = 1, len(trials)
    while left < right:
        mid = (left + right) // 2
        if reached(mid):
            right = mid
        else:
            left = mid + 1
    return left


def evaluate(exe, library, out, instance):
    folder = out / f"instance-{instance:03d}"
    folder.mkdir()

    def run(name, *args):
        path = folder / (name + ".csv")
        proc = subprocess.run(
            [str(exe), str(library), str(instance), str(path), *map(str, args)],
            cwd=folder, env={**os.environ, "OMP_NUM_THREADS": "1"},
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        )
        log = proc.stdout.decode("utf-8", errors="replace")
        (folder / (name + ".log")).write_text(log, encoding="utf-8")
        if proc.returncode:
            raise RuntimeError(f"{name} on instance {instance}: {log[-1200:]}")
        return path, log

    def grid(size):
        path, log = run(f"grid-{size}", "grid", size)
        values = np.loadtxt(path, delimiter=",", skiprows=1, usecols=(0, 1), ndmin=2)
        feasible = int(re.search(r"feasible=(\d+)", log)[1])
        return values, feasible / size**2

    if instance in (1, 50, 100):
        run("map", "map", 151)
    coarse, _ = grid(501)
    fine, fraction = grid(1001)
    # Идеал (-1,-1) известен для скалярных функций. Масштаб замораживаем до запусков.
    scale = np.max(fine, axis=0) + 1
    if np.any(scale <= 0):
        raise ValueError(f"Degenerate reference on instance {instance}")

    def normalize(points):
        return (points + 1) / scale

    runs = []
    for budget in (1000, 3000, 10000):
        for local in (0.0, 1.5):
            path, log = run(f"solve-{budget}-{local:g}", "solve", budget, local)
            data = np.genfromtxt(path, delimiter=",", skip_header=1, ndmin=2)
            trials = normalize(data[:, 3:5])
            valid = np.isfinite(trials).all(axis=1)
            calls = tuple(map(int, re.search(r"calls=([\d,]+)", log)[1].split(",")))
            if (len(trials) > budget or calls[0] != len(trials) or
                    calls[2:] != (sum(valid), sum(valid)) or data[-1, 0] != 49):
                raise ValueError(f"Trial accounting failed: {path}, {calls}, {len(trials)}")
            front = frontier(trials[valid])
            runs.append((budget, local, trials, front, calls, hv(front)))

    coarse_hv, fine_hv = hv(normalize(coarse)), hv(normalize(fine))
    smallest_effect = min(abs(a[-1] - b[-1]) for a, b in zip(runs[::2], runs[1::2]))
    size = 1001
    if abs(fine_hv - coarse_hv) >= max(smallest_effect, 1e-10):
        fine, fraction = grid(2001)
        size = 2001
    reference = normalize(fine)
    reference_hv = hv(reference)
    previous_hv = coarse_hv if size == 1001 else fine_hv
    rows = []
    for budget, local, trials, front, calls, volume in runs:
        rows.append(dict(instance=instance, budget=budget, r_local=local, trials=len(trials),
                         feasible=calls[2], g1_calls=calls[0], g2_calls=calls[1],
                         f1_calls=calls[2], f2_calls=calls[3], pareto=len(front),
                         grid=size, feasible_fraction=fraction, reference_hv=reference_hv,
                         grid_hv_change=abs(reference_hv - previous_hv),
                         scale1=scale[0], scale2=scale[1], hv=volume,
                         relative_hv=volume/reference_hv, igd_plus=igd_plus(reference, front),
                         trials_to_95=first_hit(trials, 0.95 * reference_hv)))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("exe", type=Path)
    parser.add_argument("library", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--instances", type=int, nargs="+", default=list(range(1, 101)))
    parser.add_argument("--jobs", type=int, default=1)
    args = parser.parse_args()
    if len(set(args.instances)) != len(args.instances) or not all(1 <= i <= 100 for i in args.instances):
        parser.error("Expected distinct instance IDs in [1,100]")
    exe, library = args.exe.resolve(strict=True), args.library.resolve(strict=True)
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    with ThreadPoolExecutor(max_workers=args.jobs) as pool, (out / "metrics.csv").open("w", newline="") as f:
        writer = None
        for rows in pool.map(lambda i: evaluate(exe, library, out, i), args.instances):
            if writer is None:
                writer = csv.DictWriter(f, fieldnames=rows[0])
                writer.writeheader()
            writer.writerows(rows)
            f.flush()
            print(f"Instance {rows[0]['instance']} complete", flush=True)


if __name__ == "__main__":
    main()
