#!/usr/bin/env python3
"""SwiftPhylip correctness gate and timer.

A candidate passes only if, for every case, its exit status, its stdout and every file and
directory it leaves in the working directory are byte-identical to those of the live reference
(vcf2phylip v2.9, run fresh in this evaluation) on:
  * fresh random fuzz cases from gen_data.gen_case (new seed from os.urandom unless --seed);
  * the fixed adversarial files in data/adversarial/ (two option sets each);
  * every real development set, under three option sets.

Speed: single process runs, each timed with perf_counter around the whole process; outputs are
deleted before each repetition, outside the timer; repetitions continue in a fixed window
(candidate: >= 5 runs and >= 3 s, at most 40; reference: >= 1 run and >= 10 s, at most 3).
The primary number is the geometric mean, over the real development sets and the two timed
option sets (default PHYLIP; all four formats -f -n -b), of median(reference)/median(candidate).

Commands
  evaluate.py gate   --cand CMD [--fuzz N] [--seed S] [--threads T]
  evaluate.py full   --cand CMD [--fuzz N] [--seed S] [--threads T]   (gate + timing on dev sets)
  evaluate.py final  --cand CMD                                       (dev + held-out, 1 and 2 threads)
CMD is a command prefix, e.g. "solver/build/swiftphylip"; relative paths resolve against this folder.
The candidate gets SWIFTPHYLIP_THREADS=T in its environment.
"""
import argparse
import hashlib
import json
import math
import os
import random
import shlex
import shutil
import statistics
import subprocess
import sys
import tempfile
import time

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)
REF_CMD = [sys.executable, "-I", os.path.join(ROOT, "reference", "vcf2phylip.py")]

DEV = ["data/dev/kg_all2504_r100k.vcf", "data/dev/kg_s250_r300k.vcf.gz", "data/dev/pinf_sc50_first_half.vcf.gz"]
HELD = ["data/heldout/kg_all2504_r700k.vcf", "data/heldout/kg_s500_r900k.vcf.gz", "data/heldout/pinf_sc50_second_half.vcf"]
ADV = ["data/adversarial/adv_progress.vcf", "data/adversarial/adv_wide_crlf.vcf.gz", "data/adversarial/adv_medium.vcf"]
TIMED = {"default": [], "all": ["-f", "-n", "-b"]}


def extra_config(path):
    return ["-p", "-b", "-w", "-f", "-o", "__THIRD__", "--output-prefix", "q", "--output-folder", "o/x", "-m", "17"]


def sha(b):
    return hashlib.sha256(b).hexdigest()


def file_sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        while True:
            b = f.read(1 << 22)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def env_for(threads):
    e = {"PATH": "/usr/local/bin:/usr/bin:/bin", "LC_ALL": "C.UTF-8", "LANG": "C.UTF-8", "PYTHONUTF8": "1",
         "PYTHONIOENCODING": "utf-8", "PYTHONHASHSEED": "0", "HOME": "/tmp", "TZ": "UTC"}
    if threads is not None:
        e["SWIFTPHYLIP_THREADS"] = str(threads)
    return e


def resolve(cmd):
    parts = shlex.split(cmd)
    out = []
    for i, p in enumerate(parts):
        q = os.path.join(ROOT, p)
        if not os.path.isabs(p) and (os.sep in p) and os.path.exists(q):
            out.append(q)
        else:
            out.append(p)
    return out


WORK = None


def new_dir():
    return tempfile.mkdtemp(dir=WORK)


def stage(d, fname, data=None, path=None):
    os.mkdir(os.path.join(d, "in"))
    dst = os.path.join(d, "in", fname)
    if data is not None:
        with open(dst, "wb") as f:
            f.write(data)
    else:
        os.symlink(os.path.join(ROOT, path), dst)
    return "in/" + fname


def clean(d):
    for n in os.listdir(d):
        if n == "in":
            continue
        p = os.path.join(d, n)
        if os.path.isdir(p) and not os.path.islink(p):
            shutil.rmtree(p)
        else:
            os.unlink(p)
    ind = os.path.join(d, "in")
    for n in os.listdir(ind):
        p = os.path.join(ind, n)
        if not os.path.islink(p) and n not in KEEP_IN:
            if os.path.isdir(p):
                shutil.rmtree(p)
            else:
                os.unlink(p)


KEEP_IN = set()
STOP_ON_FAIL = False


def collect(d, input_rel, keep_bytes):
    out = {}
    for base, dirs, files in os.walk(d):
        rel = os.path.relpath(base, d)
        if rel != ".":
            out[rel + "/"] = "DIR"
        for n in files:
            r = os.path.normpath(os.path.join(rel, n))
            if r == input_rel:
                continue
            p = os.path.join(base, n)
            if keep_bytes:
                with open(p, "rb") as f:
                    out[r] = f.read()
            else:
                out[r] = file_sha(p)
    return out


def run_once(cmd, d, input_rel, args, threads, timeout=3600):
    t0 = time.perf_counter()
    p = subprocess.run(cmd + ["-i", input_rel] + args, cwd=d, env=env_for(threads),
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=timeout)
    t = time.perf_counter() - t0
    return p, t


def result(p, d, input_rel, keep_bytes):
    return {"exit": p.returncode, "stdout": p.stdout if keep_bytes else sha(p.stdout),
            "files": collect(d, input_rel, keep_bytes)}


def diff(a, b):
    if a["exit"] != b["exit"]:
        return "exit status %r vs %r" % (a["exit"], b["exit"])
    if a["stdout"] != b["stdout"]:
        if isinstance(a["stdout"], bytes):
            x, y = a["stdout"], b["stdout"]
            i = next((k for k in range(min(len(x), len(y))) if x[k] != y[k]), min(len(x), len(y)))
            return "stdout differs at byte %d: ref %r vs cand %r" % (i, x[max(0, i - 30):i + 30], y[max(0, i - 30):i + 30])
        return "stdout differs"
    ka, kb = set(a["files"]), set(b["files"])
    if ka != kb:
        return "file sets differ: ref-only %s cand-only %s" % (sorted(ka - kb), sorted(kb - ka))
    for k in sorted(ka):
        if a["files"][k] != b["files"][k]:
            x, y = a["files"][k], b["files"][k]
            if isinstance(x, bytes):
                i = next((j for j in range(min(len(x), len(y))) if x[j] != y[j]), min(len(x), len(y)))
                return "file %s differs at byte %d (len %d vs %d): ref %r vs cand %r" % (
                    k, i, len(x), len(y), x[max(0, i - 20):i + 20], y[max(0, i - 20):i + 20])
            return "file %s differs" % k
    return None


def compare_case(ref, cand, fname, data, path, args, threads, keep_bytes):
    res = []
    for cmd, th in ((ref, None), (cand, threads)):
        d = new_dir()
        inp = stage(d, fname, data=data, path=path)
        p, _ = run_once(cmd, d, inp, args, th)
        res.append(result(p, d, inp, keep_bytes))
        shutil.rmtree(d)
    return diff(res[0], res[1])


def third_sample(path):
    import gzip
    op = gzip.open if path.lower().endswith(".gz") else open
    with op(os.path.join(ROOT, path), "rt") as f:
        for line in f:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[11].replace("./", "")
    return "x"


def args_for(path, cfg):
    if cfg == "extra":
        a = extra_config(path)
        return [third_sample(path) if x == "__THIRD__" else x for x in a]
    return list(TIMED[cfg])


def gate(cand, nfuzz, seed, threads, log):
    ref = REF_CMD
    fails = []
    rng = random.Random(seed)
    import gen_data
    t0 = time.time()
    for i in range(nfuzz):
        fname, data, args = gen_data.gen_case(rng)
        r = compare_case(ref, cand, fname, data, None, args, threads, True)
        if r:
            fails.append({"case": i, "fname": fname, "args": args, "why": r})
            if len(fails) <= 3:
                log("FUZZ FAIL case %d %s %s: %s" % (i, fname, args, r))
                fd = os.path.join(ROOT, "results", "fails")
                os.makedirs(fd, exist_ok=True)
                with open(os.path.join(fd, "seed%d_case%d_%s" % (seed, i, fname)), "wb") as f:
                    f.write(data)
                with open(os.path.join(fd, "seed%d_case%d_%s.args.json" % (seed, i, fname)), "w") as f:
                    json.dump(args, f)
            if STOP_ON_FAIL:
                return fails, len(fails)
    log("fuzz: %d cases, seed %d, %d failures, %.1f s" % (nfuzz, seed, len(fails), time.time() - t0))
    nf = len(fails)
    for path in ADV:
        for args in ([], ["-f", "-n", "-b", "-w", "-o", "S2", "-m", "2"]):
            r = compare_case(ref, cand, os.path.basename(path), None, path, args, threads, False)
            log("adversarial %s %s: %s" % (os.path.basename(path), args, "OK" if not r else "FAIL " + r))
            if r:
                fails.append({"adv": path, "args": args, "why": r})
                if STOP_ON_FAIL:
                    return fails, nf
    for path in DEV:
        args = args_for(path, "extra")
        r = compare_case(ref, cand, os.path.basename(path), None, path, args, threads, False)
        log("dev identity %s %s: %s" % (os.path.basename(path), args, "OK" if not r else "FAIL " + r))
        if r:
            fails.append({"dev": path, "args": args, "why": r})
    return fails, nf


def timed(cmd, path, args, threads, min_runs, window, max_runs, check_against=None):
    d = new_dir()
    inp = stage(d, os.path.basename(path), path=path)
    times = []
    first = None
    start = time.perf_counter()
    while True:
        clean(d)
        p, t = run_once(cmd, d, inp, args, threads)
        times.append(t)
        if first is None:
            first = result(p, d, inp, False)
        if len(times) >= max_runs or (len(times) >= min_runs and time.perf_counter() - start >= window):
            break
    shutil.rmtree(d)
    return times, first


def timing(cand, sets, threads, log, ref_cache=None, ref_runs=(1, 10.0, 3), cand_runs=(5, 3.0, 40)):
    rows = []
    fails = []
    for path in sets:
        for cfg in TIMED:
            args = args_for(path, cfg)
            key = (path, cfg)
            if ref_cache is not None and key in ref_cache:
                rt, rres = ref_cache[key]
            else:
                rt, rres = timed(REF_CMD, path, args, None, *ref_runs)
                if ref_cache is not None:
                    ref_cache[key] = (rt, rres)
            ct, cres = timed(cand, path, args, threads, *cand_runs)
            why = diff(rres, cres)
            if why:
                fails.append({"timed": path, "cfg": cfg, "why": why})
            row = {"set": os.path.basename(path), "cfg": cfg, "threads": threads,
                   "ref_median": statistics.median(rt), "ref_best3": min(rt[:3]), "ref_n": len(rt),
                   "cand_median": statistics.median(ct), "cand_best3": min(ct[:3]), "cand_n": len(ct),
                   "identical": not why}
            row["speedup_median"] = row["ref_median"] / row["cand_median"]
            row["speedup_best3"] = row["ref_best3"] / row["cand_best3"]
            rows.append(row)
            log("time %-34s %-7s thr=%s ref %.3fs (n=%d) cand %.4fs (n=%d) x%.1f %s" % (
                row["set"], cfg, threads, row["ref_median"], row["ref_n"], row["cand_median"], row["cand_n"],
                row["speedup_median"], "OK" if not why else "MISMATCH " + why))
    return rows, fails


def geomean(xs):
    return math.exp(sum(math.log(x) for x in xs) / len(xs)) if xs else float("nan")


def main():
    global WORK, STOP_ON_FAIL
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["gate", "full", "final"])
    ap.add_argument("--cand")
    ap.add_argument("--fuzz", type=int, default=400)
    ap.add_argument("--seed", type=int)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--label", default="")
    ap.add_argument("--stop-on-fail", action="store_true", help="end the gate at the first failure")
    ap.add_argument("--out", default=os.path.join(ROOT, "results"))
    a = ap.parse_args()
    STOP_ON_FAIL = a.stop_on_fail
    os.makedirs(a.out, exist_ok=True)
    WORK = tempfile.mkdtemp(prefix="swiftgate_")
    stamp = time.strftime("%Y%m%d-%H%M%S")
    logf = open(os.path.join(a.out, "eval_%s_%s%s.log" % (a.mode, stamp, ("_" + a.label) if a.label else "")), "w")

    def log(s):
        print(s, flush=True)
        logf.write(s + "\n")
        logf.flush()

    log("evaluate.py %s cand=%r threads=%s" % (a.mode, a.cand, a.threads))
    summary = {"mode": a.mode, "cand": a.cand, "threads": a.threads, "label": a.label}
    ok = True
    if a.mode in ("gate", "full"):
        cand = resolve(a.cand)
        seed = a.seed if a.seed is not None else int.from_bytes(os.urandom(6), "big")
        fails, nfuzzfail = gate(cand, a.fuzz, seed, a.threads, log)
        summary.update({"seed": seed, "fuzz": a.fuzz, "fuzz_failures": nfuzzfail, "failures": fails})
        if a.mode == "full":
            rows, tf = timing(cand, DEV, a.threads, log)
            fails += tf
            summary["timing"] = rows
            summary["geomean_speedup_median"] = geomean([r["speedup_median"] for r in rows])
            summary["geomean_speedup_best3"] = geomean([r["speedup_best3"] for r in rows])
            log("PRIMARY geomean speed-up (median) over dev sets x {default, all}: %.2f" % summary["geomean_speedup_median"])
        ok = not fails
        log("GATE %s (%d failures)" % ("PASS" if ok else "FAIL", len(fails)))
    elif a.mode == "final":
        cand = resolve(a.cand)
        cache = {}
        allrows = []
        fails = []
        for th in (1, 2):
            for name, sets in (("dev", DEV), ("heldout", HELD)):
                rows, tf = timing(cand, sets, th, log, ref_cache=cache, ref_runs=(3, 0.0, 3), cand_runs=(5, 3.0, 40))
                for r in rows:
                    r["group"] = name
                allrows += rows
                fails += tf
        summary["final"] = allrows
        for th in (1, 2):
            for name in ("dev", "heldout"):
                rs = [r for r in allrows if r["threads"] == th and r["group"] == name]
                summary["geomean_%s_t%d_median" % (name, th)] = geomean([r["speedup_median"] for r in rs])
                summary["geomean_%s_t%d_best3" % (name, th)] = geomean([r["speedup_best3"] for r in rs])
                log("FINAL %s threads=%d geomean speed-up: median %.2f, best-of-three %.2f" % (
                    name, th, summary["geomean_%s_t%d_median" % (name, th)], summary["geomean_%s_t%d_best3" % (name, th)]))
        ok = not fails
        summary["failures"] = fails
        log("FINAL identity %s" % ("PASS" if ok else "FAIL"))
    with open(os.path.join(a.out, "eval_%s_%s%s.json" % (a.mode, stamp, ("_" + a.label) if a.label else "")), "w") as f:
        json.dump(summary, f, indent=1, default=lambda o: o.hex() if isinstance(o, bytes) else str(o))
    shutil.rmtree(WORK, ignore_errors=True)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
