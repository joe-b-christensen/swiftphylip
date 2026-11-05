#!/usr/bin/env python3
"""Reference scaling on chromosome 22 (all 2504 samples, PHYLIP + binary NEXUS with -b, as in the application) for the
extrapolation to the full chromosome, plus byte-identity of SwiftPhylip on the same slices.
Slices: records 200000 .. 200000+n of the chr22 file (not used anywhere else).
Usage: scaling.py OUT_JSON SWIFT_BINARY CHR22_SWIFT_SECONDS"""
import gzip, hashlib, json, os, subprocess, sys, tempfile, time, shutil

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CHR22 = os.path.join(ROOT, "data/raw/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz")
ENV = {"PATH": "/usr/bin:/bin", "LC_ALL": "C.UTF-8", "PYTHONUTF8": "1", "PYTHONIOENCODING": "utf-8", "PYTHONHASHSEED": "0"}
SIZES = [1000, 4000, 10000, 20000]
work = tempfile.mkdtemp(prefix="scal_")
start = 200000
p = subprocess.Popen(["zcat", CHR22], stdout=subprocess.PIPE)
head, recs = [], []
for line in p.stdout:
    if line.startswith(b"#"):
        head.append(line); continue
    recs.append(line)
    if len(recs) >= start + max(SIZES): break
p.kill()
recs = recs[start:]
pts = []
for n in SIZES:
    f = os.path.join(work, "slice%d.vcf" % n)
    with open(f, "wb") as o:
        o.writelines(head); o.writelines(recs[:n])
    res = {}
    for name, cmd in (("ref", [sys.executable, "-I", os.path.join(ROOT, "reference/vcf2phylip.py")]), ("swift", [sys.argv[2]])):
        d = os.path.join(work, name + str(n)); os.mkdir(d)
        t = time.perf_counter()
        r = subprocess.run(cmd + ["-i", f, "-b", "--output-folder", d], env=ENV, stdout=subprocess.PIPE)
        el = time.perf_counter() - t
        out = os.path.join(d, "slice%d.min4.phy" % n)
        outb = os.path.join(d, "slice%d.min4.bin.nexus" % n)
        res[name] = (el, hashlib.sha256(open(out, "rb").read() + open(outb, "rb").read()).hexdigest(), hashlib.sha256(r.stdout.replace(d.encode(), b"D")).hexdigest(), r.returncode)
    pts.append({"records": n, "seconds": res["ref"][0], "swift_seconds": res["swift"][0], "identical_phy": res["ref"][1] == res["swift"][1],
                "identical_stdout": res["ref"][2] == res["swift"][2], "exit": [res["ref"][3], res["swift"][3]]})
    print(pts[-1], flush=True)
xs = [q["records"] for q in pts]; ys = [q["seconds"] for q in pts]
mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
slope = sum((a - mx) * (b - my) for a, b in zip(xs, ys)) / sum((a - mx) ** 2 for a in xs)
icpt = my - slope * mx
out = {"reference_points": pts, "fit_slope": slope, "fit_intercept": icpt, "chr22_records": 1103547,
       "chr22_swift_seconds": float(sys.argv[3]), "extrapolated_reference_seconds_chr22": slope * 1103547 + icpt,
       "note": "PHYLIP + binary NEXUS (-b, as in the application), 2504 samples; the full-chromosome reference time is an extrapolation of a least-squares line through the measured points"}
json.dump(out, open(sys.argv[1], "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != "reference_points"}, indent=1))
shutil.rmtree(work)
