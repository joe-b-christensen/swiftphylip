#!/usr/bin/env python3
"""Neighbour-joining tree from an IBS distance matrix and agreement with 1000 Genomes labels.

Usage: nj_analysis.py DIST_PREFIX PANEL OUTDIR
Reads DIST_PREFIX.names.txt and DIST_PREFIX.dist.f64 (from ibs_dist), the 1000 Genomes panel
file, and writes OUTDIR/tree.nwk, OUTDIR/application_results.json and the figure.
Primary result: pure-clade coverage per super-population.
"""
import json
import math
import sys
import time

import numpy as np


def nj(D):
    """Classic Saitou-Nei NJ. Returns edges [(u, v, length)], tips 0..n-1, internal nodes n.."""
    n0 = D.shape[0]
    D = D.astype(np.float64).copy()
    ids = list(range(n0))
    r = D.sum(axis=1)
    edges = []
    nxt = n0
    n = n0
    while n > 2:
        r = D.sum(axis=1)
        Q = (n - 2) * D - r[:, None] - r[None, :]
        np.fill_diagonal(Q, np.inf)
        k = int(np.argmin(Q))
        i, j = divmod(k, n)
        if i > j:
            i, j = j, i
        dij = D[i, j]
        li = 0.5 * dij + (r[i] - r[j]) / (2.0 * (n - 2))
        lj = dij - li
        u = nxt
        nxt += 1
        edges.append((ids[i], u, li))
        edges.append((ids[j], u, lj))
        du = 0.5 * (D[i] + D[j] - dij)
        # new node u takes slot i; slot j removed by moving the last slot into it
        D[i, :] = du
        D[:, i] = du
        D[i, i] = 0.0
        ids[i] = u
        last = n - 1
        if j != last:
            D[j, :] = D[last, :]
            D[:, j] = D[:, last]
            D[j, j] = 0.0
            ids[j] = ids[last]
        D = np.ascontiguousarray(D[:last, :last])
        ids = ids[:last]
        n -= 1
    edges.append((ids[0], ids[1], D[0, 1]))
    return edges, nxt


def adjacency(edges, nnodes):
    adj = [[] for _ in range(nnodes)]
    for u, v, l in edges:
        adj[u].append((v, l))
        adj[v].append((u, l))
    return adj


def rooted_order(adj, root):
    parent = [-1] * len(adj)
    plen = [0.0] * len(adj)
    order = [root]
    seen = [False] * len(adj)
    seen[root] = True
    k = 0
    while k < len(order):
        x = order[k]
        k += 1
        for y, l in adj[x]:
            if not seen[y]:
                seen[y] = True
                parent[y] = x
                plen[y] = l
                order.append(y)
    return parent, plen, order


def newick(adj, ntip, names, root):
    parent, plen, order = rooted_order(adj, root)
    children = [[] for _ in range(len(adj))]
    for x in order[1:]:
        children[parent[x]].append(x)
    s = {}
    for x in reversed(order):
        if x < ntip:
            s[x] = names[x]
        else:
            s[x] = "(" + ",".join("%s:%.8g" % (s[c], plen[c]) for c in children[x]) + ")"
    return s[root] + ";"


def pure_clade_coverage(adj, ntip, labels, groups, root):
    """Largest one-side-of-an-edge tip set whose tips all carry label g, divided by |g|."""
    parent, plen, order = rooted_order(adj, root)
    G = len(groups)
    gi = {g: k for k, g in enumerate(groups)}
    cnt = np.zeros((len(adj), G), dtype=np.int64)
    for t in range(ntip):
        cnt[t, gi[labels[t]]] = 1
    for x in reversed(order[1:]):
        cnt[parent[x]] += cnt[x]
    total = cnt[root].copy()
    best = np.zeros(G, dtype=np.int64)
    for x in order[1:]:
        for side in (cnt[x], total - cnt[x]):
            sz = side.sum()
            nz = np.nonzero(side)[0]
            if len(nz) == 1 and sz > best[nz[0]]:
                best[nz[0]] = sz
    return {g: (int(best[gi[g]]), int(total[gi[g]]), float(best[gi[g]] / total[gi[g]])) for g in groups}


def patristic_nn(adj, ntip):
    """For each tip, the nearest other tip by patristic distance (lowest index on ties)."""
    nn = np.zeros(ntip, dtype=np.int64)
    N = len(adj)
    for t in range(ntip):
        dist = np.full(N, np.inf)
        dist[t] = 0.0
        stack = [t]
        while stack:
            x = stack.pop()
            for y, l in adj[x]:
                if dist[y] == np.inf:
                    dist[y] = dist[x] + l
                    stack.append(y)
        d = dist[:ntip].copy()
        d[t] = np.inf
        nn[t] = int(np.argmin(d))
    return nn


def equal_angle(adj, ntip, root):
    parent, plen, order = rooted_order(adj, root)
    children = [[] for _ in range(len(adj))]
    for x in order[1:]:
        children[parent[x]].append(x)
    ntips = [0] * len(adj)
    for x in reversed(order):
        ntips[x] = 1 if x < ntip else sum(ntips[c] for c in children[x])
    pos = np.zeros((len(adj), 2))
    wedge = {root: (0.0, 2 * math.pi)}
    for x in order:
        a0, a1 = wedge[x]
        a = a0
        for c in children[x]:
            w = (a1 - a0) * ntips[c] / ntips[x]
            wedge[c] = (a, a + w)
            mid = a + w / 2
            L = max(plen[c], 0.0)
            pos[c] = pos[x] + L * np.array([math.cos(mid), math.sin(mid)])
            a += w
    return pos, parent


def main():
    prefix, panel, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
    t0 = time.time()
    names = [l.rstrip("\n") for l in open(prefix + ".names.txt")]
    n = len(names)
    D = np.fromfile(prefix + ".dist.f64", dtype=np.float64).reshape(n, n)
    lab, pop = {}, {}
    for line in open(panel):
        f = line.rstrip("\n").split("\t")
        if f[0] == "sample":
            continue
        pop[f[0]] = f[1]
        lab[f[0]] = f[2]
    sup = [lab[s] for s in names]
    pops = [pop[s] for s in names]
    t1 = time.time()
    edges, nnodes = nj(D)
    t2 = time.time()
    adj = adjacency(edges, nnodes)
    root = nnodes - 1
    groups = ["AFR", "AMR", "EAS", "EUR", "SAS"]
    cov = pure_clade_coverage(adj, n, sup, groups, root)
    weighted = sum(cov[g][0] for g in groups) / n
    popnames = sorted(set(pops))
    pcov = pure_clade_coverage(adj, n, pops, popnames, root)
    nn = patristic_nn(adj, n)
    nn_ok = np.array([sup[i] == sup[nn[i]] for i in range(n)])
    nn_by = {g: float(np.mean([nn_ok[i] for i in range(n) if sup[i] == g])) for g in groups}
    pop_nn = float(np.mean([pops[i] == pops[nn[i]] for i in range(n)]))
    t3 = time.time()
    neg = sum(1 for e in edges if e[2] < 0)
    res = {
        "n_tips": n,
        "n_edges": len(edges),
        "negative_branch_lengths": neg,
        "primary_pure_clade_coverage": {g: {"largest_pure_clade": cov[g][0], "n": cov[g][1], "coverage": cov[g][2]} for g in groups},
        "primary_weighted_mean_coverage": weighted,
        "nn_agreement_overall": float(nn_ok.mean()),
        "nn_agreement_by_superpop": nn_by,
        "population_pure_clade_coverage": {p: {"largest_pure_clade": pcov[p][0], "n": pcov[p][1], "coverage": pcov[p][2]} for p in popnames},
        "population_nn_agreement": pop_nn,
        "seconds_load": t1 - t0, "seconds_nj": t2 - t1, "seconds_metrics": t3 - t2,
    }
    with open(outdir + "/application_results.json", "w") as f:
        json.dump(res, f, indent=1)
    with open(outdir + "/tree.nwk", "w") as f:
        f.write(newick(adj, n, names, root) + "\n")
    pos, parent = equal_angle(adj, n, root)
    np.save(outdir + "/tree_layout.npy", pos)
    with open(outdir + "/tree_parent.json", "w") as f:
        json.dump({"parent": parent, "labels": sup, "pops": pops, "names": names, "ntip": n}, f)
    print(json.dumps({k: v for k, v in res.items() if k != "population_pure_clade_coverage"}, indent=1))


if __name__ == "__main__":
    main()
