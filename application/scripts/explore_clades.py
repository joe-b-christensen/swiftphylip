#!/usr/bin/env python3
"""Exploratory: composition of the smallest clade (one side of an edge)
holding at least 95% of each super-population, and population mix of tips outside the largest pure clade."""
import json, sys, collections
import numpy as np
d = sys.argv[1]
tp = json.load(open(d + "/tree_parent.json"))
parent, sup, pops, ntip = tp["parent"], tp["labels"], tp["pops"], tp["ntip"]
nn = len(parent)
root = [i for i, p in enumerate(parent) if p < 0][0]
children = [[] for _ in range(nn)]
for i, p in enumerate(parent):
    if p >= 0: children[p].append(i)
order = [root]; k = 0
while k < len(order):
    order += children[order[k]]; k += 1
G = ["AFR", "AMR", "EAS", "EUR", "SAS"]
cnt = np.zeros((nn, 5), int)
for t in range(ntip): cnt[t, G.index(sup[t])] = 1
for x in reversed(order[1:]): cnt[parent[x]] += cnt[x]
tot = cnt[root]
out = {}
for gi, g in enumerate(G):
    best = None
    for x in order[1:]:
        for side in (cnt[x], tot - cnt[x]):
            if side[gi] >= 0.95 * tot[gi]:
                if best is None or side.sum() < best.sum(): best = side.copy()
    out[g] = {"clade_size": int(best.sum()), "composition": {h: int(best[j]) for j, h in enumerate(G)}}
# populations outside the largest pure clade, per super-population
def pure_best(gi):
    bestsz, bestnode, bestside = 0, None, None
    for x in order[1:]:
        for which, side in ((0, cnt[x]), (1, tot - cnt[x])):
            if side.sum() == side[gi] and side[gi] > bestsz: bestsz, bestnode, bestside = side[gi], x, which
    return bestnode, bestside
desc = {}
def tips_under(x):
    st, res = [x], []
    while st:
        y = st.pop()
        if y < ntip: res.append(y)
        st += children[y]
    return res
for gi, g in enumerate(G):
    node, side = pure_best(gi)
    inside = set(tips_under(node))
    if side == 1: inside = set(range(ntip)) - inside
    outside = [t for t in range(ntip) if sup[t] == g and t not in inside]
    out[g]["outside_largest_pure_clade_by_population"] = dict(collections.Counter(pops[t] for t in outside).most_common())
json.dump(out, open(d + "/exploratory_clades.json", "w"), indent=1)
for g in G: print(g, out[g])
