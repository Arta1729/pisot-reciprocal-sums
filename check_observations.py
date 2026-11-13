#!/usr/bin/env python3
"""
check_observations.py -- the remaining numerical statements of the paper that are not tables.

  (a) Observation (iii): max |E_n| over 3<=n<=59 (3<=n<=119 for k=30) for the eleven recurrences of Table 1.
  (b) Observation (iv): Tribonacci, max deviation of E_n/(2|kappa|mu^(n-1)) from cos((n-1)theta+arg kappa), 3<=n<=45.
  (c) Section 6: (1/2 - x_k) 4^k for k-generalized Fibonacci numbers.
  (d) Remark on growth in k: (1-mu) k^3 for k = 10, 20, 40, 80, 160, and N_0^root (1-mu) for k = 8, 30.
  (e) Remark on the hypothesis of the density theorem: for F^(k) (3<=k<=40), P^(k) (3<=k<=15), (3;1,0,1) and
      (4;1,0,0,1) the roots of maximal modulus form exactly one non-real conjugate pair, separated from the
      next modulus by a gap (the smallest gap is printed).
Requires mpmath and numpy.
"""
import numpy as np
from mpmath import mp, mpf, mpc, polyroots, findroot, nstr, cos, arg, pi, re

def seq_ints(c, N):
    k = len(c); u = {i: 0 for i in range(2 - k, 1)}; u[1] = 1
    for n in range(2, N + 1):
        u[n] = sum(c[j] * u[n - 1 - j] for j in range(k))
    return u

def E_n(c, n, u, M):
    S = sum(mpf(1) / u[m] for m in range(n, M + 1))
    return 1 / S - (u[n] - u[n - 1])

def prepare(c, nmax, dps):
    mp.dps = dps
    k = len(c); M = nmax + 50
    u = seq_ints(c, M)
    while u[M] <= 10 ** (dps + 25):
        M *= 2; u = seq_ints(c, M)
    return u, M

def binet(c):
    k = len(c)
    rts = sorted(polyroots([1] + [-x for x in c], maxsteps=5000, extraprec=4 * mp.prec), key=lambda z: -abs(z))
    al = rts[0].real
    ds = []
    for a in rts:
        dp = mpf(1)
        for b in rts:
            if b is not a: dp *= (a - b)
        ds.append(a ** (k - 2) / dp)
    return al, rts, ds

# (a)
print("== (a) max |E_n| over 3<=n<=59 (119 for k=30) ==")
rows = [("Fibonacci", [1, 1]), ("Pell", [2, 1]), ("(2;3,1)", [3, 1]), ("Tribonacci", [1, 1, 1]), ("(3;6,1,1)", [6, 1, 1]),
        ("8-bonacci", [1] * 8), ("30-bonacci", [1] * 30), ("(4;3,2,2,1)", [3, 2, 2, 1]),
        ("(5;2,2,1,1,1)", [2, 2, 1, 1, 1]), ("Narayana", [1, 0, 1]), ("(4;1,0,0,1)", [1, 0, 0, 1])]
for name, c in rows:
    hi = 119 if len(c) == 30 else 59
    u, M = prepare(c, hi, 400 if len(c) == 30 else 150)
    print(" %-14s %s" % (name, nstr(max(abs(E_n(c, n, u, M)) for n in range(3, hi + 1)), 4)))

# (b)
print("== (b) Tribonacci: deviation from the leading cosine, 3<=n<=45 ==")
c = [1, 1, 1]; u, M = prepare(c, 45, 150)
al, rts, ds = binet(c)
i = [j for j in range(1, 3) if rts[j].imag > 0][0]
cx = rts[i]; kap = ds[i] * (al - cx) ** 2 / (al ** 2 - cx); mu = abs(cx); th = arg(cx)
dev = max(abs(E_n(c, n, u, M) / (2 * abs(kap) * mu ** (n - 1)) - cos((n - 1) * th + arg(kap))) for n in range(3, 46))
print(" max deviation =", nstr(dev, 3))

# (c)
print("== (c) (1/2 - x_k) 4^k ==")
for k in (2, 3, 4, 6, 12, 30, 40):
    c = [1] * k; u, M = prepare(c, 60, 300)
    x = E_n(c, 2, u, M)
    print(" k=%2d  %s" % (k, nstr((mpf(1) / 2 - x) * 4 ** k, 5)))

# (d), (e): largest non-dominant roots via numpy, refined in mpmath
def top_roots(c, how_many=4):
    mp.dps = 50
    k = len(c); co = [1] + [-x for x in c]
    approx = sorted(np.roots(co), key=lambda z: -abs(z))[1:1 + how_many]
    f = lambda z: sum(co[i] * z ** (k - i) for i in range(k + 1))
    return [findroot(f, mpc(z.real, z.imag), tol=1e-40, maxsteps=100) for z in approx]

print("== (d) (1-mu) k^3 for F^(k) ==")
print(" ", ", ".join("k=%d: %s" % (k, nstr((1 - abs(top_roots([1] * k)[0])) * k ** 3, 4)) for k in (10, 20, 40, 80, 160)))
import verify_pisot as V
for k in (8, 30):
    c = tuple([1] * k)
    P = V.constants(c); N0 = V.N0_from(P, P['C0'], P['C1'], P['q'])
    print(" k=%d  N0_root=%d  N0_root*(1-mu)=%s" % (k, N0, nstr(N0 * (1 - P['mu']), 3)))

print("== (e) single conjugate pair of maximal modulus ==")
bad = []; mingap = None
cases = [("F%d" % k, [1] * k) for k in range(3, 41)] + [("P%d" % k, [2] + [1] * (k - 1)) for k in range(3, 16)] \
        + [("(3;1,0,1)", [1, 0, 1]), ("(4;1,0,0,1)", [1, 0, 0, 1])]
for name, c in cases:
    t = top_roots(c)
    m = [abs(z) for z in t]
    gap = (m[0] - m[2]) if len(m) > 2 else mpf(1)
    pair = abs(t[0].imag) > 1e-15 and abs(t[0] - t[1].conjugate()) < 1e-30 and gap > 1e-12
    if not pair: bad.append(name)
    if mingap is None or gap < mingap[0]: mingap = (gap, name)
print(" failures:", bad, "  smallest gap to the next modulus: %s (%s)" % (nstr(mingap[0], 3), mingap[1]))
