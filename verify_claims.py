#!/usr/bin/env python3
"""
verify_claims.py -- numerical claims of the manuscript that are not produced by
verify_pisot.py or check_error_asymptotics.py.

Checks
  1. Section 6 (iii): max |E_n| over 3 <= n <= 59 (3 <= n <= 119 for 30-bonacci)
     for the eleven recurrences of Table 1, to four decimals.
  2. Table 3 caption: the floor rule of Corollary 1.2 holds for 2 <= n <= 59
     in all seven second-order rows.
  3. Section 6 (i) / Table 1: |E_n| < 1/2 for 2 <= n <= 450 for 30-bonacci.
  4. Remark on growth: (1 - mu) k^3 for k-generalized Fibonacci numbers,
     k = 10, 20, 40, 80, 160, against the values quoted in the paper.
  5. Remark on the density hypothesis: for F^(k) (4 <= k <= 40), P^(k)
     (4 <= k <= 15) and (4;1,0,0,1), the roots of maximal modulus among the
     non-dominant roots form exactly one non-real conjugate pair.

Arithmetic: mpmath multiprecision floating point (NOT interval arithmetic).
Working precision is chosen from the size of u_M, and the truncation point M of
the reciprocal sums is chosen so that the neglected tail changes E_n by far less
than 1e-30 (the error in E_n is about u_n times the relative tail).  The results are
computer-assisted, in the same sense as in Section 6 of the paper.

Usage:   python verify_claims.py            # all checks (a few minutes)
         python verify_claims.py --quick    # skips k = 160 and the n <= 450 scan
Exit status is 0 if every check passes, 1 otherwise.
"""
import sys
import math
from mpmath import mp, mpf, polyroots, floor, log10

QUICK = "--quick" in sys.argv
FAILED = []


def report(ok, label, detail=""):
    print(("[PASS] " if ok else "[FAIL] ") + label + (("  " + detail) if detail else ""), flush=True)
    if not ok:
        FAILED.append(label)


def sequence(c, N):
    """u_n for n = 1..N (exact integers), with u_{2-k..0}=0, u_1=1."""
    k = len(c)
    u = [0] * (k - 1) + [0, 1]          # u[i] = u_{i-(k-1)}
    off = k - 1
    for n in range(2, N + 1):
        u.append(sum(c[j] * u[off + n - 1 - j] for j in range(k)))
    return (lambda m: u[off + m])


def dominant_root(c):
    """Real dominant root alpha of z^k - c_1 z^{k-1} - ... - c_k (Newton on [1, 2+max|c|])."""
    mp.dps = 60
    co = [1] + [-x for x in c]
    lo, hi = mpf(1), mpf(2 + max(abs(x) for x in c))
    f = lambda x: sum(co[i] * x ** (len(c) - i) for i in range(len(c) + 1))
    for _ in range(300):
        mid = (lo + hi) / 2
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def errors(c, nmax):
    """Return (E, u) with E[n] = S_n^{-1} - (u_n - u_{n-1}) for 2 <= n <= nmax, plus u_1 case n=1."""
    alpha = dominant_root(c)
    la = float(log10(alpha))
    M = 2 * nmax + int(math.ceil(60 / la)) + 40   # error in E_n is ~ u_n * (tail/S_n), so M-n must exceed n + 60/la
    mp.dps = int((M * la)) + 70
    u = sequence(c, M)
    T = {}
    s = mpf(0)
    for m in range(M, 0, -1):
        s += mpf(1) / u(m)
        T[m] = s
    E = {n: 1 / T[n] - (u(n) - u(n - 1)) for n in range(1, nmax + 1)}
    return E, u, T


def check_max_E():
    cases = [
        ("Fibonacci", (1, 1), 59, "0.2646"),
        ("Pell", (2, 1), 59, "0.0778"),
        ("(2;3,1)", (3, 1), 59, "0.0295"),
        ("Tribonacci", (1, 1, 1), 59, "0.2182"),
        ("(3;6,1,1)", (6, 1, 1), 59, "0.0216"),
        ("8-bonacci", (1,) * 8, 59, "0.2490"),
        ("30-bonacci", (1,) * 30, 119, "0.2500"),
        ("(4;3,2,2,1)", (3, 2, 2, 1), 59, "0.0617"),
        ("(5;2,2,1,1,1)", (2, 2, 1, 1, 1), 59, "0.1056"),
        ("Narayana (3;1,0,1)", (1, 0, 1), 59, "0.3841"),
        ("(4;1,0,0,1)", (1, 0, 0, 1), 59, "0.4411"),
    ]
    print("\n1. max |E_n| over 3 <= n <= nmax (Section 6 (iii))")
    for name, c, nmax, expected in cases:
        E, u, T = errors(c, nmax)
        mx = max(abs(E[n]) for n in range(3, nmax + 1))
        got = "%.4f" % float(mx)
        report(got == expected, "max|E_n| %-20s" % name, "computed %s, paper %s" % (got, expected))


def check_parity():
    print("\n2. Floor rule, 2 <= n <= 59 (Table 3)")
    rows = [(1, 1), (2, 1), (3, 1), (3, -1), (6, -1), (4, -2), (5, -3)]
    for c in rows:
        E, u, T = errors(c, 59)
        bad = []
        for n in range(2, 60):
            t = u(n) - u(n - 1)
            pred = t - (n % 2) if c[1] > 0 else t - 1
            if floor(1 / T[n]) != pred:
                bad.append(n)
        report(not bad, "floor rule (2;%d,%d)" % c, "failures at n=%s" % bad if bad else "")


def check_30bonacci():
    print("\n3. 30-bonacci: |E_n| < 1/2 for 2 <= n <= 450")
    nmax = 120 if QUICK else 450
    E, u, T = errors((1,) * 30, nmax)
    bad = [n for n in range(2, nmax + 1) if abs(E[n]) >= mpf(1) / 2]
    mx = max(abs(E[n]) for n in range(2, nmax + 1))
    report(not bad, "30-bonacci |E_n|<1/2, 2<=n<=%d" % nmax, "max|E_n| = %s" % mp.nstr(mx, 12))


def roots_by_modulus(c):
    mp.dps = 60
    co = [1] + [-x for x in c]
    rr = sorted(polyroots(co, maxsteps=5000, extraprec=800), key=lambda z: -abs(z))
    return rr


def check_growth():
    print("\n4. (1 - mu) k^3 for k-generalized Fibonacci numbers")
    expected = {10: "31.8", 20: "40.2", 40: "41.3", 80: "40.7", 160: "40.2"}
    for k, ex in expected.items():
        if QUICK and k == 160:
            continue
        rr = roots_by_modulus((1,) * k)
        mu = abs(rr[1])
        got = mp.nstr((1 - mu) * k ** 3, 3)
        report(got == ex, "(1-mu)k^3, k=%d" % k, "computed %s, paper %s" % (got, ex))
    print("   (4*pi^2 = %s)" % mp.nstr(4 * mp.pi ** 2, 5))


def check_single_pair():
    print("\n5. Roots of maximal modulus: exactly one non-real conjugate pair")

    def pair_info(c):
        rr = roots_by_modulus(c)
        mu = abs(rr[1])
        L = [z for z in rr[1:] if abs(abs(z) - mu) < mpf(10) ** -30]
        return len(L) == 2 and all(abs(z.imag) > mpf(10) ** -20 for z in L)

    bad = [k for k in range(4, 41) if not pair_info((1,) * k)]
    report(not bad, "F^(k), 4<=k<=40", "failures: %s" % bad if bad else "")
    bad = [k for k in range(4, 16) if not pair_info((2,) + (1,) * (k - 1))]
    report(not bad, "P^(k), 4<=k<=15", "failures: %s" % bad if bad else "")
    report(pair_info((1, 0, 0, 1)), "(4;1,0,0,1)")


if __name__ == "__main__":
    check_max_E()
    check_parity()
    check_30bonacci()
    check_growth()
    check_single_pair()
    print()
    if FAILED:
        print("%d check(s) FAILED: %s" % (len(FAILED), "; ".join(FAILED)))
        sys.exit(1)
    print("All checks passed.")
