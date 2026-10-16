"""
verify_pisot.py -- computations behind the tables in
"Effective Bounds for Reciprocal Sums of Higher-Order Linear Recurrences".

Rigorous parts (python-flint ball arithmetic, arb/acb):
  * isolation of the roots of Phi_k, enclosures of alpha, mu, d_1;
  * lower bound for M_r = min_{|z|=r}|Phi_k| (grid + Lipschitz margin, Remark 3.3);
  * upper bounds for the L^1-means I_r, J_r (grid + Lipschitz margin, Remark 3.3);
  * the thresholds N_0 (all inequalities evaluated in ball arithmetic and accepted only if certainly true);
  * the discrete winding number.
Floating-point parts (marked): per-point lower bounds |Phi(z_j)| are rounded to double precision after
subtracting a relative slack of 1e-12, and the grid sums use math.fsum with a final relative slack 1e-9.
Non-rigorous parts: the scan over r (numpy), the high-precision check of the identity (mpmath), the
root-separation study.
"""
import math, sys, time
import numpy as np
from flint import fmpz_poly, arb, acb, fmpq, ctx
import mpmath as mp

ctx.prec = 200

# ---------------------------------------------------------------- basic objects
def phi_poly(c):
    k = len(c); co = [0]*(k+1); co[k] = 1
    for i, ci in enumerate(c, 1): co[k-i] = -ci
    return fmpz_poly(co)

def constants(c):
    k = len(c); P = phi_poly(c)
    rts = [z for z, m in P.complex_roots()]
    assert len(rts) == k
    rts.sort(key=lambda z: float(z.abs_upper()), reverse=True)
    alpha = rts[0].real
    assert alpha > 1, "dominant root not certainly > 1"
    mu = max((z.abs_upper() for z in rts[1:]), key=lambda x: float(x))
    assert mu < 1, "non-dominant roots not certainly inside the unit disc"
    dphi = k*alpha**(k-1) - sum((k-i)*c[i-1]*alpha**(k-i-1) for i in range(1, k))
    d1 = alpha**(k-2)/dphi
    assert d1 > 0
    return k, alpha, mu, d1

def u_list(c, N):
    """returns list U with U[n]=u_n for n=0..N (u_0=0,u_1=1)"""
    k = len(c); v = [0]*(k-1) + [1]          # v[i]=u_{i-(k-2)} : v[k-2]=u_0, v[k-1]=u_1
    for n in range(2, N+1):
        v.append(sum(c[j]*v[-1-j] for j in range(k)))
    return v[k-2:]

# ---------------------------------------------------------------- threshold (ball arithmetic)
def N0_ball(k, alpha, d1, r, I, J, nmax=10**6):
    """least n>=2 with delta_n<1, H_n<1, B(n)<1/2, all certainly true in ball arithmetic"""
    half = arb(1)/2
    def ok(n):
        d = I*r**(n+k-1)/(d1*alpha**n)
        if not d < 1: return False
        H = d*(alpha-1)/(alpha*(1-d)*(1-r/alpha**2))
        if not H < 1: return False
        B = J*r**(n+k-2) + d1*alpha**(n-1)*(alpha-1)*H/(1-H)
        return bool(B < half)
    if ok(2): return 2
    hi = 4
    while not ok(hi):
        hi *= 2
        if hi > nmax: return None
    lo = hi//2
    while hi-lo > 1:
        m = (lo+hi)//2
        if ok(m): hi = m
        else: lo = m
    assert all(ok(n) for n in range(hi, hi+25))      # sanity (monotone by Lemma 5.1)
    return hi

# ---------------------------------------------------------------- grid certificate on |z|=r
def grid_cert(c, r_q, N):
    """r_q: fmpq. Returns dict with certified M_lb, I_plus, J_plus, winding data."""
    k = len(c); old = ctx.prec; ctx.prec = 100
    r = arb(r_q)
    zs = []
    P_abs_lb = np.empty(N); mids = np.empty(N, dtype=complex)
    twopi = 2*arb.pi()
    for j in range(N):
        z = r*acb(0, twopi*j/N).exp()
        p = acb(1)
        for i in range(k):                    # Horner for z^k - c1 z^{k-1} - ... - ck
            p = p*z - c[i]
        lbj = p.abs_lower()
        P_abs_lb[j] = float(lbj.mid()) - float(lbj.rad())
        mids[j] = complex(float(p.real.mid()), float(p.imag.mid()))
    P_abs_lb *= (1 - 1e-12)
    L = k*r**(k-1) + sum((k-i)*abs(c[i-1])*r**(k-i-1) for i in range(1, k))
    eps = float(L.upper())*math.pi*float(r.upper())/N*(1 + 1e-9)
    m = float(P_abs_lb.min())
    ok = m > eps
    M_lb = m - eps if ok else 0.0
    out = {"ok": ok, "M_lb": M_lb, "min_grid": m, "eps": eps}
    if ok:
        denom = P_abs_lb - eps
        rf = float(r.upper())
        th = 2*np.pi*np.arange(N)/N
        one_minus = np.abs(1 - rf*np.exp(1j*th)) + math.pi*rf/N
        out["I_plus"] = math.fsum(1.0/denom)/N*(1+1e-9)
        out["J_plus"] = math.fsum(one_minus/denom)/N*(1+1e-9)
        ang = np.angle(np.roll(mids, -1)/mids)
        out["winding"] = int(round(float(ang.sum())/(2*np.pi)))
        out["max_step"] = float(np.abs(ang).max())
    ctx.prec = old
    return out

def to_arb(x): return arb(x)

# ---------------------------------------------------------------- float scan to choose r (not rigorous)
def scan_best_r(c, k, alpha, mu, d1, mode, N=400000, npts=25):
    a = float(alpha.mid()); m_ = float(mu.upper()); d = float(d1.mid())
    th = np.arange(N)*2*np.pi/N
    best = (10**18, None)
    for t in np.linspace(0.02, 0.98, npts):
        r = m_ + (1-m_)*t
        z = r*np.exp(1j*th)
        p = z**k - sum(c[j]*z**(k-1-j) for j in range(k))
        ap = np.abs(p)
        L = k*r**(k-1)+sum((k-1-j)*abs(c[j])*r**(k-2-j) for j in range(k-1))
        eps = L*math.pi*r/N*1.01
        if ap.min() <= 2*eps: continue
        if mode == "sup":
            I = 1/(ap.min()-eps); J = (1+r)*I
        else:
            I = np.mean(1/(ap-eps)); J = np.mean((np.abs(1-z)+math.pi*r/N)/(ap-eps))
        n = N0_float(k, a, d, r, I, J)
        if n is not None and n < best[0]: best = (n, r)
    return best

def N0_float(k, a, d1, r, I, J):
    def ok(n):
      try:
        dl = I*r**(n+k-1)/(d1*a**n)
        if dl >= 1: return False
        H = dl*(a-1)/(a*(1-dl)*(1-r/a**2))
        if H >= 1: return False
        return J*r**(n+k-2)+d1*a**(n-1)*(a-1)*H/(1-H) < 0.5
      except OverflowError:
        return True
    if ok(2): return 2
    hi = 4
    while not ok(hi):
        hi *= 2
        if hi > 10**6: return None
    lo = hi//2
    while hi-lo > 1:
        m = (lo+hi)//2
        if ok(m): hi = m
        else: lo = m
    return hi

# ---------------------------------------------------------------- Perron closed form (Theorem 1.2, corrected)
def perron_N0(c, nmax=100000):
    k = len(c); c1 = c[0]; S = sum(abs(x) for x in c[1:])
    assert c[-1] != 0 and c1 > 1 + S
    r = (arb(S+1)/c1)**(arb(1)/(k-1))
    am = arb(c1) if all(x >= 0 for x in c[1:]) else arb(c1-1)
    Ms = 1 - r**k
    half = arb(1)/2
    def ok(n):
        dt = r**(n+k-1)*(c1+1+r)**(k-1)/(Ms*am**(n+k-2))
        if not dt < 1: return False
        Ht = dt/((1-dt)*(1-r/am**2))
        if not Ht < 1: return False
        B = r**(n+k-2)*(1+r)/Ms + r**(n+k-1)/(Ms*(1-r/am**2)*(1-dt)*(1-Ht))
        return bool(B < half)
    n = 2
    hi = 2
    while not ok(hi):
        hi = hi*2 if hi > 2 else 4
        if hi > nmax: return None
    lo = max(2, hi//2)
    if ok(2): return 2
    while hi-lo > 1:
        m = (lo+hi)//2
        if ok(m): hi = m
        else: lo = m
    assert all(ok(n) for n in range(hi, hi+25))
    return hi

# ---------------------------------------------------------------- identity check (mpmath, not rigorous)
def identity_check(c, alpha_f, Nmax):
    mp.mp.dps = int(Nmax*math.log10(alpha_f)) + 80
    start = 2*Nmax + 400
    U = u_list(c, start+400)
    S = sum(mp.mpf(1)/U[m] for m in range(start, start+401))
    errs = {}
    for n in range(start-1, 0, -1):
        S += mp.mpf(1)/U[n]
        if n <= Nmax: errs[n] = abs(1/S - (U[n]-U[n-1]))
    mp.mp.dps = 30
    first = min(n for n in range(1, Nmax+1) if all(errs[m] < mp.mpf(1)/2 for m in range(n, Nmax+1)))
    mx = float(max(errs[n] for n in range(3, Nmax+1)))
    return first, mx, float(errs[1]), float(errs[2])

N0_SUP_MID = {(1,1):11,(2,1):5,(3,1):3,(1,1,1):16,(6,1,1):2,(1,)*8:97,(1,)*30:4961,(3,2,2,1):10,(2,2,1,1,1):16,(1,0,1):29,(1,0,0,1):117}
ROWS = [("Fibonacci", [1,1]), ("Pell", [2,1]), ("A", [3,1]), ("Tribonacci", [1,1,1]), ("A", [6,1,1]),
        ("8-bonacci", [1]*8), ("30-bonacci", [1]*30), ("B", [3,2,2,1]), ("B", [2,2,1,1,1]),
        ("Narayana", [1,0,1]), ("Pisot", [1,0,0,1])]

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    NGRID = int(sys.argv[2]) if len(sys.argv) > 2 else 10**6
    if which in ("all", "table"):
        sel = [int(x) for x in sys.argv[3].split(",")] if len(sys.argv) > 3 else range(len(ROWS))
        for name, c in [ROWS[i] for i in sel]:
            t0 = time.time()
            k, alpha, mu, d1 = constants(c)
            res = {}
            r_mid = fmpq(int(round((1+float(mu.upper()))/2*10**6)), 10**6)
            g = grid_cert(c, r_mid, NGRID)
            assert g["ok"] and g["winding"] == k-1 and g["max_step"] < math.pi/2, (name, c, g)
            r = arb(r_mid)
            Msup = arb(g["M_lb"])
            res["sup_mid"] = N0_ball(k, alpha, d1, r, 1/Msup, (1+r)/Msup)
            M32 = (alpha-r)*(r-mu)**(k-1)
            res["lb32"] = N0_ball(k, alpha, d1, r, 1/M32, (1+r)/M32)
            res["L1_mid"] = N0_ball(k, alpha, d1, r, arb(g["I_plus"]), arb(g["J_plus"]))
            for mode, key in (("sup", "sup_best"), ("L1", "L1_best")):
                nb, rb = scan_best_r(c, k, alpha, mu, d1, mode)
                rq = fmpq(int(round(rb*10**6)), 10**6)
                gb = grid_cert(c, rq, NGRID)
                if gb["ok"] and gb["winding"] == k-1 and float(mu.upper()) < rb:
                    rr = arb(rq)
                    if mode == "sup":
                        Mb = arb(gb["M_lb"]); res[key] = N0_ball(k, alpha, d1, rr, 1/Mb, (1+rr)/Mb)
                    else:
                        res[key] = N0_ball(k, alpha, d1, rr, arb(gb["I_plus"]), arb(gb["J_plus"]))
                    res[key+"_r"] = float(rq)
            print(name, c, "alpha=%s mu<=%s d1=%s r=%s M_lb=%.5f I+=%.4f J+=%.4f 1/M=%.3f" % (
                  alpha.str(10, radius=False), mu.str(6, radius=False), d1.str(8, radius=False), float(r_mid),
                  g["M_lb"], g["I_plus"], g["J_plus"], 1/g["M_lb"]), res, "%.0fs" % (time.time()-t0), flush=True)
    if which in ("all", "perron"):
        for c in ([3,1], [6,1,1], [3,0,1], [3,0,0,1], [11]+[0]*8+[1], [5,-1,0,1], [5,1,1,1]):
            print("Perron closed form", c, perron_N0(c), flush=True)
    if which in ("all", "identity"):
        extra = [[3,0,1], [3,0,0,1], [5,-1,0,1]]
        sel = [int(x) for x in sys.argv[3].split(",")] if len(sys.argv) > 3 else None
        allrows = ROWS + [("extra", e) for e in extra]
        for name, c in (allrows if sel is None else [allrows[i] for i in sel]):
            k, alpha, mu, d1 = constants(c)
            N0c = N0_SUP_MID.get(tuple(c), 120) + 10
            Nmax = max(40, N0c)
            print(name, c, "identity: first n, max err(3..Nmax), err1, err2 =", identity_check(c, float(alpha.mid()), Nmax), flush=True)
