#!/usr/bin/env python3
"""
verify_pisot.py -- code accompanying "Reciprocal sums of Pisot linear recurrences: effective rounding
and the asymptotics of the error" (A. Suna, P. K. Ray).

Reproduces the computer-assisted tables of the paper (root-sum constants, N_0^root, least starting
index n_*, root-free Perron index N_0^exp, and the parity indices N_par).

For each recurrence (k; c_1..c_k):
  1. roots of Phi_k at 200 digits; a-posteriori Smith inclusion disks (rad_i = k|p(z_i)|/prod|z_i-z_j|);
  2. constants  C0 >= sum_j |d_j| ,  C1 >= sum_j |d_j||1-alpha_j| ,  q >= mu  (inflated by 1+1e-30);
  3. N0_root = least n>=2 with delta_n<1, H_n<1, B(n)<1/2   (Corollary on N_0, Section 4);
  4. for every n < N0_root: enclosure of S_n = sum_{m>=n} 1/u_m (exact finite sum plus the tail enclosure
     T_M/(1+delta_M) <= tail <= T_M/(1-delta_M)), then test whether the enclosure of S_n^{-1}
     lies strictly inside (u_n-u_{n-1}-1/2, u_n-u_{n-1}+1/2).
Arithmetic: mpmath multiprecision floating point, NOT interval arithmetic; rounding errors are controlled by
the working precision (>= 60 digits beyond what is needed) and the 1e-30 inflations, not formally tracked.
The results are therefore computer-assisted, not formally certified.

Usage:  python3 verify_pisot.py            # all tables, also written to tables_generated.txt
        python3 verify_pisot.py --radii    # inclusion radii and the margins 1-q
"""
import sys, math
from mpmath import mp, mpf, mpc, polyroots, polyval, log10, findroot

CASES = [
 ("Fibonacci",(1,1)),("Pell",(2,1)),("(A)",(3,1)),("Tribonacci",(1,1,1)),("(A)",(6,1,1)),
 ("8-bonacci",(1,)*8),("30-bonacci",(1,)*30),("(B)",(3,2,2,1)),("(B)",(2,2,1,1,1)),
 ("Narayana",(1,0,1)),("Pisot",(1,0,0,1))]
INF = mpf(1)+mpf(10)**-30

def seq(c, N):
    k=len(c); u=[0]*(k-1)+[0,1]   # u_{2-k..0}=0, u_1=1 ; index shift: u[i] = u_{i-(k-1)}
    off=k-1
    for n in range(2,N+1):
        u.append(sum(c[j]*u[off+n-1-j] for j in range(k)))
    return lambda m: u[off+m]

def constants(c):
    mp.dps=200; k=len(c)
    coeffs=[1]+[-x for x in c]
    rts=polyroots(coeffs,maxsteps=2000,extraprec=800)
    p=lambda z: polyval(coeffs,z)
    dp=lambda z: polyval([ (k-i)*coeffs[i] for i in range(k)],z)
    rad=[]
    for i,z in enumerate(rts):
        den=mpf(1)
        for j,w in enumerate(rts):
            if j!=i: den*=abs(z-w)
        rad.append(k*abs(p(z))/den)
    ia=max(range(k),key=lambda i: rts[i].real if abs(rts[i].imag)<mpf(10)**-100 else -1)
    alpha=rts[ia].real
    assert alpha>1 and all(abs(rts[i])+rad[i]<1 for i in range(k) if i!=ia), "not Pisot?"
    oth=[i for i in range(k) if i!=ia]
    q=max(abs(rts[i])+rad[i] for i in oth)
    d=lambda z: z**(k-2)/dp(z)
    D=sum(abs(d(rts[i])) for i in oth); D1=sum(abs(d(rts[i]))*abs(1-rts[i]) for i in oth)
    d1=1/sum((j*c[j-1])*alpha**(1-j) for j in range(1,k+1))
    # consistency check d1 = alpha^{k-2}/Phi'(alpha)
    assert abs(d1-d(alpha))<mpf(10)**-100*abs(d1)
    return dict(k=k,alpha=alpha,q=q*INF,d1=d1/INF,C0=D*INF,C1=D1*INF/q,rad=max(rad),roots=rts,ia=ia,mu=max(abs(rts[i]) for i in oth))
# note: |e_{m-1}-e_m| <= sum|d_j||1-a_j||a_j|^{m-1} <= (D1/q) q^m

def Bfun(P,n,C0,C1,q):
    a=P['alpha'];d1=P['d1']
    dl=C0*q**n/(d1*a**n)
    if dl>=1: return None
    H=dl*(a-1)/(a*(1-dl)*(1-q/a**2))
    if H>=1: return None
    B=C1*q**n+C0*(a-1)**2*q**n/(a**2*(1-dl)*(1-H)*(1-q/a**2))
    return dl,H,B

def N0_from(P,C0,C1,q,nmax=10**6):
    lo=2
    # all three quantities decrease -> bisection on first n where all hold
    def ok(n):
        r=Bfun(P,n,C0,C1,q); return r is not None and r[2]<mpf(1)/2
    hi=2
    while not ok(hi):
        hi*=2
        if hi>nmax: return None
    lo=hi//2 if hi>2 else 2
    while lo<hi:
        mid=(lo+hi)//2
        if ok(mid): hi=mid
        else: lo=mid+1
    return lo

def sharp_index(c,P,N0):
    k=len(c); a=P['alpha']; d1=P['d1']; q=P['q']; C0=P['C0']
    M=2*N0+60
    mp.dps=int(2*M*float(log10(a)))+80
    u=seq(c,M+2)
    # tail enclosure at M
    rho=C0*q**M/(d1*a**M)
    assert rho<mpf(10)**-6
    A=1/(d1*a**(M-1)*(a-1))
    T_lo=A/(1+rho); T_hi=A/(1-rho)
    S_lo=T_lo; S_hi=T_hi
    fail=[]; undecided=[]
    for n in range(M-1,0,-1):
        t=mpf(1)/u(n)
        S_lo+=t; S_hi+=t
        if n<N0:
            tgt=u(n)-u(n-1)
            lo=1/S_hi; hi=1/S_lo
            if lo>tgt-mpf(1)/2 and hi<tgt+mpf(1)/2: pass
            elif hi<=tgt-mpf(1)/2 or lo>=tgt+mpf(1)/2: fail.append(n)
            else: undecided.append(n)
    nstar = (max(fail)+1) if fail else 1
    return nstar,fail,undecided,M

def perron(c):
    """Root-free index N_0^exp (Perron case): depends only on (k; c_1..c_k); alpha_- = c_1 - 1."""
    k=len(c); S=sum(abs(x) for x in c[1:]); c1=c[0]
    assert c1>1+S
    r=(mpf(S+1)/c1)**(mpf(1)/(k-1)); Ms=1-r**k
    am=c1-1
    th=1-r/am**2
    C0=r**(k-1)/Ms; C1=(1+r)*r**(k-2)/Ms
    n=2
    while True:
        dl=C0*(c1+2*mpf(S)/am)*(r/am)**n
        if dl<1:
            H=dl/((1-dl)*th)
            if H<1:
                B=C1*r**n+C0*r**n/(th*(1-dl)*(1-H))
                if B<mpf(1)/2: return r,n
        n+=1


def kappa2_and_Theta(c,P):
    """k=2 helper: kappa_2, alpha_2, and N_par (Remark on effective parity)."""
    mp.dps=200
    a=P['alpha']; z=[P['roots'][i] for i in range(2) if i!=P['ia']][0].real
    kap=-(a-z)/(a**2-z); q=P['q']; C0,C1=P['C0'],P['C1']; d1=P['d1']
    def Theta(n):
        r=Bfun(P,n,C0,C1,q)
        if r is None: return None
        dl,H,B=r
        return ((a-1)/a)**2*C0*q**n/(1-dl)*(dl/(1-q**2/a**3)+H/((1-H)*(1-q/a**2)))
    n=2
    while True:
        th=Theta(n); r=Bfun(P,n,C0,C1,q)
        if th is not None and r[2]<mpf(1)/2 and abs(kap)*abs(z)**(n-1)>th: return kap,z,n
        n+=1


def floor_rule_check(c, nmax=59):
    """k=2: check floor(S_n^{-1}) = u_n-u_{n-1}-[n odd] (c_2>0) or u_n-u_{n-1}-1 (c_2<0) for 2<=n<=nmax."""
    from mpmath import floor
    al=constants(c)['alpha']
    D=2*int(nmax*float(log10(al)))+80
    mp.dps=D+60
    M=nmax+1
    u=seq(c,M)
    while u(M)<=mpf(10)**(D+25):
        M*=2; u=seq(c,M)
    bad=[]
    for n in range(2,nmax+1):
        S=sum(mpf(1)/u(m) for m in range(n,M+1))
        t=u(n)-u(n-1)
        fl=int(floor(1/S))
        exp=t-(1 if n%2 else 0) if c[1]>0 else t-1
        if fl!=exp: bad.append(n)
    mp.dps=200
    return bad


def radii_report():
    """Per recurrence: largest inclusion radius, margin 1-q (q before inflation), disjointness of disks."""
    mp.dps=200
    out=[]
    for name,c in CASES:
        P=constants(c); mp.dps=200
        rts=P['roots']; k=len(c)
        coeffs=[1]+[-x for x in c]
        rad=[]
        for i,z in enumerate(rts):
            den=mpf(1)
            for j,w in enumerate(rts):
                if j!=i: den*=abs(z-w)
            rad.append(k*abs(polyval(coeffs,z))/den)
        disjoint=all(abs(rts[i]-rts[j])>rad[i]+rad[j] for i in range(k) for j in range(i+1,k))
        out.append((name,k,max(rad),1-P['q']/INF,disjoint))
        print("%-11s k=%-3d max_radius=%s  1-q=%s  disjoint=%s"%(name,k,mp.nstr(max(rad),3),mp.nstr(1-P['q']/INF,4),disjoint),flush=True)
    return out

if __name__=="__main__":
    import sys
    if len(sys.argv)>1 and sys.argv[1]=="--radii":
        radii_report(); sys.exit(0)
    out=[]
    def pr(*a):
        s=" ".join(str(x) for x in a); print(s,flush=True); out.append(s)
    mp.dps=60
    pr("== Perron root-free index N0^exp (alpha_- = c_1 - 1) ==")
    PERRON=[(3,1),(6,1,1),(3,0,1),(3,0,0,1),(5,-1,0,1)]
    for c in PERRON: pr(c, "r=%s"%mp.nstr(perron(c)[0],4), "N0_exp=%d"%perron(c)[1])
    pr("== root-sum certificate and least starting index n_* ==")
    pr("%-11s %-4s %8s %8s %6s  %s"%("family","k","N0_root","n_star","M","fails<N0 / undecided"))
    rows=[]
    for name,c in CASES+[("(A)",x) for x in PERRON[2:]]:
        P=constants(c); C0,C1,q=P['C0'],P['C1'],P['q']
        N0=N0_from(P,C0,C1,q); mp.dps=200
        ns,fail,und,M=sharp_index(c,P,N0)
        mp.dps=200
        pr("%-11s k=%-3d %8d %8d %6d  fails=%s und=%s"%(name,len(c),N0,ns,M,fail[:8],und[:8]))
        pr("    c=%s alpha=%s mu=%s d1=%s C0=%s C1=%s"%(c if len(c)<9 else "(1)^%d"%len(c),mp.nstr(P['alpha'],8),mp.nstr(P['mu'],8),mp.nstr(P['d1'],8),mp.nstr(C0,6),mp.nstr(C1,6)))
    pr("== effective parity index N_par (k=2) and direct floor-rule check, 2<=n<=59 ==")
    for c in [(1,1),(2,1),(3,1),(3,-1),(6,-1),(4,-2),(5,-3)]:
        try:
            P=constants(c); kap,z,n=kappa2_and_Theta(c,P)
            bad=floor_rule_check(c)
            pr(c,"alpha_2=",mp.nstr(z,8),"kappa_2=",mp.nstr(kap,8),"N_par=",n,"floor-rule violations:",bad)
        except AssertionError: pr(c,"not Pisot")
    open("tables_generated.txt","w").write("\n".join(out))
