# check_error_asymptotics.py
# 1) Table 4: exact rounding error E_n vs the leading term sum_j kappa_j alpha_j^(n-1) (Theorem 5.1).
# 2) Table 5: E_n - leading term vs the quadratic term of Theorem 5.5.
# 3) Observation (viii): sign statistics of E_n for Tribonacci, n = 20..330 (Theorem 5.6(b)).
# 4) Table 6: k-generalized Fibonacci numbers, E_n vs the root sum at small n (Observation (ix)).
# 5) Proposition 6.2: exact formula for E_n, 3<=n<=k+1, for k=3,5,8,12.
# Requires mpmath and numpy. The whole script takes several minutes.
from mpmath import mp, mpf, mpc, polyroots, polyval, re, fabs
mp.dps=60
def run(c, ns):
    k=len(c)
    u={i:0 for i in range(2-k,1)}; u[1]=1
    N=max(ns)+400
    for n in range(2,N+1):
        u[n]=sum(c[j-1]*u[n-j] for j in range(1,k+1))
    coeffs=[1]+[-x for x in c]
    rts=polyroots(coeffs,maxsteps=500,extraprec=200)
    al=max(rts,key=lambda z:abs(z)); al=re(al)
    others=[z for z in rts if abs(z-al)>mpf(10)**-20]
    dp=lambda b: sum(1 for _ in [0]) and polyval([ (k-i)*coeffs[i] for i in range(k)],b)
    d=lambda b: b**(k-2)/dp(b)
    mu=max(abs(z) for z in others)
    print(c,"alpha",mp.nstr(al,8),"mu",mp.nstr(mu,6))
    for n in ns:
        S=sum(mpf(1)/u[m] for m in range(n,N+1))
        E=1/S-(u[n]-u[n-1])
        M=sum(d(b)*b**(n-1)*(al-b)**2/(al**2-b) for b in others)
        print(" n=%d E=%s Main=%s ratio=%s"%(n,mp.nstr(E,8),mp.nstr(re(M),8),mp.nstr(E/re(M),8)))
run([1,1],[3,4,5,8,12])
run([2,1],[3,4,6])
run([1,1,1],[3,6,10,16,24])
run([3,1],[3,5])
run([1,0,1],[4,10,20,30])
run([5,-1,0,1],[3,6])
run([3,-1],[4,6,10])

from mpmath import mp, mpf, polyroots, polyval, re, nstr, arg, pi
mp.dps=110
def setup(c,N):
    k=len(c); u={i:0 for i in range(2-k,1)}; u[1]=1
    for n in range(2,N+1): u[n]=sum(c[j-1]*u[n-j] for j in range(1,k+1))
    co=[1]+[-x for x in c]
    rts=polyroots(co,maxsteps=800,extraprec=400)
    al=re(max(rts,key=lambda z:abs(z)))
    oth=[z for z in rts if abs(z-al)>mpf(10)**-30]
    Dp=lambda b: polyval([(k-i)*co[i] for i in range(k)],b)
    d=lambda b: b**(k-2)/Dp(b)
    return k,u,al,oth,d,Dp
def terms(al,oth,d,n,Dp,k):
    d1=al**(k-2)/Dp(al)
    M=sum(d(b)*b**(n-1)*(al-b)**2/(al**2-b) for b in oth)
    Q=0
    for bi in oth:
        for bj in oth:
            lam=al*(al-1)**2*((al-1)/((al**2-bi)*(al**2-bj))-1/(al**3-bi*bj))
            Q+=d(bi)*d(bj)/d1*lam*(bi*bj/al)**n
    return re(M),re(Q)
def run2(c,ns):
    N=max(ns)+500
    k,u,al,oth,d,Dp=setup(c,N)
    print(c)
    for n in ns:
        S=sum(mpf(1)/u[m] for m in range(n,N+1)); E=1/S-(u[n]-u[n-1])
        M,Q=terms(al,oth,d,n,Dp,k)
        print(" n=%d E-M=%s Q=%s (E-M-Q)=%s"%(n,nstr(E-M,6),nstr(Q,6),nstr(E-M-Q,4)))
run2([1,1],[8,12,16])
run2([1,1,1],[10,16,24,32])
run2([1,0,1],[20,30,40])
run2([3,-1],[6,10])
# density for tribonacci (needs more digits for n up to 330)
mp.dps=260
c=[1,1,1]; N=330+500
k,u,al,oth,d,Dp=setup(c,N)
pos=0;tot=0;mism=0
for n in range(20,331):
    S=sum(mpf(1)/u[m] for m in range(n,N+1)); E=1/S-(u[n]-u[n-1])
    M,Q=terms(al,oth,d,n,Dp,k)
    tot+=1; pos+= (E>0)
    if (E>0)!=(M>0): mism+=1
print("trib n=20..330: frac E>0 =",pos/tot,"sign mismatch with Main:",mism)
z=[b for b in oth if b.imag>0][0]
print("theta/pi=",nstr(arg(z)/pi,12))


# ---- part 4: Table 6 ----
import numpy as np
from mpmath import mp, mpf, polyroots, polyval, re, nstr
mp.dps=150
def data(c,N):
    k=len(c); co=[1]+[-x for x in c]
    u={i:0 for i in range(2-k,1)}; u[1]=1
    for n in range(2,N+1): u[n]=sum(c[j-1]*u[n-j] for j in range(1,k+1))
    rts=polyroots(co,maxsteps=2000,extraprec=500)
    al=re(max(rts,key=lambda z:abs(z))); oth=[z for z in rts if abs(z-al)>mpf(10)**-25]
    Dp=lambda b: polyval([(k-i)*co[i] for i in range(k)],b)
    return k,co,u,al,oth,Dp
def E_of(u,n,N): 
    S=sum(mpf(1)/u[m] for m in range(n,N+1)); return 1/S-(u[n]-u[n-1])
def Main(k,al,oth,Dp,n): return re(sum(b**(k-2)/Dp(b)*b**(n-1)*(al-b)**2/(al**2-b) for b in oth))
print("Table 6: Uniformity test: exact E_n vs full root-sum Main_n for k-bonacci")
for name,c,ns in [("5b",[1]*5,[2,3,5,7,10,20]),("8b",[1]*8,[2,3,5,8,10,15,30,60]),("12b",[1]*12,[2,5,14,30,60,120]),("30b",[1]*30,[2,10,32,60,150])]:
    N=max(ns)+400 if len(c)<20 else max(ns)+1500
    if len(c)>=20: mp.dps=400
    k,co,u,al,oth,Dp=data(c,N)
    for n in ns:
        E=E_of(u,n,N); M=Main(k,al,oth,Dp,n)
        print(" %-4s n=%3d E=%s Main=%s  diff=%s"%(name,n,nstr(E,6),nstr(M,6),nstr(E-M,3)))


# ---- part 5: Proposition 6.2 ----
from mpmath import mp, mpf, nstr
mp.dps=120
for k in (3,5,8,12):
    N=k+700
    F={i:0 for i in range(2-k,1)}; F[1]=1
    for n in range(2,N+1): F[n]=sum(F[n-j] for j in range(1,k+1))
    R=sum(mpf(1)/F[m] for m in range(k+2,N+1))-mpf(2)**(1-k)
    print("k=%d  R/(3*4^-k)=%s"%(k,nstr(R/(3*mpf(4)**-k),8)))
    for n in range(3,k+2):
        S=sum(mpf(1)/F[m] for m in range(n,N+1))
        E=1/S-(F[n]-F[n-1])
        pred=-mpf(4)**(n-3)*R/(1+mpf(2)**(n-3)*R)
        assert abs(E-pred)<mpf(10)**-80,(k,n,E,pred)
    S=sum(mpf(1)/F[m] for m in range(k+2,N+1)); E=1/S-(F[k+2]-F[k+1])
    print("   n=3..k+1 formula exact to 1e-80; E_{k+2}=%s"%nstr(E,8))


# ---- part 6: Observation (ix) scan and the Remark 5.5 comparisons ----
print("Observation (ix): worst deviation |E_n - root sum|, k-bonacci, 2<=n<=60")
mp.dps=200
for k in list(range(2,13))+[16,20]:
    c=[1]*k; N=70+500
    kk,co,u,al,oth,Dp=data(c,N)
    w=0;wn=None
    for n in range(2,61):
        d=abs(E_of(u,n,N)-Main(kk,al,oth,Dp,n))
        if d>w: w=d;wn=n
    print(" k=%d worst=%s at n=%d"%(k,nstr(w,3),wn),flush=True)
print("Remark 5.5: root sum vs exact E at the argmax indices of Observation (i)")
mp.dps=150
for name,c,n in [("Fib",[1,1],3),("Trib",[1,1,1],4),("Nar",[1,0,1],3),("P4",[1,0,0,1],5),("B4",[3,2,2,1],4),("B5",[2,2,1,1,1],3),("A361",[6,1,1],3),("8b",[1]*8,10)]:
    N=n+400
    kk,co,u,al,oth,Dp=data(c,N)
    print(" %-5s n=%d E=%s Main=%s"%(name,n,nstr(E_of(u,n,N),5),nstr(Main(kk,al,oth,Dp,n),5)))
