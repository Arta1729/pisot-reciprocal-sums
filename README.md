# Code for "Reciprocal sums of Pisot linear recurrences: effective rounding and the asymptotics of the error"

Artatrana Suna and Prasanta Kumar Ray.

This archive contains the Python scripts used for the numerical results of the paper (Tables 1-5 and the
numerical observations of Section 6). All computations use **mpmath multiprecision floating-point arithmetic**.
No step uses interval arithmetic, so every numerical result is *computer-assisted*, not formally certified
(see "What is and is not certified" below).

## Requirements

Python 3.9 or newer and `mpmath` 1.3 (the version used for the paper):

```bash
pip install -r requirements.txt
```

## Files and what they reproduce

| File | Reproduces | Approx. run time |
|------|------------|------------------|
| `verify_pisot.py` | Table 1 (constants, `N_0^root`, least starting index `n_*`), Table 2 (root-free index `N_0^exp` of the Perron-type theorem), Table 3 (`kappa_2`, `alpha_2`, `N_par`) | about 1 minute |
| `check_error_asymptotics.py` | Table 4 (exact `E_n` against the root sum), Table 5 (residual against the quadratic term), Tribonacci sign/density statistics, uniformity in `k` for `k`-generalized Fibonacci numbers, closed forms for `E_n`, `3 <= n <= k+1`, and `E_{k+2}` | about 4 minutes |
| `verify_claims.py` | Remaining numerical claims of Section 6: the maxima of `|E_n|` for the eleven recurrences, the floor rule for `2 <= n <= 59` in the seven second-order rows, `|E_n| < 1/2` for 30-bonacci up to `n = 450`, the values of `(1-mu)k^3`, and the check that the non-dominant roots of maximal modulus form exactly one non-real conjugate pair for `F^(k)` (`4 <= k <= 40`), `P^(k)` (`4 <= k <= 15`) and `(4;1,0,0,1)` | about 2 minutes |

Run times are approximate (measured on a standard cloud machine).

## Usage

```bash
python verify_pisot.py                # full run: Tables 1-3; writes tables_generated.txt
python verify_pisot.py --radii        # inclusion radii and margins 1-q for the eleven recurrences
python check_error_asymptotics.py     # error-expansion tables and observations
python verify_claims.py               # remaining Section 6 claims; exit status 0 if all pass
python verify_claims.py --quick       # skips k = 160 and the n <= 450 scan of 30-bonacci
```

`verify_pisot.py` takes no other arguments. Rows are, in the order of Tables 1 and 2: Fibonacci, Pell, (2;3,1),
Tribonacci, (3;6,1,1), 8-bonacci, 30-bonacci, (4;3,2,2,1), (5;2,2,1,1,1), Narayana (3;1,0,1), (4;1,0,0,1);
the three additional Perron-condition recurrences (3;3,0,1), (4;3,0,0,1), (4;5,-1,0,1) are appended automatically.

## What is and is not certified

* **Controlled by precision margins and safety factors, not by interval arithmetic:** roots of the
  characteristic polynomial (200 digits), a-posteriori inclusion radii `k|Phi(z_i)|/prod|z_i - z_j|`
  (these are floating-point residuals; a printed radius of 0.0 means the residual underflowed at the working
  precision), the constants `C_0`, `C_1`, `q`, `d_1` (inflated or deflated by a relative `1e-30`), the thresholds
  `delta_n < 1`, `H_n < 1`, `B(n) < 1/2`, and the enclosures of `S_n` for `n < N_0` (working precision at least
  60 digits beyond what is needed).
* **Not a proof:** the checks in `check_error_asymptotics.py` and `verify_claims.py` (they compare exact
  rounding errors with the formulas of the paper at finitely many indices).
* The proofs in the paper hold for every Pisot recurrence; the scripts only certify finitely many cases
  numerically, as described in Section 6 of the paper.

## Citation and licence

See `CITATION.cff`. Released under the MIT licence (see `LICENSE`).
