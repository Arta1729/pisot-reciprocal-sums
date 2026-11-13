# Code for "Reciprocal sums of Pisot linear recurrences"

Scripts accompanying **"Reciprocal sums of Pisot linear recurrences: effective rounding and the
asymptotics of the error"** by Artatrana Suna and Prasanta Kumar Ray. Version 1.5.

For an integer recurrence with Pisot characteristic polynomial, `S_n = sum_{m>=n} 1/u_m` and
`E_n = 1/S_n - (u_n - u_{n-1})`. The scripts compute the constants and starting indices of the paper and
check the error expansion numerically.

## Requirements

Python 3.9 or newer (tested with 3.12). Install with

```bash
pip install -r requirements.txt     # mpmath 1.3.0, numpy
```

## Files and what they reproduce

| File | Paper items |
|---|---|
| `verify_pisot.py` | Table 1 (constants, `N_0^root`, least starting index `n_*`), Table 2 (root-free Perron index `N_0^exp`), Table 3 (`N_par`, `kappa_2`, floor rule for 2<=n<=59), and the inclusion-radius statements of Section 7 (`--radii`) |
| `check_error_asymptotics.py` | Table 4 (exact `E_n` vs root sum), Table 5 (residual vs quadratic term), Tribonacci sign statistics, uniformity in `k`, closed form of `E_n` for 3<=n<=k+1 |
| `check_observations.py` | maximum of `|E_n|`, Tribonacci deviation from the leading cosine, `(1/2-x_k)4^k`, `(1-mu)k^3` and `N_0^root(1-mu)`, and the single-conjugate-pair hypothesis of the density theorem |

## Usage

```bash
python3 verify_pisot.py            # writes tables_generated.txt (Tables 1-3)
python3 verify_pisot.py --radii    # inclusion radii and the margins 1-q (Section 7)
python3 check_error_asymptotics.py
python3 check_observations.py      # needs verify_pisot.py in the same directory
```

Rows of Tables 1 and 2 are: Fibonacci, Pell, (2;3,1), Tribonacci, (3;6,1,1), 8-bonacci, 30-bonacci,
(4;3,2,2,1), (5;2,2,1,1,1), Narayana (3;1,0,1), (4;1,0,0,1). `verify_pisot.py` also treats the Perron
examples (3;3,0,1), (4;3,0,0,1), (4;5,-1,0,1).

Reference outputs from a run with Python 3.12.3, mpmath 1.3.0, numpy 2.5.4 are in `expected_output/`.
Measured running times on one core: a few seconds for `verify_pisot.py`, about 15 s for
`check_error_asymptotics.py`, about 10 s for `check_observations.py`.

## What is rigorous and what is not

* **Proved in the paper, evaluated by the code:** the inequalities `delta_n < 1`, `H_n < 1`, `B(n) < 1/2`
  that define `N_0` and the enclosure formulas for `S_n` are theorems. The code evaluates them.
* **Not interval arithmetic:** all arithmetic is mpmath multiprecision floating point. Rounding errors are
  controlled by the working precision (at least 60 digits beyond what is needed) and by the relative
  inflations `1 + 1e-30`; they are not formally tracked. The values of `n_*`, the Perron indices and all
  other numbers are therefore **computer-assisted, not formally certified**.
* **Numerical evidence only:** the maximum-error list, the Tribonacci sign statistics, the behaviour in `k`
  and the single-conjugate-pair check (it uses double-precision roots refined by Newton's method).

## Citation and licence

See `CITATION.cff`. Released under the MIT licence (see `LICENSE`).
