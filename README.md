# verify_pisot.py

Code accompanying the paper **"Effective Bounds for Reciprocal Sums of Higher-Order Linear Recurrences"**
by Artatrana Suna and Prasanta Kumar Ray.

The script reproduces Tables 1-3 of the paper (constants, thresholds N_0, and the root-free bound of
Theorem 1.3) and the numerical checks quoted in Sections 5 and 6.

## Requirements

Python 3.9 or newer. Install the packages with

```bash
pip install -r requirements.txt
```

(`python-flint`, `numpy`, `mpmath`; the results in the paper were produced with python-flint 0.9.)

## Usage

```bash
python3 verify_pisot.py table 1000000        # Tables 1 and 2: all 11 recurrences, N = 10^6 grid points
python3 verify_pisot.py table 1000000 0,6    # only rows 0 and 6 (0-based; 6 = 30-bonacci)
python3 verify_pisot.py perron               # Table 3: root-free bound of Theorem 1.3
python3 verify_pisot.py identity             # numerical check of the identity; gives n_min
python3 verify_pisot.py identity 0 0,1,9     # identity check for selected rows only
```

Rows are numbered from 0 in the order of Tables 1 and 2: Fibonacci, Pell, (2;3,1), Tribonacci, (3;6,1,1),
8-bonacci, 30-bonacci, (4;3,2,2,1), (5;2,2,1,1,1), Narayana (3;1,0,1), (4;1,0,0,1). For `identity`, rows 11-13
are three additional recurrences satisfying Perron's condition: (3;3,0,1), (4;3,0,0,1), (4;5,-1,0,1).

Example output (`perron`):

```
Perron closed form [3, 1] 6
Perron closed form [6, 1, 1] 5
Perron closed form [3, 0, 1] 12
Perron closed form [3, 0, 0, 1] 18
```

Running times on a laptop: roughly 10-20 seconds per row for small k, about one minute per radius for the
30-bonacci row (the full `table` run takes several minutes), and a few minutes for `identity`.

## What is rigorous and what is not

* **Rigorous (ball arithmetic, python-flint):** isolation of the roots of the characteristic polynomial,
  enclosures of alpha, mu and d_1, and the threshold inequalities delta_n < 1, H_n < 1, B(n) < 1/2 (each accepted
  only if certainly true).
* **Double precision with explicit slack:** the per-point lower bounds |Phi_k(z_j)| on the grid, the sums defining
  I+ and J+, and the winding number. The slack factors 1 - 1e-12 and 1 + 1e-9 exceed the accumulated rounding error
  by several orders of magnitude, but these steps are not interval-certified (see Section 6 of the paper).
* **Not a proof:** the scan over the radius r (it only chooses r, which is then certified as above), and the
  high-precision floating-point check of the nearest-integer identity for 1 <= n <= N_0 + 10 (mpmath).

## Citation and licence

See `CITATION.cff`. Released under the MIT licence (see `LICENSE`).
