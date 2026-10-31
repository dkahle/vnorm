# vnorm — Second Follow-up to Codex's Review: Findings Log

**Purpose:** Record of the follow-up to Codex's review of `df60784..a3a527b` (the
first follow-up, logged in `2026-10-30-1900-codex-review-followup.md`): each
finding's verdict, the fixes with their commits, defects found while fixing them,
and the validation.

**Repository:** `dkahle/vnorm`, branch `main` (local; not pushed).
**Toolchain:** R 4.6.1 (arm64 macOS), CmdStan 2.38.0 (cmdstanr 0.9.0), mpoly
1.1.2.900, jsonlite 2.0.0; Python 3.12 and Apple clang as independent decimal
readers.

---

## 1. Verdicts

Both findings reproduced as reported.

| # | Codex finding | Reproduction | Commits |
|---|---|---|---|
| 1 (P2) | The Jacobian precheck takes numerical rank deficiency at three fixed points for proof of dependence, and its rank check can throw | `generic_jacobian_check(mp(c("x^100-1", "y^100-1")))` gave `full_rank = FALSE, density_vanishes = TRUE`; `rejection_sampler(10, p, sd = .01, w = c(.99, 1.01))` stopped; for `q` (below), "system is computationally singular" | `b2c46b1`, `bbfebbe` |
| 2 (P2) | `stan_real()` certifies 15-digit literals with R's parser | `stan_real(-0x1.1f6af1a526408p-32)` gave `-0.000000000261404920558413`; Python, clang and jsonlite read that as `-0x1.1f6af1a526407p-32` | `c7ba6a8` |

Here `q = mp(c("0.000000001 x^2-1", "y^2-1", "0.000000001 x^2+y^2-2"))`.

## 2. Fixes

### 2.1 Pseudoinverse at full rank — `b2c46b1`

The throw in finding 1 came from `pdvnorm_pseudoinverse()`.
- **What it did.** It found J's rank from the SVD. At full rank, for non-square J,
  it then formed `JJ'` or `J'J` and called `solve()`. That squares J's condition
  number, so `solve()` refused full-rank J once its condition number passed about
  1e8.
- **Other failures it caused.**
  - `pdvnorm(c(.5, .5), q, sd = .1)` failed.
  - The rejection kernel failed at every point of
    `mpolyList(1e-9 x, y, 1e-9 x + y)`. The kernel used `MASS::ginv()` until
    `fd346aa` in the first follow-up switched it to this helper, so the sampler
    failure is a regression from that round. The `solve()` branches date from the
    original helper (`7adcb68`).
- **Fix.** The pseudoinverse now comes from the SVD at every rank, the formula the
  rank-deficient case already used. It is accurate to about cond(J) eps, where the
  normal equations were accurate to cond(J)^2 eps.

### 2.2 Exact dependence check — `bbfebbe`

`dependent_equations()` replaces `generic_jacobian_check()` on both paths. It
rejects a system only when the equations are proven dependent.

1. **Fast path.** At the same three fixed points, J's rows are scaled to unit size.
   If the smallest singular value is clearly nonzero (above sqrt(eps) relative),
   that settles independence. This uses singular values alone, with no inverse.
   Badly scaled systems such as `x^100 - 1`, whose gradient is about 1e-20 at
   x = 0.6, pass here.
2. **Exact stage.** Otherwise the rank at generic points is computed exactly, by
   fraction-free elimination on the derivative polynomials
   (`polynomial_matrix_rank()`, mpoly arithmetic, exact zero tests). Measured costs:
   - Typical systems: under 20 ms.
   - Dense 3 x 3 independent system of degree 6: 0.4 s.
   - Dependent 3 x 3 system with a product equation (degree 12): 2.4 s.
3. **Work cap.** Elimination without division grows fast: a dependent 4 x 4 example
   would take about 26 s. Past 1e5 term products (a few seconds) the check gives up,
   calls the result inconclusive, and rejects nothing.

On the rejection path, a proven-dependent system is stopped only if the kernel's own
support rule gives density 0 at three points of the proposal window. The old check
applied `pdvnorm()`'s rule at points of [-1, 1]^n instead. A residual that leaves
J's column space only by rounding stays within `pdvnorm()`'s tolerance, and such a
system is still sampled. Example: `mp(c("0.1 x + 0.2 y - 0.3", "0.3 x + 0.6 y - 0.9"))`
is exactly rank 1, and 0.1 × 0.9 ≠ 0.3 × 0.3 in doubles. An exact test of
[J | g] would have rejected it.

### 2.3 Literal check — `c7ba6a8`

`stan_real()` now checks its 15-digit candidate with `jsonlite::parse_json()`, which
reads numbers with C's `strtod()` and rounds correctly, as stanc and C++ do.
jsonlite was already in Imports.

R's parser misreads some 15-digit decimals that need large powers of ten, so the
old check failed in both directions. The sweep below found 172 of 20,013 literals
that read back as another double. It also found 161 that took 17 digits because R
misread a 15-digit form that was exact; those are now short. Common decimals
(`0.1`, `2.3`, `0.001`, `1e-22`, ...) are unchanged.

### 2.4 Row scaling in `pdvnorm()`'s rank cutoff — `86a2236` (found during validation)

Live validation of finding 1 showed Stan's `lp__` and `pdvnorm()` differing by up
to 2.1 for `mp(c("x^100-1", "y^100-1"))`, on 371 of 1000 draws.
- **Cause.** `pdvnorm_pseudoinverse()` and `pdvnorm_outside_column_space()` cut off
  singular values below 2 eps times the largest. That test depends on how each
  equation happens to be scaled.
- **Example.** At (1.59, 2.78), J = diag(8.7e21, 9.1e45) has full rank, but the
  cutoff dropped the smaller singular value. With it went the x component of
  J^-1 g: log density 3.51 instead of 2.24.
- **Consequences.** Stan's `lp__` had the right value. The rejection kernel, which
  uses the same helpers, sampled the truncated density. Multiplying an equation by
  1e-30 has the same effect at any point. This dates from `7adcb68`.
- **Fix.** For m ≤ n, the rows are scaled by their largest entries before the SVD,
  and the pseudoinverse is scaled back. Rescaling an equation changes neither J's
  rank nor J^+ g for g in J's column space. The support rule and the density
  therefore change only where the cutoff misjudged the rank. Overdetermined systems
  keep their rows as they are, since there the row scale weights the least squares.

## 3. Notes and remaining limitations

- **Dependence in floating point.** The exact stage tests the coefficients mpoly
  arithmetic produces.
  - Integer and dyadic coefficients cancel exactly.
  - For other decimals, rounding can hide a dependence. The system is then not
    rejected, and the Stan path fails or hangs as it did before the check existed.
  - Rounding can also cancel to exactly zero a minor whose true value lies below the
    rounding level. That happens only for a J that is singular to working precision
    everywhere.
- **Work cap.** Large dependent systems past the cap are not caught.
- **jsonlite's correctness.** It rests on the platform's `strtod()`, which rounds
  correctly in glibc, macOS and the Windows UCRT. Where it did not, a wrong literal
  could pass the check, as R's parser let one pass before.
- **Inits.** cmdstanr also writes inits at 15 significant digits. That moves a start
  point, not the target. `rvnorm()` passes user inits through unchanged.
- **Overdetermined systems near ill-conditioned points.** The Stan programs still
  use `(J'J) \ (J'g)`, so near the variety of `q` (condition number about 1e9),
  `lp__` and `pdvnorm()` differ by up to 1.1e-6.
- **Earlier items.** The first follow-up's §3 items still stand.

## 4. Validation

Live runs used an installed copy of the final tree (template binaries in use), in a
fresh R process.

| Check | Result |
|---|---|
| Codex's literal: program for `a * x - a * 2^50`, `fixed_param` at x = 2^50 (exact init file), normalized residual / density | new 0 / 3.9894228 (= `pdvnorm()`); old literal 0.2226724 / 0.3343713 |
| `stan_real()` sweep, 20,013 doubles (32-bit mantissas, exponents −300..300, weighted toward 1e-33..1e-9 and 1e12..1e33, plus typed decimals), read by Python's `float()` | old: 172 wrong, 161 needlessly long; new: 0 wrong |
| Same two literals compiled by clang | 15-digit → `...407p-32`; 17-digit → `a` |
| `mp(c("x^100-1", "y^100-1"))`, Stan path, `sd = 0.01`: max \|`lp__` − `pdvnorm`\| | 7.9e-14 (2.1 before `86a2236`) |
| Same, Stan draws vs rejection draws in [0.95, 3]^2 (two-sample KS, x / y; 4 seeds) | p = 0.08 / 0.60, 0.60 / 0.11, 0.63 / 0.78, 0.35 / 0.40 |
| `q` on the Stan path, init near (sqrt(1e9), 1): max \|`lp__` − `pdvnorm`\| | 1.1e-6 (normal equations in the program) |
| Rejection draws: `x^100` system in [0.97, 1.03]^2 / `mpolyList(1e-9 x, y, 1e-9 x + y)` vs N(0, 0.1^2) (KS x, y) | min `pdvnorm` 5.86 / p = 0.48, 0.70 |
| Dependent systems: `(x+y+z, (x+y+z)^2)` Stan / rejection; `(x+y, 2x+2y)` Stan / rejection | error / error; error / sampled |
| First follow-up's checks rerun (template and generated `g`, cmdstanr writer, list `lp__`, user-compiled `lp__`, kernel offsets, KS, `geom_variety`, `variety_solve`) | 0.0156, 0.328, 0.0156, 1.5e-13, 5.6e-13, same support with offset sds 1.4e-14 to 1.5e-12, unchanged p-values, lines at ±1, `g = 1` |
| Round-3 checks rerun: V1 `lp__` diff; V1 spreads; V2; template; circle KS; `variety_solve(t)`; tibble names; rejection class; `thin` | 2.8e-10; (0.100, 0.011); 2.8e-12; 3.5e-13; p = 0.945; `t = 1`; unchanged; data.frame; 1000 rows (all as before) |

## 5. Verification

- `devtools::test()` on the final tree: **692 passed, 0 failed, 0 warnings, 1
  skipped** (the env-gated install test). That is up from 670.
- Each fix's new tests were run against the pre-fix copies of its R files, sourced
  into the loaded namespace, and failed there. There are two exceptions. The kernel
  support test in `bbfebbe` guards the design and passes either way. The
  rescaled-equation comparison in `86a2236` also passes either way, because the old
  code truncated both sides alike; its two direct comparisons fail on the old code.
- CRAN-style `R CMD check` (`check-local.R` settings, tests and examples run), on
  the final tree and on the tree before `86a2236`: **Status: OK** (0 errors, 0
  warnings, 0 notes).
  - testthat summary: `[ FAIL 0 | WARN 2 | SKIP 3 | PASS 690 ]`.
  - The 2 warnings are cmdstanr's "E-BFMI not computed" notice from two `rvnorm()`
    tests that draw 2 samples per chain through the installed templates.
  - The skips are one `skip_on_cran()` test and the two install tests.
  - The previous round's check had the same 2 warnings and 3 skips.
