# vnorm — Third Follow-up to Codex's Review: Findings Log

**Purpose:** Record of the follow-up to Codex's review of `a3a527b..0eb69e3` (the
second follow-up, logged in `2026-10-30-2224-codex-review-followup-2.md`): each
finding's verdict, the fixes with their commits, and the validation.

**Repository:** `dkahle/vnorm`, branch `main` (local; not pushed).
**Toolchain:** R 4.6.1 (arm64 macOS), CmdStan 2.38.0 (cmdstanr 0.9.0), mpoly 1.1.2.900.

---

## 1. Verdicts

Both findings reproduced as reported.

| # | Codex finding | Reproduction | Commit |
|---|---|---|---|
| 1 (P2) | `86a2236`'s row scaling also changed the Newton projection's step | `project_onto_variety_newton(c(0, 0), mp("1 + x - 2 y + 4 x y"))` warned "Tolerance not met (residual = 1.00e+00)" and returned about (0, 0). The Moore–Penrose step is (0.12, 0.24); the row-scaled helper gave about 1e-16. | `b982f65` |
| 2 (P2) | `bbfebbe` took zeros from floating-point elimination for proof of dependence | `dependent_equations()` was `FALSE` for Codex's five-variable system `p` and `TRUE` for `1e-20 * p` | `ccaa2c9` |

## 2. Fixes

### 2.1 Newton projection — `b982f65`

- **Cause.** The Newton step is H^+ ∇ on the Hessian of g², which can be singular.
  When ∇ leaves H's column space, row scaling turns the least-squares step into a
  weighted one. At the origin above, the scaled gradient was orthogonal to the
  scaled column space, so the step was 0.
- **My earlier judgment.** I considered this case in `86a2236` and judged it a
  measure-zero edge case. The origin is an ordinary starting point.
- **Fix.** `svd_pseudoinverse()` is the plain SVD pseudoinverse with the tolerance
  cutoff, and the projection uses it. `pdvnorm_pseudoinverse()` applies the row
  scaling around it, for the density and the rejection kernel only. Those need J⁺g
  only where g lies in J's column space (m ≤ n), or don't scale at all (m > n).

### 2.2 Exact dependence in integers — `ccaa2c9`

Rounding and underflow in mpoly's double arithmetic produce exact zeros too.
Codex's 1e-20-scaled system underflowed. The exact stage now runs only on integers
below 2^53, where double arithmetic is exact, and is inconclusive (reporting
independence) whenever it cannot stay there.

- **`integer_jacobian_rows()`.**
  - Every double is an odd integer times a power of two (`binary_parts()`).
  - Each equation is divided by the gcd of its coefficients, which changes no rank:
    `0.1 x + 0.2 y` becomes `x + 2 y`, since 0.2 is exactly twice 0.1.
  - J is differentiated term by term, because mpoly's `deriv()` rounds
    `coef * power`.
  - It returns `NULL` once any coefficient reaches 2^53.
- **`polynomial_matrix_rank()`.** It returns `NA` before any update whose
  coefficients could reach 2^53. Each coefficient of `p a - f b` is at most
  |p||a| + |f||b|, with |.| the sum of absolute coefficients. Its time cap stays.
- **What is still proven dependent.**
  - Systems with integer coefficients.
  - Decimals in exact ratios.
  - Rescalings that keep those ratios, including any common factor on equal
    coefficients (so `1e-20 * mp(c("x + y + z", "(x + y + z)^2"))`) and powers of
    two.
- **What is now inconclusive.**
  - Decimals without exact ratios (0.1 and 0.3).
  - Integers whose products pass 2^53.
  - Codex's system at every scale. Its coefficients 1 and 1 + 1e-8 already need
    integers near 2^52.
- **A spurious warning.** R's `%%` warns of "probable complete loss of accuracy in
  modulus" for quotients past 4.5e15, though it is exact below 2^53. The gcd
  therefore computes `a - floor(a / b) * b` itself.

## 3. Notes and remaining limitations

- **Inconclusive dependent systems.** On the Stan path these fail or hang as they
  did before the check existed. The rejection path's window check follows its
  kernel, as before.
- **Rescaling.** Rescaling by a factor whose products with the coefficients round
  (0.1 times `x + 3 y`) can break an exact dependence. The check then reports
  independence. A verdict of dependence is always exact for the coefficients as
  given, so underflow and lost precision can no longer produce one.
- **Earlier items.** The second follow-up's §3 items still stand. The work-cap item
  now reads: large dependent systems, or ones whose integers grow past 2^53, are
  not caught.

## 4. Validation

Live runs used an installed copy of the final tree, in a fresh R process.

| Check | Result |
|---|---|
| `project_onto_variety_newton(c(0, 0), mp("1 + x - 2 y + 4 x y"))` | (0.06698716, 0.61602512), residual 3.8e-9, no warning (before: residual 1, warning) |
| Codex's `1e-20 * p` on the Stan path (2000 draws): max \|`lp__` − `pdvnorm`\| | 7.8e-7 (J's condition number is 5e8) |
| Same, KS of each coordinate vs N(0, 0.1^2), the exact target (g = J x, so J^-1 g = x) | Seeds 24–31 (40 tests): the p-values pass a KS test against U(0, 1) (p = 0.61), 7.5% below 0.05; joint χ²₅ tests of \|x\|²/0.01 give p = 0.09 to 0.97. Seeds 1 and 21–23 (20 tests): smallest p = 0.002 (x5, seed 23), where ESS was 1850–2090 of 2000 and a 500-draw subsample gave p = 0.55 |
| Dependent systems: `(x+y+z, (x+y+z)^2)` Stan, the same rescaled by 1e-20 Stan, rejection; `(x+y, 2x+2y)` Stan / rejection | error, error, error; error / sampled |

## 5. Verification

- `devtools::test()` on the final tree: **718 passed, 0 failed, 0 warnings, 1
  skipped** (the env-gated install test). That is up from 692.
- Each fix's new tests were run against the pre-fix copies of its R files and failed
  there:
  - The Newton test's warning and residual expectations.
  - The rescaling test's 1e-20 case.
  - The 2^53 bound in `polynomial_matrix_rank()`.
  - `rvnorm()` on Codex's rescaled system.
- **Shadowed targeted runs.** Some targeted test runs in the earlier rounds used
  `testthat::test_file(load_package = "none")` in an R session that still held a
  global prototype of `polynomial_matrix_rank()`, which shadowed the package's copy.
  - The prototype matched the package's algorithm at the time, so no reported
    result changes.
  - The full-suite runs and `R CMD check` resolve the namespace first, or run in
    fresh processes, so they were unaffected.
  - That session's globals were cleared before this round's runs.
- CRAN-style `R CMD check` (`check-local.R` settings, tests and examples run) on the
  final tree: **Status: OK** (0 errors, 0 warnings, 0 notes).
  - testthat summary: `[ FAIL 0 | WARN 2 | SKIP 3 | PASS 716 ]`.
  - The 2 warnings and 3 skips are the same as in the second follow-up: cmdstanr's
    E-BFMI notice, one `skip_on_cran()` test, and the install tests.
