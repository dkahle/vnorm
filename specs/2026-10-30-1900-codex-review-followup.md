# vnorm — Follow-up to Codex's Review of Round 3: Findings Log

**Purpose:** Record of the follow-up to Codex's review of the 33 unpushed commits
(`1fba652`..`df60784`, mostly the round-3 work in
`2026-10-30-1700-package-review-round-3.md`): each finding's verdict, the fixes with
their commits, defects found while fixing them, and the validation.

**Repository:** `dkahle/vnorm`, branch `main` (local; not pushed).
**Toolchain:** R 4.6.1 (arm64 macOS), CmdStan 2.38.0 (cmdstanr 0.9.0), mpoly
1.1.2.900, jsonlite 2.0.0.

---

## 1. Verdicts

All five findings reproduced as reported.

| # | Codex finding | Reproduction | Commits |
|---|---|---|---|
| 1 (P1) | Rejection sampling ignores `pdvnorm()`'s rank/support rule | `pdvnorm()` was 0 for 10 of 10 draws of `mp(c("x+y+z", "(x+y+z)^2"))` | `fd346aa`, `57c74ed` |
| 2 (P2) | `digits = 15` still rounds coefficients | Residual `-0.3125` instead of `0`; density 0.0302 instead of 3.99 | `9e039b7`, `26e9960`, `95bee2a`, `e687d15` |
| 3 (P2) | One-variable batch projections come back 1 x k | `dim()` was `1 3`, in all four functions | `56da04d` |
| 4 (P2) | `variety_solve()` misreads variables named like template quantities | `g = 0` for `mp("g - 1")` | `cd8c958` |
| 5 (P2) | Polynomial-list cache key depends on coefficients | The keys for `p` and `q` differed | `bc7ac9f` |
| Tests | `rvnorm_or_skip()` turns sampling errors into skips | Agreed | `32778d2` |

## 2. Fixes

### 2.1 Rejection kernel and dependent systems — `fd346aa`, `57c74ed`

- **Rejection kernel.** It now uses `pdvnorm_pseudoinverse()` and
  `pdvnorm_outside_column_space()`, the same rule with the same tolerances.
  - Live check: log kernel minus `log pdvnorm()` is constant to within 1.5e-12 over
    200 random points for each of four systems.
  - The kernel is `-Inf` exactly where `pdvnorm()` is 0, including 14 points on the
    singular axes of `V(x^2 - 1, y^2 - 1)`.
  - Consequences:
    - Square systems with a singular `J` no longer call `solve(J)`, which crashed.
    - Linearly dependent systems, such as `(x + y, 2 x + 2 y)`, keep `g` in `J`'s
      column space and a positive density, and are now sampled.
    - Nonlinearly dependent systems have density 0 almost everywhere. They are
      rejected up front by `generic_jacobian_check()`, which evaluates `J` and `g` at
      three fixed generic points. Before, the sampler exhausted `max_batches` and
      then advised enlarging `w` or `sd`.
- **Stan path.** If the Jacobian's generic rank is below `min(m, n)`, then
  `J \ g`, `J'((JJ') \ g)` and `(J'J) \ (J'g)` are singular everywhere.
  - Codex's system failed with "Location parameter[1] is nan", then "No chains
    finished successfully".
  - `mp(c("x + y", "2 x + 2 y"))` was still running after ten minutes when
    interrupted.
  - `rvnorm()` now stops with an error that names the cause, before any model is
    built. `homo = FALSE` (no `J`) and `code_only = TRUE` are unaffected.
  - This is not a convention change. These inputs never sampled, and the
    `pdvnorm()` rule approved in round 3 is unchanged.
- **Docs.** `?rvnorm` Details describe both paths.

### 2.2 Precision — `9e039b7`, `26e9960`, `95bee2a`, `e687d15`

Facts established while fixing:

- **R's parser is not exact at 17 digits on this platform** (arm64, where long
  double is double). It misreads 17-significant-digit decimals: 36,591 of 200,000
  random doubles via `as.numeric()`, and 3,721 of 20,000 via `parse()`. So
  evaluation that goes through strings cannot be exact at any `digits` setting.
- **mpoly drops near-1 coefficients.** It prints any coefficient within `sqrt(eps)`
  of 1 as 1.
- **cmdstanr rounds Stan data.** `cmdstanr::write_stan_json()` writes 15
  significant digits (jsonlite's `digits = NA`): 8,993 of 10,000 random doubles
  came back changed.
- **Stan reads 17 digits exactly.** Stan read all 10,000 of the same doubles back
  exactly from a JSON file written with 17 significant digits. stanc and C++ also
  read 17-digit literals exactly.

Fixes:

- **`9e039b7`: `poly_function()` builds the function from the terms.**
  Coefficients are inlined as numbers, and the result does not depend on
  `options(digits)`.
  - Speed:
    - Matrix input is evaluated in one vectorized pass: a 201 x 201 grid takes
      about 1 ms, down from 0.16 s.
    - A single point takes about 0.7 µs, up from 0.35 µs.
  - It accepts a `varorder` with extra variables. mpoly's `as.function()` refused
    one, which broke two things:
    - `rejection_sampler(n, mp("y - x^2"))`, because `|grad g|^2 = 4 x^2 + 1` has
      no `y`.
    - `geom_variety()` grids for polynomials in `x` alone.
- **`26e9960`: the Lagrange projection builds its Lagrangian with mpoly arithmetic
  on `x0`.** The formula string rounded `x0` to 15 significant digits, which missed
  the projection by 2e-5 at coordinates around 1.2e11.
- **`95bee2a`: `create_stan_code()` prints expressions itself** (`stan_poly_expr()`,
  `stan_real()`).
  - Each coefficient and window bound is written with 15 significant digits when R
    reads that back exactly, and with 17 otherwise.
  - Output for short coefficients is unchanged, with two exceptions:
    - Polynomial lists keep their own term order (`y-x`, not `-x+y`).
    - Derived coefficients that are not short now print exactly (`0.2^2` as
      `0.040000000000000008`).
- **`e687d15`: `rvnorm()` writes its Stan data itself.** It uses jsonlite at 17
  significant digits and passes the file path. This covers template coefficients,
  user-compiled models, `si` and `w`. jsonlite moved to Imports; it was already
  installed with cmdstanr.

### 2.3 Batch projections — `56da04d`

`projected_rows()` reshapes the `apply()` output to one row per input point, with
the input's dimnames. New tests cover all four functions with matrix and data-frame
input, and check that two-variable rows match pointwise projections.

### 2.4 `variety_solve()` — `cd8c958`

Template draws are recognized by the template's own `x`/`y`/`z` columns: when all of
them are present, the coordinates come from those columns. Generated code samples
the variables under their own names and has no other `x`/`y`/`z` columns.

### 2.5 Cache key — `bc7ac9f`

Equations are sorted by their canonical monomial keys. Equations with the same
structure keep their input order, an order shared by the key, the code and
`get_coefficients_data()`. Keys change for systems whose equation order changes.
The cache is per session (`compile_stan_code()` writes to `tempdir()`), so nothing
persistent is invalidated.

### 2.6 Tests — `32778d2`

`rvnorm_or_skip()` is removed. CI is unaffected: it installs cmdstanr but not
CmdStan, so `skip_if_no_cmdstan()` still skips there.

## 3. Notes and remaining limitations

- **`stan_real()`'s 15-digit test.** It relies on R reading decimals of at most 15
  significant digits exactly, which takes one correctly rounded scaling. If that
  check fails, the 17-digit form is used, and stanc reads it exactly. Literals at
  extreme magnitudes are long: `1e-20 / 3` prints as
  `0.0000000000000000000033333333333333333`.
- **Reading draws back can be off by one ulp.** The CSV goes through R's parser,
  which here was off by 0.0156 at x around 1.2e14. The validation tolerances below
  allow for this. It affects comparisons, not sampling.
- **`generic_jacobian_check()` uses three fixed points.** A system whose singular
  locus contained all three would be misclassified. That is implausible for
  ordinary systems, but such a system could be constructed.
- **Round-3 §4 items still stand.** The exception is the template-path
  `variety_solve()` collision, now fixed. Generated code still fails to compile for
  variables named `g`, `ndg`, and so on. On the template path, `rvnorm()`'s tibble
  still keeps such a variable and drops the transformed `g`.
- **Upstream mpoly.** The recommended `as.function()` precision fix in mpoly still
  stands, but vnorm no longer depends on it. A second item for mpoly:
  `as.function.mpoly()` requires `varorder` to match the polynomial's own
  variables exactly.

## 4. Validation

Live runs used an installed copy of the final tree (template binaries in use), in a
fresh R process. `a = 123456789012345.67`, and `g` is the program's transformed
residual compared with R's `x - a`. One ulp at that magnitude is 0.0156; a
15-digit coefficient is 0.3125 off.

| Check | Result |
|---|---|
| Template `1_1_vn`, `mp("x") - a`, `sd = 1`: max \|g - (x - a)\| | 0.0156 (cmdstanr's writer: 0.328) |
| Generated code, same: max \|g - (x - a)\| | 0.0156 |
| Generated polynomial list, coefficients `1/3` and `0.1 + 0.2`, `sd = 1e-3`: max \|`lp__` - `pdvnorm`\| | 1.5e-13 |
| `compile_stan_code(p)`, then `rvnorm(q, user_compiled = TRUE)` after a coefficient change that reordered equations: same key; max \|`lp__` - `pdvnorm`\| | TRUE; 8.0e-13 |
| `variety_solve(mp("g - 1"))` (template path) | `g = 1` |
| Rejection kernel vs `log pdvnorm()`: constant-offset sd (linear dependent / overdetermined / underdetermined / singular axes) | 1.1e-14 / 3.2e-16 / 7.8e-14 / 1.5e-12; same support everywhere |
| Rejection draws, `(x + y, x - 2 y)` and `(x, y, x + y)`, `sd = 0.2` vs N(0, 0.2^2) (KS, n = 4000) | p = 0.79, 0.78 / 0.44, 0.37 |
| `geom_variety(poly = mp("x^2 - 1"))` | lines at x = -1, 1 |
| Dependent `mp(c("x+y+z", "(x+y+z)^2"))`, Stan path / rejection path | errors naming the cause |
| Round-3 checks rerun: V1 `lp__` diff; V1 spreads; V2; template; circle KS; `variety_solve(t)`; tibble names; rejection class; `thin` | 2.8e-10; (0.100, 0.011); 2.8e-12; 3.5e-13; p = 0.945; `t = 1`; unchanged; data.frame; 1000 rows |

## 5. Verification

- `devtools::test()`: **670 passed, 0 failed, 0 warnings, 1 skipped** (the env-gated
  install test). That is up from 570.
- Each fix's new regression tests were run against the pre-fix code (with the
  fix's R files stashed) and failed there. The one exception is the extended stanc
  syntax check, which guards the new long literals and passes either way.
- CRAN-style `R CMD check` (`check-local.R` settings, tests and examples run):
  **Status: OK** (0 errors, 0 warnings, 0 notes).
