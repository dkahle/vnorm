# vnorm — Package Review Round 3: Findings Log

**Purpose:** Record of the 2026-10-30 full package review: what was checked, what was
found, and what was done. Like the August log, this records work already landed; each
item names its commit.

**Repository:** `dkahle/vnorm`, branch `main` (local; not pushed).
**Review scope:** line-by-line code review against the July/August findings, full test
suite + CRAN-style R CMD check, live CmdStan validation of the generated programs
against `pdvnorm()`, and distributional checks of the samplers.
**Toolchain:** R 4.6.1, CmdStan 2.38.0 (cmdstanr 0.9.0), mpoly 1.1.2.900, roxygen2 8.0.0.

---

## 1. State at review start

No commits since the round-2 log (`e631bad`, 2026-08-15); clean tree. Baseline
`devtools::test()`: 501 passed, 0 failed, 1 skipped (env-gated install test).

## 2. Correctness fixes landed

### 2.1 Generated `mpolyList` Stan code used term-appearance variable order — `21af3b3`

`create_stan_code()` declared parameters and Jacobian columns in `mpoly::vars()` order,
so for `mp(c("y^2 + x^2 - 1", "y - x"))` a vector `sd` or a `Sigma` over variables
(`homo = TRUE`) applied to `(y, x)`, while `pdvnorm()` uses sorted order. Same class as
the July/August varorder bugs; this path (all `mpolyList` sampling) was missed.
Live check: `Sigma = diag(c(.1, .01)^2)` gave per-coordinate spreads `(0.011, 0.100)`
before and `(0.100, 0.011)` after; `lp__` vs `pdvnorm(log = TRUE)` max |diff| 670
before, 3e-10 after.

### 2.2 Polynomials evaluated with 7-significant-digit coefficients — `91bd4c8`, `17b7bb4`

Root cause is in **mpoly**: `as.function.mpoly()` builds its function by parsing
`print.mpoly()` output, and `print_term()` formats coefficients with
`format(coef, scientific = FALSE)`, i.e. `getOption("digits")` (7).
`as.function(mp("x - 1.23456789012345"))(0)` returns `-1.234568`. Consequences in
vnorm: `pdvnorm()`, the rejection sampler, the projection family, and
`geom_variety()` contours all used rounded coefficients, and generated Stan code
(`create_stan_code()`, which also prints through mpoly) baked rounded coefficients into
every `mpolyList` program and every polynomial outside the template catalog. The
precompiled templates get exact coefficients as data, so the sampler and the density
disagreed depending on the path.

- `91bd4c8` raised printing to 15 significant digits in `create_stan_code()` and writes
  integer literals of 10+ digits as reals (stanc rejects ints past 2^31 - 1).
- That commit assumed `pdvnorm()` was exact; the live `lp__` comparison showed the
  mismatch had moved (max |diff| 5e-12 before, 9e-4 after). `17b7bb4` routes every
  package `as.function()` call through `poly_function()`, which sets `digits = 15`
  around it. All paths now agree (see §5).
- **Upstream recommendation (not done here):** fix `as.function.mpoly()` in mpoly so
  evaluation does not depend on the display option (e.g. format coefficients with
  `digits = 17` or `deparse()` them). `poly_function()` can then be dropped.

### 2.3 Projection family derivatives — `d53d73a`, `cb10d4f`, `910a3c3`, `ca62562`

- `project_onto_variety_newton()` built its Hessian with `as.function()` and no
  `varorder`, so each row used its own appearance order (message silenced): wrong in
  every entry for `y^2 - x^3 - x` (Newton stopped at residual 0.67), malformed for
  `x^2 + y`, row-inconsistent for `x y - 1`. One-variable polynomials crashed
  `project_onto_variety()`, the gradient-descent variant (in a Hessian it never used),
  and Newton (`deriv()` returns a bare mpoly for one variable, so `lapply()` walked its
  terms). New helpers `gradient_function()` / `hessian_function()` (R/helpers.R).
- Gradient descent and Newton looped while `g^2 > tol` but warned when `|g| >= tol`,
  so default runs on the unit circle warned on success; the tests hid this with
  `suppressWarnings()`. Loops now test `|g|`.
- Newton errored on every linear polynomial in >= 2 variables (Hessian of `g^2` is
  `2 grad(g) grad(g)'`, singular). Directions now use the package pseudoinverse
  (unchanged at full rank; the exact orthogonal projection for linear `g`).
- `for (. in 1:n_correct)` ran two corrector steps when `n_correct = 0`.

### 2.4 `variety_solve()` — `a398370`, `e91b0e0`

- The zero-dimensionality guard compared equations with `length(vars)` (the variables
  to average), so averaging a subset failed with "not zero-dimensional".
- For a single polynomial in variables other than x/y/z, the precompiled path's fit has
  x/y/z column names, so `variety_solve(mp("t^2 - 2 t + 1"))` failed whenever templates
  were installed. Draw columns are mapped back with `check_and_replace_vars()`.

### 2.5 Rejection sampler — `e0ade47`, `ce70da8`

- Scalar `sd`, `homo = TRUE`, `mpolyList`: the kernel `g'(JJ')^{-1}g` needs `JJ'`
  invertible and crashed on every overdetermined system — including
  `rvnorm(rejection = TRUE)`, which only accepts scalar `sd`. All homoskedastic kernels
  now share `J^+ g` (scalar/vector/matrix scales give identical draws for equivalent
  scales; kernel verified proportional to `pdvnorm()` exactly). Also fixed: crash on
  one-variable systems with a linear member (`as.function()` rejects a bare constant
  mpoly with a `varorder`) and `diag()` reading a 1 x 1 matrix scale as an identity size.
- Proposal batches shrank to the number of draws still needed, so finishing took about
  `log(n) / rate` batches: 1000 draws on the unit sphere at `sd = .001` hit
  `max_batches` with 852 accepted. Batches are now sized from the running acceptance
  rate (cap 1e5); the box-misses-variety case still fails as fast as before. Draws for
  a given seed change beyond the first batch.

### 2.6 Smaller defects

- `rvnorm(rejection = TRUE)` returned a matrix for `"simple"` output; documented and
  Stan-path type is a data frame — `6b3e0e8`.
- Compiled-model cache key (`generate_model_name()`) depended on in-term variable order
  and on input order of same-degree terms, so `user_compiled = TRUE` reported "model not
  found" for `y x + ...` after compiling `x y + ...` — `b8893fa`.
- `rvnorm()` tibble output sorted template quantities (`dgx`, `dgy`, `g`, `ndg`) ahead
  of `x`, `y`; variables now come first, then model quantities in Stan order, then
  diagnostics — `f335676`.
- `pdvnorm()` failed ("incorrect number of dimensions") for a single point of a
  one-variable system — `09779af`.
- Docs: `rvnorm()` `n` was described as per chain (it is the total), `poly` as mpoly
  only, and several arguments by reference to rstan's `stan()`; `compile_stan_code()`'s
  cache note had been absorbed into `@param homo` and it lacked `@return` — `2347861`.

## 3. Judgment calls (user-approved 2026-10-30)

| Decision | Commit |
|---|---|
| `pdvnorm()` systems, `homo = TRUE`: density 0 where `J` loses rank and `g` leaves its column space, for square/underdetermined systems (revisits the August "convention" note). `V(x^2-1, y^2-1)` scored 15.9 at the origin, same as the solutions; singular loci showed as lines of isolated bright pixels. Overdetermined least squares unchanged. | `fe35d08` |
| `rvnorm(n, thin = k)` returns `n` draws (runs `ceiling(n * thin / chains)` iterations); previously `n / k`. | `b242ceb` |
| `rvnorm(rejection = TRUE)` default `w = 1.25` (as `rejection_sampler()`); `w = 1` cut 18% of the unit circle's mass at `sd = .05`. | `725719a` |
| Removed the open research note ("Question: what is the relative likelihood...") from `?rvnorm` Details. | `849e82e` |

## 4. Known limitations and notes (reviewed, not changed)

- **Overdetermined least-squares convention can put mass far from the variety.**
  `V(x^2 - 1, x - 1) = {1}`, yet the density at `x = -0.5` is `e^-0.78` of the peak: the
  two equations' Newton steps point in opposite directions and largely cancel. Sampler,
  rejection kernel, and `pdvnorm()` agree; this is a property of the convention, possibly
  worth a sentence in the paper.
- `project_onto_variety()` (straight-line homotopy) fails with "failed to meet
  `error_tol` at `dt_min`" for points about 1000x farther out than the variety's scale
  (all motion happens in a t-window of width ~1e-6 near 0); the fixed-step path returns a
  wrong point with a tolerance warning. Algorithmic, not a coding bug.
- Variable names that collide with generated Stan identifiers (`g`, `ndg`, `J`, `si`,
  `w`, `dg*`, `b*`) fail to compile on the generated path; on the template path the
  tibble output keeps the variable and drops the transformed `g`. A clear validation
  error would help.
- `geom_variety()` requires variables named x and y; other names give ggplot2's generic
  "Computation failed in `stat_variety()`" warning.
- Projection functions match matrix/data-frame columns by position (per `varorder`);
  `pdvnorm()` matches names. Possible enhancement.
- `rvnorm(rejection = TRUE, output = "stanfit")` errors via `match.arg()`;
  `remove_stan_files()` does not remove Windows `.exe` binaries.
- mpoly: `as.function()` also rejects a bare constant mpoly when given a `varorder`
  (worked around via `gradient_mpolylist()`).

## 5. Validation

All live runs used CmdStan 2.38.0. The `lp__` checks rely on `lp__` being exactly
`normal_lpdf(0 | g/ndg, si)` (or `multi_normal_lpdf`) for these unconstrained programs,
which is `pdvnorm(log = TRUE)`; draws were written with `sig_figs = 18` (12 for V1).
"After" rows ran against an installed copy of the final tree (template binaries in use).

| Check | Before | After |
|---|---|---|
| V1 `mp(c("y^2 + x^2 - 1", "y - x"))`, `Sigma = diag(c(.1, .01)^2)`: max \|`lp__` - `pdvnorm`\| | 670 | 2.8e-10 |
| V1 spreads (sd x, sd y) near the solution, requested (0.1, 0.01) | (0.011, 0.100) | (0.100, 0.011) |
| V2 `x^4 + y^4 - 1.23456789012345` (generated code), `sd = 1e-4`: max \|diff\| | 5e-12 (both rounded) | 2.8e-12 (both exact) |
| Template `2_2_vn`, `x^2 + y^2 - 1.23456789012345`, `sd = 1e-3`: max \|diff\| | — | 3.5e-13 |
| Circle radial distribution of Stan draws vs closed form (KS, thinned) | — | p = 0.945 |
| Rejection sampler, circle radial (KS, n = 4000) | — | p = 0.91 |
| Rejection sampler, sphere `sd = .002` radial (KS) | — | p = 0.84 |
| Rejection sampler, `mp("x")`, `sd = 1e-3` (0.13% acceptance), 20 seeds x 3000 draws vs N(0, 1e-3) | — | 0/20 KS rejections at 5% |
| Rejection kernel vs `log pdvnorm()`, `V(x^2 - 1, x - 1)`, 6 points | — | constant offset (exactly proportional) |
| `variety_solve(mp("t^2 - 2 t + 1"))` (template path) | error: missing `t` | `t = 1` |
| `rvnorm(1000, ..., thin = 5)` rows | 200 | 1000 |

## 6. Verification

- `devtools::test()`: **570 passed, 0 failed, 0 warnings, 1 skipped** (the env-gated
  install integration test), up from 501. New regression tests accompany every fix;
  CmdStan-gated ones check `lp__ == pdvnorm(log = TRUE)` and stanc syntax.
- CRAN-style `R CMD check` (`check-local.R` settings, tests and examples run):
  **Status: OK** (0 errors, 0 warnings, 0 notes).
- Every bug was reproduced at runtime before its fix; the codegen-order regression test
  was also run against the pre-fix code and failed there.
