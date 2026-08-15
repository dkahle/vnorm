# vnorm — Fix Specification

**Purpose:** This document specifies defects and improvements identified in a code review of the
`vnorm` R package (variety normal distribution: Stan-based sampling, pseudo-density evaluation,
variety projection, and ggplot2 visualization). It is written as input for an implementation
plan: each item states the problem, evidence, affected files, the required fix, and acceptance
criteria.

**Repository:** `dkahle/vnorm`, branch `main`.
**Package conventions:** tidyverse-style ergonomics without a hard tidyverse dependency; coding
style should follow `mpoly` (https://github.com/dkahle/mpoly). Polynomials are `mpoly` /
`mpolyList` objects. Sampling backend is CmdStan via `cmdstanr` + `instantiate` precompiled
templates. Tests use testthat edition 3; CmdStan-dependent tests must skip gracefully when
CmdStan is unavailable (see helpers at the top of `tests/testthat/test-rvnorm.R`).

**General constraints for all fixes:**

- Do not change exported function signatures except where an item explicitly says to.
- Every behavioral fix gets a regression test. Tests that need CmdStan must use the existing
  `skip_if_no_cmdstan()` / fake-model patterns already in the test suite.
- Run `devtools::document()` after roxygen changes; `man/` is generated, never hand-edited.
- Validation errors should use `stop(..., call. = FALSE)` with backtick-quoted argument names,
  matching the existing style.

---

## Priority 1 — Silent correctness bugs

### 1.1 Precompiled path samples the wrong distribution when monomial variable order is non-canonical

**Problem.** `make_coefficients_data()` (`R/helper-rvnorm.R`) matches user polynomial
coefficients to precompiled Stan template data slots by *printed monomial name*. The printed
name depends on the order variables appear within each term. `mp("y x + x^2 - 1")` produces the
coefficient name `"y x"` → `byx`, which does not match the template slot `bxy`. The assignment
`required_coefs[names(available_coef)] <- available_coef` then *appends* a spurious `byx` entry
and leaves `bxy = 0`. CmdStan ignores unknown data fields, so with the default
`pre_compiled = TRUE`, `rvnorm(n, mp("y x + x^2 - 1"), sd = .1)` silently samples from
`x^2 - 1 = 0` instead of the intended variety.

**Verified reproduction:**

```r
vnorm:::make_coefficients_data(mp("y x + x^2 - 1"), 2, 2)
#> $bxy: 0   (template slot never filled)
#> $byx: 1   (spurious field, silently ignored by CmdStan)
```

**Same root cause elsewhere.** `get_custom_stan_code()` (`R/compile_stan_code.R`) reorders
monomials (`reorder(m, varorder = vars)`) before building coefficient names for the generated
Stan code, but `get_coefficients_data()` (`R/helper-rvnorm.R`) in its single-`mpoly` branch
takes `coef(poly)` names raw. So `compile_stan_code()` + `rvnorm(..., user_compiled = TRUE)`
can also produce code/data name mismatches (these at least fail loudly in Stan rather than
silently).

**Required fix.**

1. Create one internal helper (e.g. `stan_coef_names(poly, varorder)`) that reorders each
   monomial to a canonical variable order **before** converting to a Stan data name, and route
   *all* coefficient-name computations through it: `make_coefficients_data()`,
   `get_coefficients_data()` (both branches), `get_listed_coeficients()` call sites, and
   `get_custom_stan_code()`.
2. In `make_coefficients_data()`, after mapping, **fail loudly** if any user coefficient name
   does not match a template slot (i.e. any name in `available_coef` not in `stan_names`).
   This is a cheap invariant check that prevents any future silent mismatch.

**Acceptance criteria.**

- `rvnorm(n, mp("y x + x^2 - 1"), ...)` and `rvnorm(n, mp("x y + x^2 - 1"), ...)` produce
  draws from the same distribution (test via `make_coefficients_data()` equality; a sampling
  smoke test may be CmdStan-gated).
- `vnorm:::make_coefficients_data(mp("y x + x^2 - 1"), 2, 2)` returns `bxy = 1` and contains
  no `byx` entry.
- Unit test: for a polynomial written in scrambled variable order, the coefficient lists from
  both orderings are identical.
- Unit test: `get_custom_stan_code()` data-block names equal `get_coefficients_data()` names
  for the same polynomial, for both `mpoly` and `mpolyList` inputs.

### 1.2 `variety_solve()` rejects valid single-polynomial input

**Problem.** `R/variety_solve.R` computes `n_eqns <- length(polylist)`, but `length()` of an
`mpoly` counts **terms**, not equations. `variety_solve(mp("x - 1"))` errors with
"The variety is not zero-dimensional." even though it is one equation in one variable.

**Also in this function:**

- The `vars` argument is documented as "variable names to average" but the draw extraction
  uses `mpoly::vars(polylist)`; `vars` only feeds the dimension check. Honor `vars` in the
  column subset (validate that all requested names exist in the draws), or remove the argument.
- The docs/argument `stanfit` describe a "Stan fit"; the object is a `CmdStanMCMC`. Fix the
  wording (do not rename the argument).

**Required fix.** `n_eqns <- if (is.mpoly(polylist)) 1L else length(polylist)`. Honor `vars`
in the extraction. Update roxygen docs.

**Acceptance criteria.**

- `variety_solve(mp("x - 1"), n = <small>)` runs without the dimension error (CmdStan-gated
  test, or assert only that the guard passes by mocking `rvnorm`).
- A multi-term single `mpoly` in 1 variable passes the guard; a genuine non-square system
  (e.g. 1 polynomial in 2 variables as an `mpolyList` of length 1) still errors.

### 1.3 Wrong objective in `project_onto_variety_lagrange()` optim branch

**Problem.** `R/project_onto_variety.R`, non-`"newton"` branch: the objective is
`function(.) sum(dLf(.))^2` — the **square of the sum** of gradient components — instead of the
sum of squares. Gradient components can cancel, so `optim` can converge to a non-stationary
point. The `la0` initialization two lines above and the Newton branch both correctly use
`sum(dLf(...)^2)`.

**Required fix.** Change the objective to `function(.) sum(dLf(.)^2)`.

**Acceptance criteria.** Test that `project_onto_variety_lagrange(x0, poly, method = "BFGS")`
(or another `optim` method) returns a point whose residual `|poly(x)|` is below tolerance for a
case where the two objectives demonstrably differ (e.g. a starting point where gradient
components cancel under the buggy objective).

### 1.4 Users cannot override `adapt_delta` / `max_treedepth` in `rvnorm()`

**Problem.** `R/rvnorm.R` builds `sample_args <- list(..., adapt_delta = .999,
max_treedepth = 20L, ..., <dots>)`. If a user passes `rvnorm(..., adapt_delta = .9)`, the list
holds two `adapt_delta` entries and `do.call(model$sample, sample_args)` fails with
"formal argument matched by multiple actual arguments" (verified).

**Required fix.** Either promote `adapt_delta` and `max_treedepth` to formal arguments of
`rvnorm()` with the current defaults (preferred — discoverable, documented), or merge `...`
over the defaults so user values win (e.g. build a defaults list and `utils::modifyList()` it
with the dots). Document that these defaults are conservative and why.

**Acceptance criteria.**

- `rvnorm(..., adapt_delta = .9)` does not error (testable with the existing
  `fake_cmdstan_model()` capture pattern in `tests/testthat/test-rvnorm.R` — assert the
  captured `sample()` args contain the user value, not the default).
- Default behavior unchanged when the user passes nothing.

---

## Priority 2 — User-visible defects

### 2.1 Fixed-step homotopy in `project_onto_variety()` can stop short of `t = 0`

**Problem.** In the `adaptive = FALSE` path, `ts <- seq(t_start, t_end, -dt)` does not include
`t_end = 0` when `dt` does not divide 1 (e.g. `dt = 0.03` ends at `t = 0.01`). The final Newton
corrections then solve the homotopy at `t ≈ 0.01`, not the target variety, and the residual
warning may fire spuriously.

**Required fix.** After building `ts`, append `t_end` if `ts[length(ts)] != t_end`. Also:

- Remove the dead `ts <- seq(1, 0, -dt)` computed unconditionally before the adaptive branch.
- In the adaptive branch, advance with the Heun point `v_heun` (already computed for error
  control) instead of the lower-order Euler point `v_euler` — it is free and more accurate.

**Acceptance criteria.** `project_onto_variety(c(1, 1), mp("x^2 + y^2 - 1"),
adaptive = FALSE, dt = 0.03)` returns a point with `|g(x)| < tol` and no tolerance warning.

### 2.2 `rejection_sampler()` can loop forever

**Problem.** The `while (n_remaining > 0)` acceptance loop has no iteration bound. If the
variety does not intersect the proposal box (a plausible user mistake in `w`), the function
hangs with no feedback.

**Required fix.** Add a cap (e.g. a `max_batches` or total-proposals argument with a generous
default). On hitting the cap, `stop()` with an informative message stating how many draws were
accepted and suggesting enlarging `w` or increasing `sd`. Optionally warn earlier if the first
few batches accept nothing.

**Acceptance criteria.** A test with a box that misses the variety (e.g.
`p = mp("x^2 + y^2 - 1")`, `w = list(x = c(10, 11), y = c(10, 11))`) errors promptly with the
informative message instead of hanging. Normal sampling behavior unchanged.

### 2.3 `pdvnorm()` normalization inconsistent between `mpoly` and length-1 `mpolyList`

**Problem.** For the same polynomial, `pdvnorm(x, mp("..."), sd)` and
`pdvnorm(x, mpolyList-of-1, sd)` return values differing by a constant factor (verified:
0.402 vs 1.604 on the circle at the same point; the shapes agree). Single-poly path uses a 1-D
normal constant; list path uses an n-D constant. Users comparing densities across calls (HDR
plots, importance weights) get inconsistent results.

**Required fix.** Pick one convention and make both paths agree. Recommended: make a length-1
`mpolyList` delegate to the single-`mpoly` path (unwrap), so equal inputs give equal outputs.
Document the chosen normalization convention in the roxygen Details for `pdvnorm()`.

**Acceptance criteria.** Test: `pdvnorm(x, p, sd)` equals
`pdvnorm(x, structure(list(p), class = "mpolyList"), sd)` for several points, both `homo`
settings where applicable, `log = TRUE` and `FALSE`.

### 2.4 Named-list `w` breaks the precompiled path in `rvnorm()`

**Problem.** `rvnorm()` selects a `_w` (windowed) precompiled template whenever `w` is
supplied, but shipped windowed templates declare scalar `real w`. A named-list `w` flows into
`c(stan_data, "w" = w)`, producing `w.x`, `w.y` data fields; CmdStan fails with a confusing
"variable w not found" error.

**Required fix.** Validate `w` early in `rvnorm()` (reuse/adapt `rejection_normalize_window()`
semantics: scalar, length-2, or named list). When `w` is not a positive scalar and
`pre_compiled`/`user_compiled` would be used, fall back to the code-generation path (which
already bakes per-variable bounds into the Stan source) — with a `message()` explaining the
fallback when `verbose = TRUE`.

**Acceptance criteria.** `rvnorm(n, p, sd, w = list(x = c(-1.5, 1.5)))` works via the
generated-code path (CmdStan-gated), and a unit test confirms the precompiled branch is not
selected for non-scalar `w`. Scalar `w` keeps using the `_w` template.

### 2.5 `rvnorm()` return-value docs wrong for `output = "simple"`

**Problem.** Roxygen `@return` promises "a matrix whose rows are the individual draws"; the
code returns a `data.frame` (the package's own tests assert `is.data.frame`).

**Required fix.** Fix the documentation to say data frame (keep behavior; a data frame is the
better tidyverse-style default). While in the docs: the `@return` also says `stanfit`; the
object is `CmdStanMCMC` — correct the wording.

**Acceptance criteria.** `man/rvnorm.Rd` regenerated; description matches actual classes.

### 2.6 `rejection_sampler()` documented but not exported

**Problem.** `R/rejection_sampler.R` has a full roxygen page (generating
`man/rejection_sampler.Rd` with a `\usage` block) but no `@export`; the pkgdown site advertises
it while users cannot call it. R CMD check flags usage entries for unexported objects.

**Required fix (decision needed, pick one):**

- **Option A (recommended):** `@export rejection_sampler`. It is user-facing quality already:
  full docs, examples, thorough input validation, 251 lines of tests.
- **Option B:** mark internal (`@noRd` or `@keywords internal` without usage) and route users
  exclusively through `rvnorm(rejection = TRUE)`; update pkgdown reference index and README
  mentions accordingly.

**Acceptance criteria.** R CMD check passes the code/documentation-mismatch check; pkgdown
reference index consistent with the choice.

---

## Priority 3 — Ergonomics and polish

### 3.1 `geom_variety()` scale collisions

`geom_variety()` returns `list(layer, scale_linetype_discrete(...))`, so (a) two
`geom_variety()` layers in one plot trigger "Scale for linetype is already present" and the
second silently replaces the first, and (b) a user's own linetype scale added before the layer
is clobbered. Mitigate: only add the scale when it is needed (e.g. an `mpolyList` with a
legend), and document the interaction in `@details`. Full fix may not be worth complexity —
minimum bar is documentation plus no warning for the common single-polynomial case.

### 3.2 Example hygiene

- `pdvnorm()` `\dontrun` examples call `library("ggfunction")`, which was removed from
  `Suggests` and may not be installable — replace those examples with `ggplot2`-only
  equivalents or restore the dependency decision explicitly.
- `rvnorm()` examples use `lissajous()`, `plug()`, `swap()` (mpoly exports) without
  `library("mpoly")`; only `mp` is re-exported. Either add `library("mpoly")` to the examples
  or re-export the handful of mpoly verbs the examples rely on (`lissajous`, `plug`, `swap`,
  `vars`) in `R/reexport-mpoly.R`.

### 3.3 `warmup` default scales with total `n`, not per-chain draws

`warmup = max(500, floor(n / 2))` while `iter_sampling = ceiling(n / chains)`: at
`n = 1e5, chains = 8`, that is 50,000 warmup iterations per chain for 12,500 sampling
iterations. Change default to scale with per-chain draws, e.g.
`max(500, ceiling(n / chains / 2))`, and document.

### 3.4 Non-reproducible default in `project_onto_variety()`

`al = rnorm(2)` makes results vary run-to-run by default. Either keep and document loudly
(patch randomization is standard in parameter continuation) with a note to `set.seed()`, or use
a fixed default patch and expose randomization as opt-in.

### 3.5 Internal style cleanups (no behavior change)

- Replace `get("print.mpoly", asNamespace("mpoly"))` and `mpoly:::coef.mpoly`
  (`R/helpers.R`) with plain generic calls — `print()`, `coef()`, `deriv()` all dispatch via
  exported S3 generics. Then remove `"coef"`, `"deriv"`, `"reorder"` from
  `utils::globalVariables()` in `R/vnorm-package.R` if no longer needed.
- `R/compile_stan_code.R`: local `vars <- vector("list", ...)` shadows `mpoly::vars()` (works
  only because R skips non-function bindings in call position). Rename to `vars_by_poly`.
- Rename internal `get_listed_coeficients` → `get_listed_coefficients` (typo) — internal only,
  update all call sites (`R/helpers.R`, `R/helper-rvnorm.R`, `R/compile_stan_code.R`,
  `R/make_stan_files.R`, `R/make_multivariety_stan_codes.R`, tests). May be subsumed by the
  shared helper in item 1.1.
- Shipped Stan templates declare `real si`; generated code declares `real<lower=0> si`.
  Regenerate templates (via `make_stan_files()`) for consistency.
- Deduplicate the "Zero contours were generated" message (emitted in both
  `StatVariety$compute_group` and `iso_to_path()`) and reword to something user-actionable,
  e.g. "no zero contour found in the plotting window".
- `.onLoad` sets `options(cmdstanr_no_ver_check = TRUE)` unconditionally; only set when the
  option is unset, to respect user configuration.
- Fix overstated comment in `R/helpers-ggplot2.R` (`break_contour_ties`): the perturbation
  guarantees a deterministic nonzero sign, not distinct signs at distinct grid points.

---

## Suggested sequencing

1. **1.1** (shared coefficient-name helper + strict match check) — highest impact, touches the
   most files; land first so later items rebase onto it.
2. **1.2, 1.3, 1.4** — small, independent, high-value; can be one PR each or one combined PR.
3. **2.1–2.6** — independent; 2.6 requires the export/internal decision.
4. **3.x** — batch as one cleanup PR; 3.5's rename overlaps 1.1, do after it.

## Verification for the whole effort

- `devtools::test()` green locally with and without CmdStan installed (skip paths exercised).
- `R CMD check` (see `check-local.R`) with no new warnings/notes; the pre-existing
  documentation-mismatch warning for `rejection_sampler` resolved by 2.6.
- `devtools::document()` output committed; no hand edits under `man/`.
- README examples re-knit only if user-facing behavior they show has changed.
