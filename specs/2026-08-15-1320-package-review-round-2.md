# vnorm — Package Review Round 2: Findings Log

**Purpose:** Record of the 2026-08-15 full package review: what was checked, what was
found, and what was done. Unlike the July 5 specs (which were pre-approval fix
specifications), this is a log of work already landed; each item names its commit.

**Repository:** `dkahle/vnorm`, branch `main`.
**Review scope:** line-by-line code review, full test suite + CRAN-style R CMD check,
live CmdStan sampling smoke tests, and numerical validation of sampler/density/projection
agreement. Toolchain: R 4.6.1, CmdStan 2.38.0 (cmdstanr 0.9.0), mpoly 1.1.2.900.

---

## 1. Status of the July 5 specs at review start

- **Spec A** (`2026-07-05-1232-package-review-fixes.md`): essentially fully implemented
  by earlier commits. Remaining loose ends closed this round (§3).
- **Spec B** (`2026-07-05-1440-pdvnorm-varorder.md`): **entirely unimplemented** at
  review start. The uncommitted working-tree edit to the `pdvnorm` example
  demonstrated the bug rather than fixing it. Implemented in full this round (§2.1).

## 2. Correctness fixes landed

### 2.1 `pdvnorm` sorted variable order (Spec B) — commit `8cc3d77`

Both evaluation paths now pin `sort(mpoly::vars(poly))`: `as.function(...,
varorder = vs, silent = TRUE)` everywhere, `deriv(poly, var = vs)` instead of
`mpoly::gradient()`, name-based matching of `x` entries/columns (with errors on
partial/extra names; names disjoint from the variables stay positional so
`as.data.frame(matrix)` defaults keep working), and documented conventions.
The `deriv(var = vs)` change also fixes a second latent bug: for polynomials like
`y^2 + x`, the gradient `(1, 2y)` has a smaller variable set than the polynomial,
and mpoly's generated function then mis-evaluated points (returned a matrix), so
`homo = TRUE` densities were wrong. Regression tests: term-order invariance (incl.
`Sigma`-matrix branch and `mpolyList` path), constant-gradient ordering,
dropped-variable gradients, name matching/validation, and the plot-symmetry guard.

### 2.2 Zero density at off-variety critical points — commit `0f36144`

`pdvnorm(..., homo = TRUE)` returned `Inf` wherever `∇g = 0`. Correct limit: the
normalized residual `g/||∇g||` diverges unless the point is on the variety, so the
density is 0 (log `-Inf`) at off-variety critical points; `Inf` remains only at
singular points on the variety (e.g. the origin for `V(xy)`). Found by numerical
validation (§5): `integrate(pdvnorm(...), ...)` hit the spurious `Inf` at `x = 0`
for `x^2 - 1`. Note: in the `mpolyList` homo path, rank-deficient Jacobians go
through the SVD pseudoinverse, which drops residual components outside the row
space — this matches the intentional overdetermined-system convention and was
left unchanged.

### 2.3 Sorted variable order in all projection functions — commit `c03778c`

`project_onto_variety()` already used sorted order, but the Lagrange,
gradient-descent, and Newton variants defaulted to term-appearance order
(`mpoly::vars(poly)`), so for e.g. `mp("y^2 - x^3 - x")` (appearance order `y, x`)
the same `x0` meant different points across the family, and results returned in
different coordinate orders. Unified on sorted order; tolerance warnings now print
the residual in scientific notation. Found when a projection-residual sweep showed
"non-convergence" that was actually coordinate swapping.

### 2.4 rvnorm output column order unified — commit pending this round

Stan paths returned draw columns in term-appearance order while the rejection path
sorted them. All `"simple"`/`"tibble"` output now uses sorted variable order
(user-approved; matches the package-wide convention). `"tibble"` output keeps
sampler diagnostics (`.chain`, `.iteration`, `.draw`, `lp__`) after the variable
columns, and `lp__` is relocated **by name** rather than by position
(commit `48d242a`), removing the fragile `cbind(df[, -1], df[, 1])`.

### 2.5 Degenerate-input validation — commit `48d242a`

Constant (variable-free) polynomials previously died with cryptic errors
("argument is of length zero", "attempt to select less than one element").
`rvnorm()`, `pdvnorm()`, `rejection_sampler()`, and all four projection functions
now stop with "`poly` must contain at least one variable." `stat_variety()` errors
at construction when `poly` is missing (matching `geom_variety()`) instead of a
ggplot2 compute-time warning. `check_and_replace_vars()`'s >3-variable error gained
`call. = FALSE`.

## 3. Spec A loose ends closed — commit `7c2148d`

- `pdvnorm` examples no longer use the unavailable `ggfunction` package; replaced
  with small ggplot2-only helpers (`pdf_curve`, `pdf_raster`) defined inside the
  example.
- The misleading "Zero contours were generated" message (emitted when NO contour is
  found) was reworded to "No zero contour found in the plotting window." and
  deduplicated behind a single helper (was in 3 places).
- Overstated tie-break comment in `helpers-ggplot2.R` corrected.
- `geom_variety()` mapping-override behavior (user `aes(linetype=)` etc. is
  intentionally overridden by the stat wiring) is now documented in the roxygen
  details (user-approved: keep behavior, document).

## 4. Cleanups — commit `48d242a` and this round

- Removed six dead `ptilde` closures and a redundant `missing(vars)` check in
  `rejection_sampler.R`; joined its split roxygen block.
- Dropped unused imports: `stats::rnorm`, `stats::coef`, `utils::head`,
  `ggplot2::ggplot`, `ggplot2::labs`; geom test files now attach ggplot2
  explicitly (they had leaned on the package's imports).
- `mean(is.constant(x)) == 1` → `all(is.constant(x))` in `pdvnorm.R`.
- Windowed single-`mpoly` codegen now emits `real w;` on its own line (matches the
  `mpolyList` branch and shipped templates); the test snapshot was updated.
- LICENSE trimmed to the 2-line CRAN form (full text stays in `LICENSE.md`) —
  commit `a3defb3`.

## 5. Numerical validation performed (scripts in session scratchpad)

All validation ran live CmdStan (2.38.0) against the current code.

| Check | Result |
|---|---|
| 1d `x^2 - 1`: rejection draws vs numerically normalized `pdvnorm` (KS, n = 4000) | KS 0.014, p = 0.43 — agree |
| 1d: Stan draws (thinned) vs density, mode-folded | KS 0.027, p = 0.44 — agree |
| 2d circle: θ uniformity / radial density vs closed form | p = 0.08 / p = 0.61 — agree |
| `pdvnorm` vs closed-form density on 4000 draws (circle) | exactly equal |
| Precompiled template `2_2_vn` vs generated code (circle) | KS p = 0.95/0.16 — agree |
| Hetero template `2_2_hvn` vs generated code | p = 0.18 — agree |
| Windowed template `2_2_vn_w` vs generated code (`xy - 1`, branch-conditional) | p = 0.19/0.12 — agree |
| Term-order invariance through precompiled path (`y x + x^2 - 1`) | coefficient data identical; sampled differences attributable to chain noise |
| `sd = .1` vs `Sigma = .01`, same seed | identical draws |
| Projection residuals (circle, cubic, 3d ellipsoid; 10 random starts each) | homotopy (adaptive + fixed): ≤ 3e-15; Lagrange/GD converge; Newton can stall at indefinite-Hessian starts but warns |

**Known behavior (documented, not a bug):** on varieties with several disconnected
components, HMC chains cannot cross between components, so per-component draw
shares reflect chain initialization (observed 3:1 on `x^2 - 1` with 4 chains; the
rejection sampler gives the correct ~1:1). Now documented in `rvnorm()` details
with the advice to raise `chains` or use `rejection = TRUE`.

## 6. Items reviewed and intentionally left alone

- `mpolyList` homo pseudoinverse convention at rank deficiency (§2.2 note).
- `geom_variety()` scale collision when two `mpolyList` layers share a plot
  (documented limitation; Spec A 3.1 minimum bar met previously).
- `make_multivariety_stan_codes.R` string patch `v[v == "**"] <- 0` — brittle but
  developer-only template generator, covered by tests.
- `dist = "unif"` acceptance in `rejection_sampler` for `mpolyList` compares a
  quadratic form against `sd` — reviewed; `rejection_validate_scale()` restricts
  the reachable combinations, and the single-`mpoly` path (the documented use) is
  dimensionally consistent.
- Top-level `alm.R`/`alm.Rd` (outside the package): untouched; not referenced by
  package or paper.

## 7. Verification

- `devtools::test()` final full run: **0 failed, 0 warnings, 500 passed,
  2 skipped** (the env-gated install-integration test and one `skip_on_cran`).
- `R CMD check` (CRAN-style via `check-local.R`, tests and examples enabled):
  clean after restoring the `stats::coef` import that the review's import trim
  had over-pruned (`get_derivative` passes `coef` as a bare function value,
  which call-site greps miss — the one review-introduced regression, caught by
  the check itself).
- README re-knit against the installed post-fix package (message rewording and
  regenerated first figure).
- Manuscript consistency review delivered separately:
  `../vnorm-paper-fix-spec.md` (repo top level, next to the PDF).
