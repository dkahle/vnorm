# vnorm — Fix Specification: `pdvnorm()` variable-order bug

**Purpose:** Input for an implementation plan. One defect, fully diagnosed and verified. Each
section states the problem, root cause, required fix, and acceptance criteria.

**Repository:** `dkahle/vnorm`, branch `main`.
**Conventions:** tidyverse-style ergonomics; coding style follows `mpoly`
(https://github.com/dkahle/mpoly). Tests are testthat edition 3. Validation errors use
`stop(..., call. = FALSE)` with backtick-quoted argument names. Run `devtools::document()`
after roxygen changes; never hand-edit `man/`.

---

## The defect

`pdvnorm()` evaluates the density with coordinates in the polynomial's **term-appearance
variable order**, not a deterministic order, and not the order of the user's input columns.
Two mathematically identical polynomials give different answers:

```r
p1 <- mp("-1 x^2 (x + 1) + y^2", varorder = c("x", "y"))   # vars(p1) is c("x", "y")
p2 <- mp("y^2 - x^2 (x + 1)",    varorder = c("x", "y"))   # vars(p2) is c("y", "x") (!)

pdvnorm(c(0.5, 0.8), p1, sd = .1)   # 2.1364...
pdvnorm(c(0.5, 0.8), p2, sd = .1)   # 0.1912...  <- WRONG: (x, y) evaluated as (y, x)
```

User-visible symptom: density surfaces / `geom_pdf_2d`-style plots come out reflected
("rotated, loop pointing down") depending on which term the user happened to write first.

### Root cause (three interacting facts)

1. **`mpoly::mp(varorder=)` only reorders variables *within each term*** (verified in `mp()`
   source); it does not establish an overall variable order. `mpoly::vars()` returns variables
   in **order of first appearance across terms**. So `mp("y^2 - ...", varorder = c("x","y"))`
   still has `vars()` equal to `c("y", "x")`. Do **not** try to fix this in mpoly; fix the
   consumer.
2. **`pdvnorm()` builds evaluation functions without a `varorder`**: `as.function(poly)` in the
   single-`mpoly` path and `as.function(poly)` for the `mpolyList` path (`R/pdvnorm.R`,
   currently around lines 365 and 487; search for `suppressMessages(as.function(`). mpoly then
   assumes `vars(poly)` order — appearance order.
3. **The diagnostic is suppressed.** mpoly emits `f(.) with . = (y, x)` when it guesses the
   order; `pdvnorm()` wraps every `as.function()` call in `suppressMessages()`, so the user
   never sees it.

### Package precedent

The rest of vnorm already pins a deterministic order — follow it:

- `rvnorm()` rejection path: `vars = sort(mpoly::vars(poly))`
- `project_onto_variety()`: default `varorder = sort(mpoly::vars(poly))`
- `geom_variety()` helpers: explicit `varorder = c("x", "y")`

`pdvnorm()` is the outlier.

---

## Required fix

All changes are in `R/pdvnorm.R` (plus docs/tests). Do not change any exported signature.

### 1. Pin sorted variable order in every evaluation function

Define once, early in `pdvnorm()` (and reuse in both paths):

```r
vs <- sort(mpoly::vars(poly))
```

Then:

- **Single-`mpoly` path:** `g_func <- as.function(poly, varorder = vs, silent = TRUE)`.
- **Gradient — subtle and load-bearing:** the current code uses
  `grad_g_obj <- mpoly::gradient(poly)`, whose *component order follows `vars(poly)`*
  (appearance order). Component `i` of the gradient must correspond to variable `vs[i]`.
  Replace with `stats::deriv(poly, var = vs)` (an `mpolyList` ordered by `vs`), and build the
  gradient function with `as.function(..., varorder = vs, silent = TRUE)`. The
  constant-gradient fast path (`unname(unlist(grad_g_obj))`) inherits the component order from
  the deriv call, so it is correct once the deriv is ordered by `vs` — but add a test for it
  (a linear polynomial written y-first).
- **`Sigma`-matrix single-poly branch:** the normal-direction vector lives in variable space;
  after the change its coordinates follow `vs`. Document that `Sigma` rows/columns correspond
  to the variables in sorted order.
- **`mpolyList` path:** replace `vars <- mpoly::vars(poly)` with `vars <- sort(mpoly::vars(poly))`;
  pass `varorder = vars` to the `as.function(poly, ...)` call for `g_fns` (it currently has
  none) — the per-polynomial `deriv(poly[[i]], var = vars)` / `as.function(..., varorder = vars)`
  calls then line up automatically. The `homo = TRUE` covariance (`q = n`) likewise follows
  sorted variable order; document it.
- Prefer passing `silent = TRUE` to `as.function()` over blanket `suppressMessages()` where
  mpoly supports it; keep `suppressMessages()` only if some call site still chats.

### 2. Match matrix/data-frame input by column name when names are present

Ergonomic half of the fix, applied in both paths where `x` is coerced to a matrix:

- If `x` has column names and they contain **all** of `vs` (as a set): reorder columns to `vs`
  before evaluation. Extra columns beyond `vs`: **error** (ambiguous input), message naming the
  unexpected columns.
- If `x` has column names that only partially overlap `vs`: **error**, message naming the
  missing variables. Do not silently fall back to positional.
- If `x` has no column names (or is a bare vector): interpret positionally **in sorted
  variable order** and proceed (this is the documented convention).
- Named numeric vectors: nice-to-have; if implemented, same rules as column names.

### 3. Documentation

- `@param x`: state the convention — "columns/entries are interpreted in sorted variable order
  (`sort(mpoly::vars(poly))`); if `x` has column names covering the variables, columns are
  matched by name instead."
- `@param Sigma` / `@details`: note that when `homo = TRUE` the covariance coordinates follow
  the sorted variable order.
- **Fix the motivating example** (`R/pdvnorm.R`, the "1d variety" example currently near lines
  116–132): keep only the natural form `mp("y^2 - x^2 (x + 1)")` — no `varorder` argument, no
  defensive term ordering — since after the fix it plots correctly. Delete the workaround line.

### 4. Out of scope

- No changes to mpoly.
- No changes to `rvnorm()`, `geom_variety()`, or projection functions (already consistent).
- The length-1 `mpolyList` → `mpoly` delegation added previously must keep working unchanged
  (it inherits the fix automatically via delegation).

---

## Acceptance criteria / tests

Add to `tests/testthat/test-pdvnorm.R`:

1. **Term-order invariance (the bug):**
   ```r
   p_xy <- mp("-1 x^2 (x + 1) + y^2")
   p_yx <- mp("y^2 - x^2 (x + 1)")
   x <- rbind(c(0.5, 0.8), c(-0.3, 0.2), c(1, 1))
   expect_equal(pdvnorm(x, p_yx, sd = .1), pdvnorm(x, p_xy, sd = .1))
   ```
   Also with `homo = FALSE`, `log = TRUE`, and scalar `Sigma`.
2. **Constant-gradient fast path ordering:** `p <- mp("2 y - 3 x")` (y appears first in vars);
   compare `pdvnorm` against a hand-computed value at a couple of points, `homo = TRUE`.
3. **`mpolyList` term-order invariance:** `mp(c("y x - 1", "x^2 + y^2 - 1"))` vs the same
   system written x-first; equal densities for `homo = TRUE` and `FALSE`.
4. **Name-based column matching:** for `p <- mp("x^2 + y^2 - 1")`,
   `pdvnorm(cbind(y = 0.8, x = 0.5), p, sd = .1)` equals
   `pdvnorm(cbind(x = 0.5, y = 0.8), p, sd = .1)` and equals the unnamed
   `pdvnorm(c(0.5, 0.8), p, sd = .1)`.
5. **Name validation:** columns `c("x", "z")` for a polynomial in `x, y` errors mentioning `y`;
   extra column `c("x", "y", "z")` errors mentioning `z`.
6. **Regression guard for the plot symptom:** evaluate `p_yx` on a small grid and check the
   density is symmetric under `y -> -y` and *not* under `x -> -x` (the curve
   `y^2 = x^2(x+1)` is symmetric about the x-axis; the buggy swapped evaluation is symmetric
   about the y-axis instead).
7. **Sigma-matrix branch:** single poly, 2x2 `Sigma`, `homo = TRUE`: term-order invariance as
   in (1).

Whole-effort verification:

- `devtools::test()` fully green (CmdStan-dependent tests skip cleanly as configured).
- `devtools::document()`; `man/pdvnorm.Rd` regenerated.
- The previously verified value `pdvnorm(c(0.5, 0.8), mp("-1 x^2 (x + 1) + y^2"), sd = .1)`
  must remain `2.136453` (the x-first form was already evaluated correctly; the fix must not
  change it).
