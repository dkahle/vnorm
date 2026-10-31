# vnorm — Pre-compiled Templates by Default: Findings Log

**Purpose:** Record of David's question whether `rvnorm()` uses the shipped Stan
templates by default ("where a pre-compiled model is on file, use it"): the
check that answered it, the gaps it found, and their fixes.

**Repository:** `dkahle/vnorm`, branch `main` (local; not pushed).
**Toolchain:** R 4.6.1 (arm64 macOS), CmdStan 2.38.0 (cmdstanr 0.9.0, instantiate
0.2.3), mpoly 1.1.2.900.

---

## 1. The check

The default works. On David's installed copy of `d634a37`:

- All 36 template binaries were in `bin/stan`, built from template sources
  identical to `src/stan`.
- With default arguments, every catalog shape ran its own binary, in 0.4–1.0 s:
  1–3 variables, total degree 1–3, `homo` TRUE and FALSE, with and without a
  scalar `w`. So did other variable names (`a^2 + b^2 - 1`, `y^2 + z^2 - 1`)
  and sparse polynomials (`x^3 - y`).
- The README's single-polynomial calls ran `2_2_vn`, `3_2_vn` and `2_2_vn_w`,
  and a `homo = FALSE` call ran `2_2_hvn`, in 1.2–1.8 s for 2,000 draws.
- Requests the catalog does not cover compiled, at 5–10 s each: degree 4, a
  system of two equations.

**Method.** A fit's executable is `fit$runset$exe_file()`. For calls that return
draws, the check wrapped the two constructors `rvnorm()` calls:
`instantiate::stan_package_model()` (via `assignInNamespace()`), and the imported
`cmdstan_model()` (in `parent.env(asNamespace("vnorm"))`). `trace()` on the
imported function hit the wrong binding and logged template runs as compiles.

## 2. Gaps and fixes

| # | Gap | Commit |
|---|---|---|
| 1 | `w = c(-a, a)`, or a list giving every variable `c(-a, a)`, compiled (4.8 s), though it is the box of `w = a`, which the `_w` templates take | `9a454aa` |
| 2 | The templates ignored `vars`: with `vars = c("x", "y", "z"), w = 2`, the default sampled x and y, and `pre_compiled = FALSE` sampled z too | `5fab350` |
| 3 | When a template binary was missing, `rvnorm()` compiled without saying so | `2a82502` |
| 4 | `mpoly::totaldeg()` errors on an `mpolyList` of length 1 (its list branch needs `length(x) > 1`), so `rvnorm(n, mpolyList(p), sd)` fails on every path | mpoly (David); then vnorm's side, `03a0bde` (§6) |

### 2.1 Equal symmetric boxes — `9a454aa`

`rvnorm_normalize_window()` returns the same box `(-a, a)` on every variable as
the scalar `a`. Boxes that differ between variables, are asymmetric, or leave a
variable unbounded still use generated code. The `pre_compiled` documentation now
says which requests the templates cover.

### 2.2 `vars` on every path — `5fab350`

- **Meaning.** `vars` names the variables to sample. It must include the
  polynomial's. Any others are coordinates the polynomial does not involve, such
  as z for the cylinder x² + y² = 1 in three dimensions, and the distribution is
  flat along them.
- **Before.** Only the generated single-polynomial program honored it.
  - The templates ignored it.
  - `"simple"` output kept only the polynomial's variables, on every path.
  - The programs for systems, `code_only`, and `rejection = TRUE` ignored it.
  - For a polynomial in one variable, the generated program's gradient norm failed
    to build once `vars` added a variable: it took `Reduce()` over an mpoly's
    terms.
- **Now.**
  - Every path samples `sort(vars)`.
  - The templates take an extra variable as a template variable whose terms have
    coefficient 0, so `x^2 + y^2 - 1` in (x, y, z) runs `3_2_vn_w`. The x/y/z
    mapping covers extra variables.
  - Generated programs give extra variables parameters, but no Jacobian column or
    gradient term.
  - `user_compiled = TRUE` gives way to the other paths, since a
    `compile_stan_code()` model samples only the polynomial's variables. With
    `verbose = TRUE`, it says so.
- **Decisions I made**, within David's "honor `vars` everywhere":
  - A `vars` that omits a variable of `poly` is an error.
  - An extra variable without a window is an error ("bound z with `w`"). The
    distribution is improper along it, and HMC cannot converge.
  - The template path now maps variables to x/y/z in sorted order, not in the
    order they appear in the terms. Output names are restored either way. Only
    the variable that a template's `dgx`, `dgy` columns refer to can change.

### 2.3 Missing binaries — `2a82502`

- An installation can hold the templates' `.stan` files without their binaries.
  Installs under R CMD check skip them by design, as do
  `VNORM_PRECOMPILE_STAN=false` and a failed compile.
- **Before.** `rvnorm()` compiled every model without saying so. Its "falling
  back" message fired only when a `.stan` file was missing, which the install hook
  never causes, and then on every call.
- **Now.** One message per session, in both cases and whatever `verbose` is: the
  models are not installed, and reinstalling vnorm with CmdStan available builds
  them.

## 3. Not changed

- **README.** Its example `rvnorm(2000, mp("x^2 + y^2 + z^2 - 1"), sd = 0.1,
  pre_compiled = TRUE)` reads as if the templates must be switched on. I left it
  to David, who has `README.Rmd` open with uncommitted edits.
- **Length-1 `mpolyList`.** Finding 4 was mpoly's first; §6 covers the rest.
- **Model quantities in tibble output.** Template draws carry `g`, `dgx`, `dgy`,
  `ndg`; generated draws carry `g`, `ndg`. As before.

## 4. Validation

Installed copies of `2a82502`, in fresh R processes. The extra-variable runs use
`x^2 + y^2 - 1` with `vars = c("x", "y", "z")`, `w = 1.5`, `sd = 0.1` and 4,000
draws. Their `lp__` should be `pdvnorm()` of (x, y) plus the log Jacobian of
Stan's three box transforms.

| Check | Result |
|---|---|
| `w = 1.5`, `c(-1.5, 1.5)`, `list(x = , y = )` | all run `2_2_vn_w` (0.7–1.2 s) |
| Extra variable, default path | runs `3_2_vn_w`. Max \|`lp__` − expected\| = 3.7e-13. KS of z vs U(−1.5, 1.5): p = 0.92. KS of the radius vs the rejection sampler's iid draws: p = 0.81. ESS of z: 1,715 |
| Extra variable, `pre_compiled = FALSE` | generated code; 1.7e-13; p = 0.28; p = 0.54; ESS 1,620 |
| Extra variable, `rejection = TRUE` | columns x, y, z; KS of z: p = 0.37 |
| `"simple"` output with `vars = c("z", "y", "x")` | columns x, y, z |
| System `(x^2 + y^2 - 1, x - y)`, `sd = 0.05`, same `vars` and `w`, 2,000 draws | max \|`lp__` − expected\| = 9.9e-14; KS of z: p = 0.36 |
| The catalog check (§1), rerun | 36 of 36 shapes run their own binary (≤ 0.7 s); the README calls run templates; `w = c(-2, 2)` runs `2_2_vn_w` (0.4 s; was 4.8 s) |
| Installed with `VNORM_PRECOMPILE_STAN=false`, two calls | one message; both sample |
| `vars = "x"`; `vars` with z and no `w` | "it lacks y"; "bound z with `w`" |

## 5. Verification

- `devtools::test()` on the final tree: **763 passed, 0 failed, 0 warnings, 1
  skipped** (the env-gated install test). That is up from 718.
- Each fix's new tests were run against the pre-fix copies of its R files and failed
  there: 2 tests for fix 1, 7 for fix 2, and 2 for fix 3. One of fix 3's fails
  there only on its new message text.
- CRAN-style `R CMD check` (`check-local.R` settings, tests and examples run) on
  the final tree: **0 errors, 0 warnings, 1 NOTE**.
  - The NOTE is "Non-standard file/directory found at top level: `README_files`".
    That directory comes from knitting `README.Rmd`; it is untracked, so a clean
    checkout does not have it. `.Rbuildignore` lists `README.html` and
    `README_cache` but not `README_files`.
  - testthat summary: `[ FAIL 0 | WARN 2 | SKIP 3 | PASS 761 ]`. The 2 warnings
    and 3 skips are the same as in the third Codex follow-up.

## 6. Follow-up: an `mpolyList` of one polynomial — `03a0bde`

- **mpoly's fix.** The installed mpoly (1.1.2.900, rebuilt 2026-10-31 14:58 UTC)
  takes the list branch of `totaldeg()` for any list that is not an mpoly, so
  `totaldeg(mpolyList(p))` is 2.
- **vnorm still failed.** With that mpoly, `rvnorm()` failed on
  `mpolyList(mp("x^2 + y^2 - 1"))` on every Stan path:
  - the templates' helpers take an mpoly ("`poly` must be an mpoly object"), and
    `swap()` failed when variables needed renaming;
  - the generated single-polynomial program, and `code_only`, called `mp()` on the
    list.
  The rejection path worked. `pdvnorm()` already treats such a list as its
  polynomial (`R/pdvnorm.R:316`).
- **Why that is right.** For one equation, the systems program's normalized
  residual J⁺g = g ∇g / |∇g|² has norm |g| / |∇g|, the single-polynomial residual.
- **Fix.** `rvnorm()` unwraps the list before it chooses a path, and
  `compile_stan_code()` does too.
  - The templates now serve such lists.
  - A vector `sd` gets the single-polynomial error, as in `pdvnorm()`.
  - `compile_stan_code()` caches a model for the list under the name that
    `rvnorm(user_compiled = TRUE)` looks up. For shapes the templates cover, it
    now asks for `custom_stan_code = TRUE`, as for the polynomial.
- **Other functions.** `variety_solve()`, `geom_variety()` and
  `rejection_sampler()` gave the same results for `p` and `mpolyList(p)`.
  `project_onto_variety()` and `project_onto_variety_newton()` are documented for
  an mpoly. They fail on every `mpolyList`, systems included, with "argument
  "namevec" is missing", so I did not change them.

Validation on an installed copy of `03a0bde` with the new mpoly, in a fresh R
process:

| Check | Result |
|---|---|
| `mpolyList(p)`, defaults; `w = 1.5`; `homo = FALSE`; variables `a, b` | `2_2_vn`; `2_2_vn_w`; `2_2_hvn`; `2_2_vn` |
| `lp__` − `pdvnorm(log = TRUE)` of the list, 2,000 draws, template and generated code | within ±7.5e-15 |
| `rejection = TRUE`, same seed | draws identical to the polynomial's |
| `compile_stan_code(mpolyList(mp("x^4 + y^4 - 1")))`, then `user_compiled = TRUE` | runs the cached model's executable |
| `compile_stan_code(mpolyList(p))`; `sd = c(0.1, 0.2)` | the polynomial's errors |

Verification:

- `devtools::test()` with the new mpoly: **771 passed, 0 failed, 0 warnings, 1
  skipped**.
- Both new tests failed against the pre-fix copies of `R/rvnorm.R` and
  `R/compile_stan_code.R`.
- CRAN-style `R CMD check` on `03a0bde`: **0 errors, 0 warnings, 1 NOTE**, the
  same `README_files` NOTE as in §5. testthat summary:
  `[ FAIL 0 | WARN 2 | SKIP 3 | PASS 769 ]`.
