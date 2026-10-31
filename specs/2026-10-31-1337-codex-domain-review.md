# vnorm — Codex's Review of Domain Noise: Findings Log

**Purpose:** Record of the follow-up to Codex's review of `d634a37..d99e67e`, the
domain-noise feature (logged in `2026-10-31-1242-domain-noise.md`): each finding's
verdict, the fixes with their commits, and the validation.

**Repository:** `dkahle/vnorm`, branch `main` (local; not pushed).
**Toolchain:** R 4.6.1 (arm64 macOS), CmdStan 2.38.0 (cmdstanr 0.9.0), mpoly 1.1.2.900.

---

## 1. Verdicts

All six findings reproduced as reported, and all are fixed.

| # | Codex finding | Reproduction | Commit |
|---|---|---|---|
| 1 (P1) | Projection accepts points off the variety | (x² − 1, 0.1x − 0.1): 61 of 100 latent points at −0.995. (x, x − 1), an empty variety: every point at 0.5 | `539661f` |
| 2 (P1) | Repeated factors bias component weights | {0, 1} gets shares 0.49 for x² − x and 0.66 for x(x − 1)², at every `pilot_sd` | `f060950` |
| 3 (P2) | Scalar-window padding of free variables | `"x"`, `vars = c("x", "y", "z")`, `w = 0.01`: "Only 1 of 20 pilot draws projected" | `5c4f213` |
| 4 (P2) | Free coordinates change a system's classification | (x, x²) with `vars = c("x", "z")`: "density is 0 off the variety" on the rejection path, all projections failing on the HMC path. The first also hit codomain `rejection = TRUE` | `45073b0` |
| 5 (P2) | Top-up runs reuse chain and iteration labels | Clear from the code | `fee6f67` |
| 6 (P2) | Remapping's temporary names collide with variables | `mp("tmp2 + 2 a - 1")`: mpoly's `swap()` error | `b52f803` |

**Architecture.** Codex recommended splitting the latent sampler from the noise.
`3ef5e27` does that internally, with draws unchanged for a given seed:
- `rvnorm_domain_latent()` returns the points, their pilot draws' chain and
  iteration, and projection counts;
- `rvnorm_gaussian_noise()` adds N(0, Σ).

The public parts of the recommendation are open questions for David (§4).

## 2. Fixes

### 2.1 Feasibility — `539661f`

- **Cause.** A vanishing Gauss–Newton step J⁺g does not put a point on the
  variety. For an overdetermined system, it also vanishes at least-squares points
  where g does not; for an empty variety, only there.
- **Fix.** `variety_feasible()` requires every equation to vanish at the projected
  point, to within a distance |g_i| ≤ 1e-8 (1 + |u|) |∇g_i|, or to rounding
  error in its terms. That covers singular points, where ∇g_i vanishes too.

### 2.2 Repeated factors — `f060950`

- **Cause.** Where `poly` is square-free, the pilot's residual |J⁺g| is the
  distance to the variety to first order. On a component where `poly` has a
  factor of multiplicity k, it is about distance / k. The pilot's tube there is
  k times wider, so the component gets k times its share. A smaller `pilot_sd`
  does not change this.
- **Fix: detect, don't correct.**
  - `project_rows_closest()` now returns each draw's residual, and
    `rvnorm_domain_latent()` compares it with the distance to the projection.
  - The ratio distance / residual is about 1 on smooth, square-free parts, about
    k on a repeated component, and at most √2 at a crossing of two smooth
    branches.
  - Above 1.5 at more than 1% of draws, `rvnorm()` stops with a message naming
    repeated factors, or a `pilot_sd` large next to the curvature.
  - A repeated equation in a system trips it too.
  - Gauss–Newton converges only linearly onto a repeated root, so the projection
    now allows 100 iterations.
- **Docs.** The `pilot_sd` and Details text now state the conditions for uniform
  points: a variety of finite size (compact, or bounded by `w`), equations
  without repeated factors, and an O(pilot_sd²) error where the variety is
  smooth. Isolated singular points do not affect it.

### 2.3 The others

- **`5c4f213`.** Free variables keep their window bounds. A scalar window becomes
  per-variable bounds when there are free variables, at the cost of the
  templates.
- **`45073b0`.** Whether a system is overdetermined, its dependence check, and the
  projection's rank requirement now count only the variables the system involves.
  This fixes the rejection sampler for codomain noise too.
- **`fee6f67`.** Each extra pilot run's chains continue the numbering, so
  (`.chain`, `.iteration`) identifies one pilot draw. Documented.
- **`b52f803`.** The temporary names get a further "t" prefix until no variable
  has them. Sorted mapping (`5fab350`) exposed this case. Other orders were
  already broken, such as a polynomial in a and tmp1.

## 3. Validation

Run on an installed copy of `f060950`, in a fresh R process, with both pilot
backends wherever they apply.

| Check | Rejection pilot | HMC pilot |
|---|---|---|
| (x² − 1, 0.1x − 0.1): max \|x_proj − 1\| | 0 | 0 |
| (x, x − 1), empty | "Only 0 of 8 pilot draws projected" | same |
| x² − x: share at 1 | 0.495 | 0.375 (components disconnected: the documented chain caveat) |
| x(x − 1)² | stops: 2,927 of 4,410 draws stretched | stops: 2,204 of 4,410 |
| `"x"`, `vars = c("x", "y", "z")`, `w = 0.01` | 20 rows; \|y, z\| ≤ 0.0099 | 20 rows; ≤ 0.0095 |
| (x, x²), `vars = c("x", "z")`, `w = 1` | \|x_proj\| ≤ 6e-11; z over (−0.99, 1.00) | \|x_proj\| ≤ 7e-11; z over (−1.00, 0.99) |
| `mp("tmp2 + 2 a - 1")`, `w = 2` | — | runs template `2_1_vn_w` |
| Ellipse, Σ = [.01 .006; .006 .01]: KS of arclength; cov(X − U) / Σ | p = 0.78; 1.01, 1.01, 1.03 | thin = 25: p = 0.039 (see below); 1.03, 1.00, 1.03 |

**Chain labels.** The live window case did not need a second pilot run, because
the pilot samples inside the widened window. The mocked test covers the
numbering across runs.

**HMC ellipse p-value.** The 0.039 came from one seed.
- Seeds 9–14 gave p = 0.64, 0.095, 0.28, 0.22, 0.055, 0.011, more small values
  than independent draws would give.
- A test that tolerates autocorrelation finds no bias. Each chain's share of
  points with |x_proj| > 1.5 counts as one observation: over 64 chains (8 seeds
  × 8 chains, `thin = 5`) the mean is 0.3609 against the exact 0.3586, with
  SE 0.0040, so z = 0.57.
- So a KS test on thinned HMC draws still overstates their independence.
- These fixes do not touch the HMC pilot or the projection of regular points.
- Exact checks belong to the rejection pilot, whose draws are independent.

## 4. Open questions for David

1. **Repeated factors: correct instead of stop?** For one polynomial, the ratio
   distance / residual estimates each draw's multiplicity k. Keeping each draw
   with probability 1/k would restore uniform weights. For systems the
   correction is a product over normal directions, which one draw cannot
   estimate, so systems would still stop.
2. **A public latent sampler.** For example `rvariety()`, or `output =
   "latent"`. Latent points could then be reused across noise covariances, and
   other backends could be added later.
3. **Exact sampling of finite root sets.** For zero-dimensional varieties,
   `variety_solve()` could find the points, and U would be uniform over them.
   That would avoid the chain caveat that gives x² − x its 0.375.

## 5. Verification

- `devtools::test()` on `f060950`: **876 passed, 0 failed, 0 warnings, 1
  skipped** (the env-gated install test). That is up from 848.
- Each fix's new tests failed against the pre-fix copies of its functions. For
  fixes whose tests call internal functions directly, the HEAD versions were
  bound in both the namespace and the attached package, with the namespace as
  their environment, so that mocks applied.
- CRAN-style `R CMD check` on `f060950`: **Status: OK**. testthat summary:
  `[ FAIL 0 | WARN 2 | SKIP 3 | PASS 874 ]`, with the usual 2 warnings and 3
  skips.
