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

## 4. David's decisions

1. **Repeated factors: correct single polynomials (`497c70a`).**
   - Each pilot draw's multiplicity k = round(distance / residual) is kept
     with probability 1/k. That draws from R's RNG only when some k > 1, so a
     square-free polynomial's draws are unchanged for a given seed.
   - With 8,000 rejection draws, the share at 1 of {0, 1} is 0.5065 for
     x² − x, 0.5039 for x(x − 1)², and 0.5020 for x(x − 1)³ (SE 0.0056).
   - A squared circle is uniform in angle.
   - Systems still stop, since their correction is a product over normal
     directions.
2. **Public latent sampler: not now.** `"tibble"` output already returns the
   latent points.
3. **Exact sampling of finite root sets: later.**

**Projecting repeated factors.** The full test suite on `c964731` flagged a
warning: 37 of 2,450 projections had failed for the squared circle. Two
numerical problems on repeated factors, both now fixed:
- **Refinement.** On (x² + y² − 1)² the gradient vanishes on the circle and
  comes from cancelling terms near it. The tangent stalled at 3.5e-10, above
  its 2e-10 tolerance. On (x² + y² − 1)³ the computed gradient at the landing
  point is rounding noise (about 6e-18).
  - **Fix.** The refinement now runs only where distance ≈ residual, i.e. on
    simple components. On repeated ones Gauss–Newton's point stands. It
    steps along the reduced factor's normal, so it lands within O(d²) of the
    closest point.
- **Gauss–Newton itself.** On the cube, once f³ falls below the expanded
  polynomial's rounding error, each step is noise over a gradient that
  vanishes like f². That threw points across the circle, up to 1.39 away.
  - **Fix.** Gauss–Newton now stops once |g| is within the rounding error of
    its terms, the floor `variety_feasible()` already uses.
- **Jump limit.** On a k-fold component the pilot's tube is k times wider,
  so the 10 `pilot_sd` limit now scales with k. A ratio above 4.5 counts as
  a jump to another branch.

Projections of 2,000 rejection-pilot draws (`pilot_sd` 0.01) had no failures.
Errors against `project_onto_variety()` on the reduced polynomial:

| Polynomial | Error |
|---|---|
| ellipse | 3e-10 |
| squared circle | 1.2e-7 |
| cubed circle | 2.4e-5 |
| squared ellipse | 1.5e-3, Gauss–Newton's O(d²) drift |

With 8,000 draws, the shares at 1 of {0, 1} were 0.5065, 0.5071 and 0.5070
for multiplicities 1–3, with no warnings.

**Found along the way (`c964731`).** For a single polynomial,
`rejection_sampler()` expanded |∇g|² into one polynomial and took its log.
- Near a repeated root that expansion cancels below 0: −3.6e-15 for the true
  9.0e-16, for x(x − 1)³ at x = 1 − 1e-4.
- The log was then NaN, with warnings, for codomain and domain rejection
  sampling alike.
- It now squares the evaluated partials.

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
- After David's decisions, on the final tree `47a5d7a`:
  - `devtools::test()`: **886 passed, 0 failed, 0 warnings, 1 skipped**.
  - `R CMD check`: **Status: OK**; testthat summary
    `[ FAIL 0 | WARN 2 | SKIP 3 | PASS 884 ]`.
  - The new tests of `497c70a`, `c964731` and `47a5d7a` failed against their
    pre-fix code.
