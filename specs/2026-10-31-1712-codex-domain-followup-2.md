# vnorm — Codex's Second Review of Domain Noise: Findings Log

**Purpose:** Record of the follow-up to Codex's review of `d99e67e..94f1aba`,
logged in `2026-10-31-1337-codex-domain-review.md`. Covers each finding's
verdict, the fixes with their commits, the numerical contract the latent
sampler now keeps, and the validation.

**Repository:** `dkahle/vnorm`, branch `main` (local; not pushed).
**Toolchain:** R 4.6.1 (arm64 macOS), CmdStan 2.38.0 (cmdstanr 0.9.0), mpoly 1.1.2.900.

---

## 1. Verdicts

All five findings reproduced as reported, and all are fixed.

| # | Codex finding | Reproduction | Commit |
|---|---|---|---|
| 1 (P1) | Translation reopens the empty-variety bug | (x − 1e9, x − 1e9 − 1): eight latent points at 1e9 + 0.5. (x, x − 1) with a free z near 1e9: the same | `67dfd14` |
| 2 (P1) | Free coordinates change least-squares weights | (x, 2x − 2): mean of x 0.80, but 0.50 with `vars = c("x", "z")`. Domain (x² − x, x³ − x²): share at 1 0.502, but 0.428 with z | `f51648e` |
| 3 (P2) | The rounding-floor stop biases a component on the window's edge | x² − x with w = 1: share at 1 0.41 (0.504 with w = 1.1) | `f0bcca5`, `9b06a0d` |
| 4 (P2) | The multiplicity cap removes components | x⁶ − x⁵: every latent point at 1, after a generic warning | `94586bd` |
| 5 (P2) | Numerical uncertainty becomes latent noise | (x − 1)³ with sd 1e-5: latent RMS error 4.6e-6, noise RMS 1.10 times the sd, no warning | `f4d429c` |

## 2. The latent sampler's numerical contract

Codex asked for a clearer contract before more heuristic correction.
`project_rows_closest()` and `rvnorm_domain_latent()` now keep these rules.

- **Free coordinates** never enter pseudoinverses, ranks, the classification of
  a system, or tolerances.
  - `pdvnorm_pseudoinverse_cols()` works on the columns of the variables the
    equations involve.
  - It serves the rejection kernel, the Gauss–Newton steps, and the tangent
    steps (`f51648e`).
- **Tolerances are distances:** `tol × scale` plus the rounding error of the
  polynomial's own coordinates, 8 eps |u|. The latent sampler passes
  `scale = pilot_sd`. Translation no longer loosens them (`67dfd14`).
- **Acceptance:** a projection is accepted only when it converged and every
  equation vanishes at the point, to within 1e-6 × scale plus each equation's
  evaluation error.
- **Precision:** each accepted projection reports the distance within which
  the equations' evaluation errors leave it.
  - Simple root: the error over |∇g|, or through J⁺ for a system.
  - Factor of multiplicity k: d (error / |g(x₀)|)^(1/k).
  - Every latent point must be within 1% of the noise's smallest sd
    (`f4d429c`).
- **Multiplicity:** 1 up to the degree, corrected by keeping each draw with
  probability 1/k. A system with a repeated equation stops (`94586bd`).
- **How each pilot draw ends:**
  - did not converge: a failure, which stops sampling above 5%;
  - converged off the variety: at an overdetermined system's least-squares
    point;
  - landed far from its draw: a branch jump;
  - imprecise: stops sampling above 1%.

  Failures, off-variety and far draws are dropped, with one warning past 1%.
  The stop conditions run after every pilot run.
- **Polishing:** projections onto simple roots take up to two Newton steps past
  the rounding floor and land within rounding error of the variety
  (`9b06a0d`).
- **Window:** points within their precision of the window's edge count as
  inside (`f0bcca5`).

## 3. Found along the way

- **Off-variety least-squares points.** For (x² − x, x³ − x²), the
  least-squares pilot puts about 3% of its draws at x ≈ 0.552. That is a
  least-squares stationary point off the variety, where g = (−0.247, −0.137).
  Their projections jump to 0 or 1 and are dropped, with the warning. The
  latent shares stay 0.50.
- **The precision estimate at the rounding floor.** It first assumed 4 eps.
  But Gauss–Newton stops at 64 eps, up to 1 + 2.2e-14 for x² − x from 1.02.
  So with w = 1, `f0bcca5` alone left the root at 1 with 0.480 of 4,000
  points, 2.5 SE short. Polishing simple roots (`9b06a0d`) fixed that: the
  share is 0.496.
- **Open: an HMC pilot on a non-reduced polynomial at tiny noise can hang.**
  `rvnorm(400, mp("(x - 1)^3"), sd = 1e-5, noise = "domain", w = c(0.999,
  1.001))` with the default HMC pilot ran over 12 minutes, until killed. A
  rerun capped at 2 minutes timed out.
  - At `pilot_sd` 1e-6, the pilot's log density on the expanded cubic is
    rounding noise.
  - The precision check can only act once the pilot returns.
  - The rejection pilot stops in seconds with the precision error.
  - A check of the evaluation floor before sampling could fail fast. It is not
    done.

## 4. Validation

Run on an installed copy of `9b06a0d`, in a fresh R process.

| Check | Result |
|---|---|
| (x − 1e9, x − 1e9 − 1), rejection and HMC; (x, x − 1) with z near 1e9 | all stop: "Only 0 of 8 pilot draws projected" |
| x − 1e9, HMC, 200 draws | latent points exactly at 1e9 |
| (x, 2x − 2) with a free z: mean of x | rejection 0.798; HMC 0.800 (least squares: 0.8) |
| (x² − x, x³ − x²), 4,000 draws, share at 1 | 0.499 for `vars` x; 0.497 for x and z |
| x² − x, 4,000 draws, share at 1 | w = 1: 0.496; w = 1.1: 0.499 |
| x⁶ − x⁵, 4,000 draws, share at 1 | 0.513 at seed 53. Over seeds 61–66, mean 0.4964 (SE 0.0032) |
| (x − 1)³, sd 1e-5, rejection | stops: "At 450 of 450 pilot draws, the projection locates the variety only to …" |
| x − 1, sd 1e-5 | latent points exactly 1; noise RMS 0.988 times the sd |
| Ellipse, Σ = [.01 .006; .006 .01] | KS of arclength p = 0.78 (rejection); cov(X − U) / Σ 1.01, 1.01, 1.03; HMC max \|g(U)\| 2.2e-16 |
| Circle in 3-D (system) | z_proj = 0; radius exact to 1.1e-16 |
| (x² + y² − 1)³ | max radius error 2.4e-5; KS of the angle p = 0.52 |
| Earlier Codex examples | overdetermined: every point at 1; (x, x − 1): stops; `tmp2` names: run |

## 5. Verification

- `devtools::test()` on `9b06a0d`: **905 passed, 0 failed, 0 warnings, 1
  skipped** (the env-gated install test). That is up from 886.
- Each fix's new tests failed against the pre-fix copies of its functions.
- Changed expectations in existing tests:
  - The top-up test now uses a realistic share of failing draws, since
    failures past 5% stop sampling.
  - The all-failing mock now stops at once with "did not converge", rather
    than after five runs with "Only 0 of 10".
- CRAN-style `R CMD check`:
  - on `f0bcca5`, before the polishing: **Status: OK**,
    `[ FAIL 0 | WARN 2 | SKIP 3 | PASS 902 ]`;
  - on `9b06a0d`: **Status: OK**, `[ FAIL 0 | WARN 2 | SKIP 3 | PASS 903 ]`.
  - The 2 warnings and 3 skips are the usual ones.
