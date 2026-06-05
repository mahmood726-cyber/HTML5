# DTA Meta-Analysis Suite

A single-file, offline browser tool for **diagnostic test accuracy (DTA) meta-analysis**.
Enter 2×2 counts (TP, FP, FN, TN) per study and the suite computes pooled sensitivity
and specificity, likelihood ratios, the diagnostic odds ratio (DOR), a Moses–Littenberg
SROC curve, a Spearman threshold-effect test, and a Deeks' funnel-plot asymmetry test,
with a Fagan nomogram and QUADAS-2 inputs in the UI.

## Offline use

Open `index.html` directly in a browser — no network access is required.

- Plotly is vendored locally as `plotly.min.js` (loaded with `<script src="plotly.min.js">`).
- Google Fonts links were removed; the page falls back to system fonts.
- `grep -E 'https?://' index.html` returns nothing — there are zero external references.

## Layout

| File | Purpose |
|------|---------|
| `index.html` | The application (UI + DOM + Plotly wiring). |
| `engine.js`  | Pure statistical functions, extracted verbatim from the original inline script. The single source of truth — the HTML loads it before its own inline script. |
| `tests.js`   | Pure Node tests for `engine.js`. Every expected value is hand-derived. |
| `plotly.min.js` | Vendored Plotly 2.35.2 (offline). |

## Tests

```
node tests.js
```

Prints per-check `PASS`/`FAIL` lines, ends with `N passed, M failed`, and exits non-zero
on any failure. Coverage includes a hand-worked 2×2 study (Se, Sp, logit(Se), logit(FPR),
DOR), the DOR-sign identity, the SROC sign-flip identity, the χ²₂ bivariate-region constant,
logit clamping at 0/1, `normalCDF(0)=0.5`, DerSimonian–Laird random effects, Wilson
intervals, the Spearman threshold test, and the Deeks' ESS / regression guards.

## Engine functions

`logit`, `invLogit`, `wilsonCI`, `calcStudyStats`, `univariateRE`, `bivariateModel`,
`calcSROC`, `thresholdEffectTest`, `normalCDF`, `deeksFunnelTest`.

## Fixes applied during revival

- **Offline**: vendored Plotly 2.35.2 locally and removed the Google Fonts `<link>`; the
  page now has zero external references.
- **Single source of truth**: extracted the pure statistics into `engine.js` and deleted the
  inline duplicates, so the HTML and the tests run the exact same code.
- **Tests added**: `tests.js` with 49 hand-derived checks (none previously existed).
- **Renamed** the extension-less `DTA` file to `index.html` for GitHub Pages.

### Statistics review (no methodology changes)

The DTA-specific correctness rules were checked independently and the existing code was
found **already correct**; nothing was changed:

- **DOR sign** — `DOR = exp(μ₁ + μ₂)` (here `PLR/NLR`, and the per-study cross-product
  ratio), the additive form, not the buggy subtraction `exp(μ₁ − μ₂)`.
- **SROC sign** — the FPR→specificity conversion uses `logit(1 − FPR) = −logit(FPR)`, the
  correct sign flip.
- **logit clamping** — `logit(p)` clamps `p` to `[1e-4, 1 − 1e-4]`, so `logit(0)`/`logit(1)`
  are finite (no `−Inf` path).
- **normalCDF(0) = 0.5** — the Abramowitz–Stegun approximation returns `0.5` at `x = 0`.

The SROC AUC is a trapezoidal integral of the Moses–Littenberg curve (labelled "AUC (SROC)"),
not an HSROC-Φ AUC, so the Φ-vs-logistic rule does not apply. The SROC region ellipses are a
simplified display aid drawn in percentage units, not a χ²₂ statistical prediction region; the
correct √χ²₂ = 2.4477 constant is documented in `tests.js` for any future region computation.
