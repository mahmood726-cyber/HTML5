# E156 Protocol — DTA Meta-Analysis Suite

- **Project:** DTA Meta-Analysis Suite (`HTML5`)
- **Revived:** 2026-06-05
- **Type:** Offline single-file browser tool (diagnostic test accuracy meta-analysis)
- **Dashboard:** https://mahmood726-cyber.github.io/HTML5/

## What changed

- Vendored Plotly 2.35.2 locally (`plotly.min.js`) and removed the Google Fonts link; zero external references remain.
- Extracted the pure statistical functions into `engine.js` (verbatim) and deleted the inline duplicates so the app and tests share one source of truth.
- Added `tests.js` with 49 hand-derived Node checks.
- Reviewed the DTA-specific statistics (DOR sign, SROC sign, logit clamping, normalCDF(0)) and confirmed the existing code is already correct — no methodology was changed.
- Renamed the extension-less `DTA` file to `index.html` and added `.nojekyll`, `.gitignore`, and `README.md`.

## Body (E156 draft — CURRENT BODY)

How accurately can a browser-only tool reproduce standard diagnostic test accuracy meta-analysis from raw 2×2 study counts without any server or network dependency? The suite ingests per-study true-positive, false-positive, false-negative, and true-negative counts for an arbitrary set of diagnostic studies. It applies a 0.5 continuity correction, logit-transforms sensitivity and specificity, pools them by DerSimonian–Laird random effects, and fits a Moses–Littenberg SROC curve. Across a worked twelve-study troponin example it returns pooled sensitivity, specificity, positive and negative likelihood ratios, a diagnostic odds ratio computed as exp(μ₁+μ₂), a trapezoidal SROC area, and Spearman threshold and Deeks' asymmetry diagnostics. Forty-nine hand-derived Node tests pin the DOR sign, the logit(1−FPR)=−logit(FPR) SROC flip, logit clamping at zero and one, and normalCDF(0)=0.5. The findings indicate the extracted engine is numerically faithful to the original inline implementation, justifying reuse of the same audited code in both the page and the test harness. Hierarchical REML and exact HSROC likelihoods remain outside this single-file scope. SUBMITTED: [ ]
