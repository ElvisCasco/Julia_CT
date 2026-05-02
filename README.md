# Julia_CT — Julia translation of Cameron &amp; Trivedi's *Microeconometrics Using Stata*

Julia translations of the Stata code from *Microeconometrics Using Stata* (Cameron &amp; Trivedi) as Quarto notebooks.

## Files

- **`Cameron & Trivedi_stata.qmd`** — original Stata code (1st edition) as a Quarto notebook using `stata_setup`.
- **`Cameron & Trivedi_stata_vJulia.qmd`** — **Julia translation** of the Stata notebook above. Covers all 18 chapters + appendices.
- **`Cameron & Trivedi_stata_2ed_V2.qmd`** — 2nd-edition Stata notebook (work in progress).
- **`Cameron & Trivedi_stata_2ed_V2_testJulia.qmd`** — Julia test translation for the 2nd edition.
- **Data folders** (`musr/`, `mus2/`, `mma10252005/`, `nhh2017cameron/`) — sample datasets from the book's companion materials.
- **`auto.dta`**, **`census.dta`**, etc. — Stata built-in / sample data.

## Running

### Stata notebook
Requires Stata 19 (or compatible) and `stata_setup` via Python:
```bash
pip install stata-setup
```
Render with:
```bash
quarto render "Cameron & Trivedi_stata.qmd"
```

### Julia notebook
Requires Julia 1.10+ with these packages:
```julia
using Pkg
Pkg.add([
    "CairoMakie", "CategoricalArrays", "CSV", "Chain", "DataFrames",
    "Downloads", "FixedEffectModels", "GLM", "HypothesisTests",
    "PrettyTables", "ReadStatTables", "RegressionTables", "Statistics",
    "StatsBase", "CategoricalArrays", "LinearAlgebra", "KernelDensity",
    "Loess", "RDatasets", "Random", "Printf", "Optim", "Distributions",
    "SpecialFunctions", "MixedModels", "Combinatorics"
])
```
Render with:
```bash
quarto render "Cameron & Trivedi_stata_vJulia.qmd"
```

## Scope

The Julia translation covers all 18 chapters of *Microeconometrics Using Stata* plus Appendices A and B.

### Chapters

1. **Stata basics** — interactive use, command syntax, do-files, scalars/matrices, macros, looping, user-written commands.
2. **Data management and graphics** — data types, input/output, manipulation, reshaping, graphical display.
3. **Linear regression basics** — OLS in levels and logs, specification analysis, prediction, sampling weights, native Julia matrix-algebra OLS.
4. **Simulation** — pseudorandom-number generators, sample-mean distribution, computing integrals, simulation for regression.
5. **GLS regression** — GLS/FGLS, heteroskedastic models, systems of linear regressions (SUR), survey data (weights/clusters/strata).
6. **Linear instrumental-variables regression** — 2SLS, weak-instrument diagnostics, LIML/JIVE-style robust inference, 3SLS.
7. **Quantile regression** — QR for medical-expenditure data, generated heteroskedastic data, and count data (Machado–Santos Silva).
8. **Linear panel-data models: Basics** — pooled, within (FE), between, RE, first-difference, long panels, panel-data management.
9. **Linear panel-data models: Extensions** — panel IV, Hausman–Taylor (analytic/bootstrap/jackknife SEs), Arellano–Bond (difference GMM), mixed linear models, clustered data.
10. **Nonlinear regression methods** — nonlinear least squares and related estimators.
11. **Nonlinear optimization methods** — `Optim.jl` replacements for Stata's `ml` / Mata routines (Newton, BHHH, BFGS).
12. **Testing methods** — Wald, LR, LM/score, conditional moment, and non-nested tests.
13. **Bootstrap methods** — pairs/residual bootstrap, percentile-t, jackknife, bootstrap inference.
14. **Binary outcome models** — logit, probit, marginal effects, goodness-of-fit.
15. **Multinomial models** — multinomial / conditional / nested logit, multinomial probit via GHK simulator, mixed logit via simulated MLE (Halton draws).
16. **Tobit and selection models** — censored Tobit, Heckman two-step, sample-selection models.
17. **Count-data models** — Poisson, NegBin, hurdle, zero-inflated, finite mixtures.
18. **Nonlinear panel models** — FE logit, RE Tobit, FE/RE Poisson, panel NegBin.

### Appendices

- **A. Programming in Julia** (replaces "Programming in Stata") — functions, control flow, modules.
- **B. Native Julia matrix algebra** (replaces Mata) — linear algebra patterns used throughout the book.

## Notes

- PDFs of the textbooks are **not included** (copyright); please obtain them from the publisher.
- Rendered HTML outputs are not included; regenerate with `quarto render`.
- Bootstrap standard errors may differ slightly from Stata due to differing RNGs (Julia's MersenneTwister vs Stata's internal RNG). Point estimates match exactly.
