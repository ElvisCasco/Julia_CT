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

The Julia translation covers all 18 chapters plus Appendices A and B. Highlights of translated techniques:

- OLS, GLS, IV (2SLS, LIML, JIVE, GMM, 3SLS)
- Quantile regression + quantile count regression (Machado–Santos Silva)
- Fixed-effects / random-effects panel (xtreg equivalents)
- Hausman-Taylor (with analytic, bootstrap, and jackknife SEs)
- Arellano-Bond (difference GMM)
- Nonlinear regression, optimization (Optim.jl replaces `ml` / Mata)
- Binary, multinomial, and count models (logit, probit, Poisson, NegBin)
- Multinomial probit via GHK simulator
- Mixed logit via simulated MLE (Halton draws)
- Tobit and sample-selection models
- Finite mixture models (Poisson / NegBin)
- Nonlinear panel models (FE logit, RE Tobit, FE/RE Poisson, panel NegBin)

## Notes

- PDFs of the textbooks are **not included** (copyright); please obtain them from the publisher.
- Rendered HTML outputs are not included; regenerate with `quarto render`.
- Bootstrap standard errors may differ slightly from Stata due to differing RNGs (Julia's MersenneTwister vs Stata's internal RNG). Point estimates match exactly.
