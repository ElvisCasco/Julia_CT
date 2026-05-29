# Julia_CT — Julia translation of Cameron &amp; Trivedi's *Microeconometrics Using Stata*

Julia translations of the Stata code from *Microeconometrics Using Stata* (Cameron &amp; Trivedi) as Quarto notebooks.

## Files

- **`Cameron & Trivedi_stata.qmd`** — original Stata code (1st edition) as a Quarto notebook using `stata_setup`.
- **`Cameron & Trivedi_stata_vJulia.qmd`** — **Julia translation** of the Stata notebook above. Covers all 18 chapters + appendices.
- **`ch01_Stata_basics.qmd` … `ch18_Nonlinear_panel-data_models.qmd`** — **per-chapter standalone notebooks** generated from the master Julia qmd. Each is self-contained and renders independently (see [Per-chapter notebooks](#per-chapter-notebooks) below).
- **`split_chapters.jl`** — script that produces the per-chapter qmds from the master Julia qmd by extracting only the chunk-1 helpers each chapter transitively uses.
- **`annotate_with_stata.jl`** — script that prepends the corresponding Stata source from `Cameron & Trivedi_stata.qmd` to each Julia chunk in a chapter qmd, as a `#= … =#` reference comment.
- **`Cameron & Trivedi_stata_2ed_V2.qmd`** — 2nd-edition Stata notebook (work in progress).
- **`Cameron & Trivedi_stata_2ed_V2_testJulia.qmd`** — Julia test translation for the 2nd edition.
- **Data folders** (`musr/`, `mus2/`, `mma10252005/`, `nhh2017cameron/`) — sample datasets from the book's companion materials.
- **`auto.dta`**, **`census.dta`**, etc. — Stata built-in / sample data.

## Per-chapter notebooks

In addition to the single master notebook (`Cameron & Trivedi_stata_vJulia.qmd`, ~38 k lines), the repo ships **eighteen standalone chapter files** — one per textbook chapter:

```
ch01_Stata_basics.qmd
ch02_Data_management_and_graphics.qmd
ch03_Linear_regression_basics.qmd
…
ch18_Nonlinear_panel-data_models.qmd
```

Each chapter file is **self-contained** and renders on its own. Its structure is:

1. **YAML preamble** (same as the master notebook).
2. **Imports + compat-shims** chunk — every `using` statement and the PrettyTables / Optim v3 compatibility shims.
3. **Helpers chunk** — only those chunk-1 helpers that this chapter actually references (transitively, via `split_chapters.jl`'s dependency closure).
4. **Chapter content** — the original example chunks from the master notebook. Each Julia chunk is prepended with the corresponding **Stata source** from `Cameron & Trivedi_stata.qmd` as a `#= … =#` reference block, e.g.

   ````markdown
   ```{julia}
   #=
   sysuse auto.dta
   summarize
   =#
   auto = DataFrame(ReadStatTables.readstat("auto.dta"))
   summarize_stata(auto)
   ```
   ````

This makes each chapter file useful as both a runnable Quarto notebook and a side-by-side Stata-↔-Julia reference.

### Regenerating the chapter files

The per-chapter qmds are derived; the master file is the source of truth. To regenerate after editing `Cameron & Trivedi_stata_vJulia.qmd`:

```bash
julia --project=@. split_chapters.jl                   # rebuild ch01…ch18
for f in ch??_*.qmd; do                                # re-attach Stata refs
  julia --project=@. annotate_with_stata.jl "$f"
done
```

(Each `annotate_with_stata.jl` run also writes a `*.qmd.bak` the first time it touches a file, so the un-annotated version stays recoverable.)

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
Render the master notebook with:
```bash
quarto render "Cameron & Trivedi_stata_vJulia.qmd"
```

Or render a single chapter:
```bash
quarto render ch04_Simulation.qmd
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

## Testing status

The master notebook has been executed end-to-end and verified for all chapters.

All packages and helper functions used anywhere in the master notebook are loaded and defined in two upfront sections, so any chapter can be run independently after these are evaluated:

- **Julia setup** — top-of-notebook cell that calls `using` on every required package (CairoMakie, CategoricalArrays, CSV, Chain, DataFrames, Downloads, FixedEffectModels, GLM, HypothesisTests, PrettyTables, ReadStatTables, RegressionTables, Statistics, StatsBase, LinearAlgebra, KernelDensity, Loess, RDatasets, Random, Printf, Optim, Distributions, SpecialFunctions, MixedModels, Combinatorics).
- **Chapter 1 — Stata basics** — defines all reusable helper functions used by later chapters (table formatting, summary statistics, regression-output helpers, panel utilities, bootstrap/jackknife wrappers, etc.).

If you want to run only a specific chapter from the master notebook, evaluate the *Julia setup* cell and the helpers in *Chapter 1* first. Alternatively, use the [per-chapter standalone notebooks](#per-chapter-notebooks) — each already includes the setup, only the helpers it needs, and the chapter's content in one self-contained file. Per-chapter render status is spot-checked (e.g. `ch04_Simulation.qmd` and `ch17_Count-data_models.qmd` have been rendered end-to-end); a full sweep of all 18 is in progress.

## Notes

- PDFs of the textbooks are **not included** (copyright); please obtain them from the publisher.
- Rendered HTML outputs are not included; regenerate with `quarto render`.
- Bootstrap standard errors may differ slightly from Stata due to differing RNGs (Julia's MersenneTwister vs Stata's internal RNG). Point estimates match exactly.
