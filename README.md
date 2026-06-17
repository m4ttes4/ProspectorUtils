# ProspectorUtils

Julia post-processing for [Prospector](https://github.com/bd-j/prospector) SED-fitting
outputs. Load an HDF5 run, pull posterior estimates, and reconstruct star-formation
histories. Optional plotting via a CairoMakie extension.

See [`CONTEXT.md`](CONTEXT.md) for the domain glossary and [`docs/adr/`](docs/adr) for
design decisions.

## Install

```julia
using Pkg
Pkg.develop(path="path/to/ProspectorUtils")
```

The core depends only on HDF5, DataFrames, JSON, StatsBase (weighted quantiles), PythonCall,
plus the `Printf` stdlib. CairoMakie is a weak dependency, loaded only when you `using CairoMakie`.

Age bins are read directly from the HDF5 output (`bestfit/agebins`) when present. Only when
they are absent does the package fall back to recomputing them from redshift via astropy's
Planck18 cosmology (through PythonCall); that path therefore requires `astropy` in the
active Python environment.

## Usage

```julia
using ProspectorUtils

p = ProspectResults("run.h5")     # load an HDF5 output
show(p)                            # summary table: 16% / 50% / 84% per parameter

get_z(p)                           # redshift
labels(p)                          # parameter names (excludes the weights column)

# Point estimates — choose how the posterior is collapsed
get_mass(p)                        # log10 surviving stellar mass (default estimator)
get_mass(p, BestFit())             # max-likelihood
estimate(p, "logmass", Median())   # any parameter, any estimator

# Posterior quantiles (weighted automatically when dynesty weights are present)
quantiles(p, "logmass")                       # [16%, 50%, 84%]
quantiles(p, "dust2"; q=[0.05, 0.5, 0.95])
weighted_median(p.chain.logmass, p.chain.weights)

# Star-formation history
lookback, sfh = get_sfh(p)                 # uses default estimator
lookback, sfh, mass = get_sfh(p; return_mass=true)
get_sfr(p; nbins=3)                        # mean log10(SFR) over recent bins
```

### Estimators

| Estimator        | Source                                            |
|------------------|---------------------------------------------------|
| `BestFit()`      | max-likelihood parameters (`bestfit`)             |
| `Median()`       | chain median                                      |
| `WeightedMedian()` | chain median weighted by dynesty weights        |

`default_estimator(p)` returns `WeightedMedian()` when the chain has a `weights` column,
otherwise `Median()`.

## Plotting

```julia
using CairoMakie          # activates the plotting extension
using ProspectorUtils

fig = Figure(); ax = Axis(fig[1, 1])
render!(ax, p, ObsSpecPlot())
render!(ax, p, BFSpecPlot() |> Calibration())
```

## Tests

```julia
using Pkg; Pkg.test()
```
