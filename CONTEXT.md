# ProspectorUtils

Julia post-processing for [Prospector](https://github.com/bd-j/prospector) SED-fitting
outputs: load an HDF5 run, extract posterior estimates, and reconstruct star-formation
histories. Plotting lives in a CairoMakie package extension.

## Language

**ProspectResults**:
The full output of one Prospector run — the MCMC `chain` plus the `obs`, `bestfit`
and `sampling` groups.
_Avoid_: result dict, run output.

**chain**:
Posterior MCMC samples as a `DataFrame`; columns are parameters, with an optional
`weights` column carrying dynesty importance weights.
_Avoid_: samples, trace.

**bestfit**:
Maximum-likelihood parameters and models. Holds `parameter`, `mfrac`, `agebins`.

**obs**:
Observational inputs: spectrum, photometry (maggies), uncertainties, masks.

**sampling**:
Sampler metadata: `theta_labels`, `chain`, `weights`.

**logmass**:
log₁₀ of the total *formed* stellar mass (∑ Mᵢ over age bins).

**mfrac**:
Surviving stellar-mass fraction (linear, ~0.6). Surviving mass = `logmass + log10(mfrac)`.

**agebins**:
Age-bin limits as `(log₁₀ t_low, log₁₀ t_high)` in years.

**SFH** (star-formation history):
Per-bin formed mass divided by bin duration → SFR per bin.

**Estimator**:
How a posterior is collapsed to a point estimate: `BestFit`, `Median`, `WeightedMedian`.
_Avoid_: point estimate method.

## Relationships

- A **ProspectResults** owns one **chain** and the **bestfit**/**obs**/**sampling** groups.
- An **Estimator** maps a **chain** parameter to a scalar; `WeightedMedian` requires the
  `weights` column.
- **agebins** + **logmass** + SFR ratios → per-bin masses → **SFH**.
- **agebins** come from the **bestfit** group when stored, else are recomputed from redshift.

## Flagged ambiguities

- `get_mass` previously computed `logmass - mfrac` (mixing a log quantity with a linear
  fraction). Resolved: surviving mass is `logmass + log10(mfrac)`.
- "logmass for the SFH" was ambiguous between best-fit and chain median. Resolved: the
  point estimate is explicit via an `Estimator`; default is `WeightedMedian` when weights
  exist, else `Median`.

## Example dialogue

> **Dev:** "When I call `get_sfh`, whose `logmass` does it use?"
> **Domain expert:** "Whatever estimator you pass. With dynesty weights it defaults to the
> weighted median of the chain; you can force `BestFit()` for the max-likelihood value."
