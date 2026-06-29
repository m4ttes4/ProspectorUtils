module DataUtils

export AbstractProspectResult, ProspectResults, ProspectorObs, ProspectorBestFit, ProspectorSampling
export Estimator, BestFit, Median, WeightedMedian, has_weights, estimate
export quantiles
export get_mass, maggies2μJy, get_z, labels, bestfit, n_bins
export get_obs_sflux, get_obs_swave, get_obs_serr, get_obs_pflux, get_obs_pwave, get_obs_perr
export get_bf_sflux, get_bf_swave, get_bf_pflux, get_bf_pwave, get_bf_sed, get_bf_cont, get_bf_calib, get_full_grid

using HDF5
using DataFrames
using JSON
using StatsBase
using Printf

const WEIGHTS_COL = "weights"

# =============================================================================
# Types
# =============================================================================

abstract type AbstractProspectResult end

"Container for a full Prospector output: MCMC chain plus the obs/bestfit/sampling groups."
mutable struct ProspectResults <: AbstractProspectResult
    chain::DataFrame
    runparams::Dict{String,Any}
    bestfit::Dict{String,Any}
    sampling::Dict{String,Any}
    obs::Dict{String,Any}
end

mutable struct ProspectorObs <: AbstractProspectResult
    obs::Dict{String,Any}
end
ProspectorObs(r::ProspectResults) = ProspectorObs(r.obs)

mutable struct ProspectorBestFit <: AbstractProspectResult
    bestfit::Dict{String,Any}
end
ProspectorBestFit(r::ProspectResults) = ProspectorBestFit(r.bestfit)

mutable struct ProspectorSampling <: AbstractProspectResult
    sampling::Dict{String,Any}
end
ProspectorSampling(r::ProspectResults) = ProspectorSampling(r.sampling)

# =============================================================================
# Point estimators
# =============================================================================

"""
How to collapse a posterior into a point estimate.
`BestFit` reads the maximum-likelihood parameters, `Median` the chain median,
`WeightedMedian` the chain median weighted by the (dynesty) importance weights.
"""
abstract type Estimator end
struct BestFit <: Estimator end
struct Median <: Estimator end
struct WeightedMedian <: Estimator end

"True when the chain carries a `weights` column (dynesty importance weights)."
has_weights(p::ProspectResults) = WEIGHTS_COL in names(p.chain)

estimate(p::ProspectResults, param::AbstractString, ::BestFit) = bestfit(p, param)
estimate(p::ProspectResults, param::AbstractString, ::Median) = median(skipmissing(p.chain[!, param]))
function estimate(p::ProspectResults, param::AbstractString, ::WeightedMedian)
    has_weights(p) || throw(ArgumentError("chain has no `$WEIGHTS_COL` column; WeightedMedian unavailable"))
    return median(p.chain[!, param], StatsBase.weights(p.chain[!, WEIGHTS_COL]))
end

# =============================================================================
# Weighted quantiles
# =============================================================================

# =============================================================================
# HDF5 reader
# =============================================================================

function _read_group(f::HDF5.File, path::String; at::Bool=false)
    if at
        haskey(f, path) || return Dict{String,Any}()
        res = Dict{String,Any}()
        for key in keys(attrs(f[path]))
            res[key] = _maybe_json(read_attribute(f[path], key))
        end
        return res
    else
        haskey(f, path) ? read(f[path]) : nothing
    end
end

function _maybe_json(x)
    (x isa AbstractString || x isa Vector{UInt8}) || return x
    str = x isa AbstractString ? x : String(x)
    return try
        JSON.parse(str)
    catch
        str
    end
end

_read_all_groups(f::HDF5.File; at::Bool=false) = (
    obs      = _read_group(f, "obs"; at=at),
    bestfit  = _read_group(f, "bestfit"; at=at),
    sampling = _read_group(f, "sampling"; at=at),
)

function _get_run_params(f::HDF5.File)::Dict{String,Any}
    try
        raw = read_attribute(f, "run_params")
        return JSON.parse(raw isa AbstractString ? raw : String(raw))
    catch e
        @warn "Failed to parse run_params: $e"
        return Dict{String,Any}()
    end
end

"""
    ProspectResults(filename; verbose=true)

Load a Prospector HDF5 output file. All datasets and attributes are read once,
the MCMC chain is assembled into a `DataFrame` (with a `weights` column when the
sampler provides importance weights).
"""
function ProspectResults(filename::String; verbose::Bool=true)
    h5open(filename, "r") do f
        attrs_nt = _read_all_groups(f; at=true)
        data_nt = _read_all_groups(f; at=false)
        run_params = _get_run_params(f)

        merge!(data_nt.obs, attrs_nt.obs)
        merge!(data_nt.bestfit, attrs_nt.bestfit)
        merge!(data_nt.sampling, attrs_nt.sampling)

        # convenience: expose the wavelength grids on the bestfit group too
        data_nt.bestfit["wavelength"] = data_nt.obs["wavelength"]
        data_nt.bestfit["phot_wave"] = data_nt.obs["phot_wave"]

        chain = _build_chain_df(data_nt.sampling, attrs_nt.sampling, run_params; verbose=verbose)
        return ProspectResults(chain, run_params, data_nt.bestfit, data_nt.sampling, data_nt.obs)
    end
end

"""
Build the MCMC chain `DataFrame` from pre-loaded sampling data. Layout depends on
the sampler (dynesty vs emcee); views avoid copies until DataFrame construction.
A `weights` column is appended when `sampling/weights` is present.
"""
function _build_chain_df(sampling_data::Dict, sampling_attrs::Dict, run_params::Dict; verbose::Bool=true)::DataFrame
    labels = get(sampling_attrs, "theta_labels", String[])
    chain = get(sampling_data, "chain", Float64[])
    isempty(chain) && return DataFrame()

    has_dyn = get(run_params, "dynesty", false)
    has_e = get(run_params, "emcee", false)
    has_dyn && has_e && @warn "Both dynesty and emcee flags set—using dynesty layout"

    if has_dyn
        verbose && @info "Building chain from Dynesty sampler"
        mat = PermutedDimsArray(chain, (2, 1))
    elseif has_e
        verbose && @info "Building chain from Emcee sampler"
        nl, nch, nwk = size(chain)
        mat = reshape(chain, nl, nch * nwk)'
    else
        mat = chain
    end

    df = if !isempty(labels) && size(mat, 2) == length(labels)
        DataFrame(mat, labels, makeunique=true)
    else
        !isempty(labels) && @warn "Column/label mismatch ($(size(mat, 2)) vs $(length(labels))); using generic labels"
        DataFrame(mat, :auto)
    end

    weights = get(sampling_data, WEIGHTS_COL, nothing)
    if weights !== nothing
        w = vec(weights)
        length(w) == nrow(df) ? (df[!, WEIGHTS_COL] = Float64.(w)) :
            (verbose && @warn "weights length $(length(w)) ≠ chain rows $(nrow(df)); skipping")
    end
    return df
end

# =============================================================================
# Accessors
# =============================================================================

get_z(p::ProspectResults) = get(p.runparams, "redshift", nothing)
"Parameter names in the chain, excluding the `weights` column."
labels(p::ProspectResults) = filter(!=(WEIGHTS_COL), names(p.chain))
bestfit(p::ProspectResults) = get(p.bestfit, "parameter", nothing)
n_bins(p::ProspectResults) = count(n -> occursin("logsfr_ratios", n), names(p.chain)) + 1

"Maximum-likelihood value of a single parameter, looked up via `theta_labels`."
function bestfit(p::ProspectResults, param::AbstractString)
    i = findfirst(==(param), p.sampling["theta_labels"])
    i === nothing && throw(KeyError(param))
    return p.bestfit["parameter"][i]
end

get_obs_sflux(p::ProspectResults) = get(p.obs, "spectrum", nothing)
get_obs_swave(p::ProspectResults) = get(p.obs, "wavelength", nothing)
get_obs_serr(p::ProspectResults) = get(p.obs, "unc", nothing)
get_obs_pflux(p::ProspectResults) = get(p.obs, "maggies", nothing)
get_obs_pwave(p::ProspectResults) = get(p.obs, "phot_wave", nothing)
get_obs_perr(p::ProspectResults) = get(p.obs, "maggies_unc", nothing)

get_bf_sflux(p::ProspectResults) = get(p.bestfit, "spectrum", nothing)
get_bf_swave(p::ProspectResults) = get(p.bestfit, "wavelength", nothing)
get_bf_pflux(p::ProspectResults) = get(p.bestfit, "photometry", nothing)
get_bf_pwave(p::ProspectResults) = get(p.bestfit, "phot_wave", nothing)
get_bf_cont(p::ProspectResults) = get(p.bestfit, "speccont", nothing)
get_bf_calib(p::ProspectResults) = get(p.bestfit, "speccal", nothing)
get_bf_sed(p::ProspectResults) = get(p.bestfit, "full_sed", nothing)
get_full_grid(p::ProspectResults) = get(p.bestfit, "full_grid", nothing)

maggies2μJy(mag::Real) = mag * 1e6 * 3631
maggies2μJy(::Nothing) = nothing

_scalar(x) = x isa AbstractArray ? first(x) : x

"""
    get_mass(p, est=BestFit()) -> Float64

log₁₀ of the *surviving* stellar mass: `logmass + log₁₀(mfrac)`, where `logmass`
is the total formed mass (estimated with `est`) and `mfrac` the surviving fraction.
"""
function get_mass(p::ProspectResults, est::Estimator=BestFit())
    return estimate(p, "logmass", est) + log10(_scalar(p.bestfit["mfrac"]))
end

"""
    quantiles(p, param; q=[0.16, 0.5, 0.84], weighted=has_weights(p))

Posterior quantiles of `param` from the chain. Weighted by importance weights when
`weighted` is true (the default whenever a `weights` column is present).
"""
function quantiles(p::ProspectResults, param::AbstractString;
                   q=[0.16, 0.5, 0.84], weighted::Bool=has_weights(p))
    col = p.chain[!, param]
    weighted || return quantile(collect(skipmissing(col)), q)
    has_weights(p) || throw(ArgumentError("chain has no `$WEIGHTS_COL` column"))
    return quantile(col, StatsBase.weights(p.chain[!, WEIGHTS_COL]), q)
end

# =============================================================================
# Display
# =============================================================================

function Base.show(io::IO, p::ProspectResults)
    println(io, "ProspectResults: ", nrow(p.chain), "×", length(labels(p)), " chain",
            has_weights(p) ? " (weighted)" : "", ", z=", get_z(p))
    params = labels(p)
    isempty(params) && return
    @printf(io, "%-22s %12s %12s %12s\n", "parameter", "16%", "50%", "84%")
    for name in params
        q16, q50, q84 = quantiles(p, name)
        @printf(io, "%-22s %12.4g %12.4g %12.4g\n", name, q16, q50, q84)
    end
    print(io, "estimator: BestFit")
end

Base.show(io::IO, p::ProspectorObs) = print(io, "ProspectorObs(", length(p.obs), " keys: ", join(keys(p.obs), ", "), ")")
Base.show(io::IO, p::ProspectorBestFit) = print(io, "ProspectorBestFit(", length(p.bestfit), " keys: ", join(keys(p.bestfit), ", "), ")")
Base.show(io::IO, p::ProspectorSampling) = print(io, "ProspectorSampling(", length(p.sampling), " keys: ", join(keys(p.sampling), ", "), ")")

end # module
