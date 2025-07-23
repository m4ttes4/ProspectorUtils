module SFHUtils
# -----------------------------------------------------------------------------
# External dependencies
# -----------------------------------------------------------------------------
using PythonCall
using Statistics
using DataFrames
using ..DataUtils

# -----------------------------------------------------------------------------
# Public API
# -----------------------------------------------------------------------------
export zred_to_agebins, build_agebins, logmass_to_masses, get_sfh, pyshow

"""
    SFHUtils – Star‑Formation History utilities

A lightweight collection of helpers mirroring the Prospector conventions for
age–bin construction and mass‑partitioning.

Public API
----------
* [`zred_to_agebins`] – build age bins from redshift.
* [`build_agebins`]    – construct Prospector‑style age bins from universe age.
* [`logmass_to_masses`] – convert a total mass + SFR ratios to per‑bin masses.

All routines avoid external dependencies except for `PythonCall` → `astropy` to
translate redshift into universe age.  If `astropy` is unavailable an
informative error is thrown.
"""


# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

pyshow() = let 
    cosmo = pyimport("astropy")
    display(cosmo)
end

"Return the age of the universe at redshift `zred` (Gyr) using Planck18 cosmology."
function _universe_age(zred::Real)
    # _cosmology === nothing &&
    #     throw(Error("Astropy cosmology is required for zred_to_agebins (missing astropy)",))
    let
        cosmo = pyimport("astropy.cosmology")
        cosmo === nothing && throw(Error("Astropy cosmology is required for zred_to_agebins (missing astropy)",))
        return pyconvert(Float64, cosmo.Planck18.age(zred).value)
    end
end

@inline _mydiff(bin::NTuple{2,T}) where {T<:Real} = 10.0^bin[2] - 10.0^bin[1]
@inline _mydiff(bin::Vector{T}) where {T<:Real} = 10.0^bin[2] - 10.0^bin[1]


"""
    zred_to_agebins(zred, nbins) -> Vector{NTuple{2,Float64}}

Build Prospector‑style age bins (log₁₀ years) for a universe at redshift `zred`.

Arguments
---------
* `zred::Real`  – Cosmological redshift.
* `nbins::Int`  – Number of bins **≥ 4**.

Notes
-----
* Cosmology: Planck18 via *astropy*.
* Oldest bin upper limit: 90 % of the universe age.
"""
function zred_to_agebins(zred::Real, nbins::Integer)
    tuniv = _universe_age(zred)               # [Gyr]
    return build_agebins(tuniv=tuniv, nbins=nbins)
end

"""
    build_agebins(; tuniv=13.7, nbins=7) -> Vector{NTuple{2,Float64}}

Return Prospector‑style age bins for a universe age `tuniv` (Gyr).

Keyword arguments
-----------------
* `tuniv::Real = 13.7` – Universe age in **Gyr**.
* `nbins::Int  = 7`    – Number of bins **≥ 4**.

Each element of the returned vector is `(log10(tₗ), log10(tᵤ))` in **years**.
"""
function build_agebins(; tuniv::Real=13.7, nbins::Integer=7)

    nbins < 4 && throw(ArgumentError("nbins must be ≥ 4 (got $nbins)"))

    tbinmax = (tuniv * 0.9) * 1e9
    lim1, lim2 = 7.4772, 8.0

    # Compute logarithmically spaced bins
    log_bins = range(lim2, stop=log10(tbinmax), length=nbins - 2) |> collect
    agelims = vcat([0, lim1], log_bins, [log10(tuniv * 1e9)])

    # Construct the agebins matrix
    agebins = [(agelims[i], agelims[i+1]) for i in 1:(length(agelims)-1)]#
    #agebins = hcat(agelims[1:end-1], agelims[2:end])
end

"""
    logmass_to_masses(logmass, logsfr_ratios, agebins) -> Vector{Float64}

Partition the total stellar mass across age bins following Prospector.

Parameters
----------
* `logmass::Real`               – log₁₀(∑ Mᵢ).
* `logsfr_ratios::AbstractVector{<:Real}` – log₁₀(SFRⱼ / SFRⱼ₊₁) for `nbins-1` bins.
* `agebins::AbstractVector{<:NTuple{2,<:Real}}` – Age bin limits from `build_agebins`.

Returns
-------
`Vector{Float64}` with the individual Mᵢ (same order as `agebins`).
"""
function logmass_to_masses(
    logmass::Real,
    logsfr_ratios::AbstractVector{<:Real},
    agebins::AbstractVector{<:NTuple{2,<:Real}}
)::Vector{Float64}

    nbins = length(agebins)
    length(logsfr_ratios) == nbins - 1 ||
        throw(DimensionMismatch("Expected $(nbins - 1) SFR ratios, got $(length(logsfr_ratios))"))

    # Convert to Float64 for type stability and clamp values
    sratios = map(x -> 10.0^clamp(float(x), -10.0, 10.0), logsfr_ratios)
    dt = map(_mydiff, agebins)

    if dt[1] == 0.0
        throw(ArgumentError("First age bin width is zero, cannot normalize mass coefficients."))
    end

    coeffs = Vector{Float64}(undef, nbins)
    coeffs[1] = 1.0

    acc = 1.0
    for j in 2:nbins
        acc *= sratios[j - 1]
        coeffs[j] = dt[j] / (dt[1] * acc)
    end

    m1 = 10.0^float(logmass) / sum(coeffs)
    return m1 .* coeffs
end


function logmass_to_masses(p::ProspectResults)::Vector{Float64}
    labels = filter(x -> occursin("logsfr_ratios", x), names(p.chain))
    ratios = map(name -> median(skipmissing(p.chain[!, name])), labels)
    agebins = zred_to_agebins(get_z(p), length(ratios))

    logmass = float(median(skipmissing(p.chain[!, "logmass"])))
    return logmass_to_masses(logmass, ratios, agebins)
end


function logmass_to_masses(
    p::ProspectResults,
    zred::Real,
    quantiles::AbstractVector{<:Real}
)::Vector{Vector{Float64}}

    logmasses = quantile(p.chain.logmass, quantiles)
    sfr_labels = filter(name -> occursin("logsfr_ratios", name), names(p.chain))
    ratios = map(name -> median(skipmissing(p.chain[!, name])), sfr_labels)
    agebins = zred_to_agebins(float(zred), length(ratios))

    return [logmass_to_masses(mass, ratios, agebins) for mass in logmasses]
end

function logmass_to_masses(p::ProspectResults, agebins::Vector{Tuple{Float64, Float64}})
    logmass = median(p.chain.logmass)
    labels = names(p.chain)[findall(x -> occursin.("logsfr_ratios", x), names(p.chain))]
    ratios = [median(p.chain[!, name]) for name in labels]
    return logmass_to_masses(logmass, ratios, agebins)
end


"""
    get_sfh(results::ProspectResults; normalize=false, safe=true, return_mass=false)
        -> (lookback::Vector{Float64}, sfh::Vector{Float64})

Compute the star formation history (SFH) from a `ProspectResults` object.

# Keyword arguments
- `normalize::Bool=false`: if true, normalize SFR to unity.
- `safe::Bool=true`: apply checks to avoid NaNs, negative SFRs, or zero-duration bins.
- `return_mass::Bool=false`: if true, return tuple `(lookback, sfh, mass)`.

# Returns
- `lookback`: vector of lookback times in Gyr (bin edges).
- `sfh`: log₁₀ of SFR values in each bin.
- optionally, `mass`: per-bin stellar mass.

# Notes
- `zred_to_agebins` is used to recompute bins at `get_z(results)`.
- Assumes age bins are given in log₁₀(years).
- Bins with zero width or invalid entries are skipped if `safe=true`.
"""
function get_sfh(results::ProspectResults; normalize::Bool=false, safe::Bool=true, return_mass::Bool=false)
    z = get_z(results)

    nbins = n_bins(results)#sum(occursin.("logsfr_ratios", names(results.chain)))+1 
    agebins = zred_to_agebins(z, nbins)
    mass = logmass_to_masses(results, agebins)


    initial_bins = first.(agebins)
    final_bins = last.(agebins)

    Δt = @. 10.0^final_bins - 10.0^initial_bins
    sfr = mass ./ Δt

    if safe
        for i in eachindex(sfr)
            @inbounds x = sfr[i]
            @inbounds sfr[i] = (isfinite(x) && x > 0) ? x : 1e-30
        end
    end
    
    lookback = 10 .^ initial_bins[1] ./ 1e9 .+ 10 .^ final_bins ./ 1e9
    lookback = vcat(1e-9, lookback)  # Prepend tiny lookback time to match bins

    if normalize
        sfr ./= sum(sfr)
    end

    sfh = log10.(vcat(sfr[1], sfr))

    return return_mass ? (lookback, sfh, mass) : (lookback, sfh)
end


# function get_sfh_quantiles(results::ProspectResults; normalize::Bool=false, save::Bool=true, return_mass::Bool=false)
#     z = get_z(results)

#     nbins = n_bins(results)#sum(occursin.("logsfr_ratios", names(results.chain)))+1 
#     agebins = zred_to_agebins(z, nbins)
#     mass = logmass_to_masses(results, z, [0.16, 0.5, 0.84])

#     initial_bins = first.(agebins)
#     final_bins = last.(agebins)

#     Δt = @. 10.0^final_bins - 10.0^initial_bins
#     sfr = map(x -> x ./ Δt, mass)

#     if safe
#         for bin in sfr
#             for i in eachindex(bin)
#                 @inbounds x = bin[i]
#                 @inbounds bin[i] = (isfinite(x) && x > 0) ? x : 1e-30
#             end
#         end
#     end

#     if normalize
#         @warn "Currently normalize the SFH quantile is broken"
#         sfr ./= sum(sfr)
#     end

#     sfh = log10.(vcat([sfr[1]], sfr))
#     return sfh
# end

end # module
