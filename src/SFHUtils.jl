module SFHUtils

using PythonCall
using StatsBase
using DataFrames
using ..DataUtils

export zred_to_agebins, build_agebins, get_agebins, logmass_to_masses
export get_sfh, get_sfr, sfh_lookback

"""
SFHUtils – Star-Formation History utilities mirroring Prospector conventions for
age-bin construction and mass partitioning.

Age bins are reconstructed from redshift using the Planck18 cosmology via
*astropy* (`zred_to_agebins`). Stored `bestfit/agebins` are intentionally not
used because they may be inconsistent with the result redshift.
"""
SFHUtils

# -----------------------------------------------------------------------------
# Age bins
# -----------------------------------------------------------------------------

"Age of the universe at redshift `zred` (Gyr) from Planck18 (astropy)."
function _universe_age(zred::Real)
    cosmo = pyimport("astropy.cosmology")
    return pyconvert(Float64, cosmo.Planck18.age(zred).value)
end

@inline _bin_width(bin) = 10.0^bin[2] - 10.0^bin[1]

"""
    build_agebins(; tuniv=13.7, nbins=7) -> Vector{NTuple{2,Float64}}

Prospector-style age bins for a universe age `tuniv` (Gyr). Each element is
`(log₁₀ t_low, log₁₀ t_high)` in years. Matches Prospector's bin construction.
"""
function build_agebins(; tuniv::Real=13.7, nbins::Integer=7)
    nbins < 4 && throw(ArgumentError("nbins must be ≥ 4 (got $nbins)"))
    tbinmax = (tuniv * 0.9) * 1e9
    lim1, lim2 = 7.4772, 8.0
    log_bins = collect(range(lim2, stop=log10(tbinmax), length=nbins - 2))
    agelims = vcat([0.0, lim1], log_bins, [log10(tuniv * 1e9)])
    return [(agelims[i], agelims[i+1]) for i in 1:(length(agelims)-1)]
end

"Build Prospector-style age bins for a universe at redshift `zred` (needs astropy)."
zred_to_agebins(zred::Real, nbins::Integer) = build_agebins(tuniv=_universe_age(zred), nbins=nbins)

"""
    get_agebins(p) -> Vector{NTuple{2,Float64}}

Age bins for a result, reconstructed from redshift (`zred_to_agebins`, requires
astropy). Stored `bestfit/agebins` are ignored because some Prospector outputs
contain bins that are not limited by the universe age at `get_z(p)`.
"""
get_agebins(p::ProspectResults) = zred_to_agebins(get_z(p), n_bins(p))

# -----------------------------------------------------------------------------
# Mass partitioning
# -----------------------------------------------------------------------------

"""
    logmass_to_masses(logmass, logsfr_ratios, agebins) -> Vector{Float64}

Partition the total formed stellar mass across age bins following Prospector.
`logmass` is log₁₀(∑ Mᵢ); `logsfr_ratios` holds the `nbins-1` log₁₀(SFRⱼ/SFRⱼ₊₁).
"""
function logmass_to_masses(logmass::Real, logsfr_ratios::AbstractVector{<:Real},
                           agebins::AbstractVector{<:NTuple{2,<:Real}})::Vector{Float64}
    nbins = length(agebins)
    length(logsfr_ratios) == nbins - 1 ||
        throw(DimensionMismatch("expected $(nbins - 1) SFR ratios, got $(length(logsfr_ratios))"))

    sratios = @. 10.0^clamp(logsfr_ratios, -10.0, 10.0)
    dt = map(_bin_width, agebins)
    dt[1] == 0.0 && throw(ArgumentError("first age bin has zero width; cannot normalize"))

    coeffs = Vector{Float64}(undef, nbins)
    coeffs[1] = 1.0
    acc = 1.0
    for j in 2:nbins
        acc *= sratios[j-1]
        coeffs[j] = dt[j] / (dt[1] * acc)
    end

    m1 = 10.0^logmass / sum(coeffs)
    return m1 .* coeffs
end

"""
    logmass_to_masses(p, est=BestFit(); agebins=get_agebins(p))

Per-bin formed masses for a result, using the point estimate `est` for `logmass`
and the SFR ratios.
"""
function logmass_to_masses(p::ProspectResults, est::Estimator=BestFit();
                           agebins=get_agebins(p))::Vector{Float64}
    ratio_labels = filter(n -> occursin("logsfr_ratios", n), labels(p))
    logmass = estimate(p, "logmass", est)
    ratios = [estimate(p, n, est) for n in ratio_labels]
    return logmass_to_masses(logmass, ratios, agebins)
end

# -----------------------------------------------------------------------------
# Star-formation history
# -----------------------------------------------------------------------------

"""
    sfh_lookback(agebins) -> Vector{Float64}

Lookback-time bin edges (Gyr) for a stairs-style SFH: a leading near-zero point
followed by the upper edge of each age bin (offset by the first lower edge),
matching Prospector's plotting convention.
"""
function sfh_lookback(agebins::AbstractVector{<:NTuple{2,<:Real}})
    offset = 10.0^first(agebins[1]) / 1e9
    return vcat(1e-9, offset .+ 10.0 .^ last.(agebins) ./ 1e9)
end

"""
    get_sfh(p, est=BestFit(); normalize=false, safe=true, return_mass=false)
        -> (lookback, sfh[, mass])

Star-formation history from a result, one value per age bin. `sfh` is
log₁₀(SFR) = log₁₀(mass / Δt); `lookback` is each bin's older-edge lookback time
(Gyr). `safe` floors non-finite/≤0 SFR to 1e-30; `normalize` scales SFR to unit
sum; `return_mass` also returns the per-bin formed mass. All outputs have length
`n_bins(p)` — the leading-point duplication a stairs plot needs is a rendering
concern handled in the plotting layer, not here.
"""
function get_sfh(p::ProspectResults, est::Estimator=BestFit();
                 normalize::Bool=false, safe::Bool=true, return_mass::Bool=false)
    agebins = get_agebins(p)
    mass = logmass_to_masses(p, est; agebins=agebins)

    Δt = @. 10.0^last(agebins) - 10.0^first(agebins)
    sfr = mass ./ Δt

    safe && map!(x -> (isfinite(x) && x > 0) ? x : 1e-30, sfr, sfr)
    normalize && (sfr ./= sum(sfr))

    lookback = @. 10.0^last(agebins) / 1e9
    sfh = log10.(sfr)
    return return_mass ? (lookback, sfh, mass) : (lookback, sfh)
end

"Log₁₀ of the mean SFR over the most recent `nbins` SFH bins."
function get_sfr(p::ProspectResults, est::Estimator=BestFit(); nbins::Integer=3)
    _, sfh = get_sfh(p, est)
    return log10(mean(10.0 .^ @view sfh[1:nbins]))
end

end # module
