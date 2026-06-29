module ProspectorUtilsCairoMakieExt

# CairoMakie-dependent plotting behaviour for ProspectorUtils. Loaded automatically
# when both ProspectorUtils and CairoMakie are in scope.

using CairoMakie
using StatsBase
using Printf
using ProspectorUtils
using ProspectorUtils.PlotUtils
using ProspectorUtils.SFHUtils: zred_to_agebins, logmass_to_masses, get_sfh, get_agebins, sfh_lookback

import ProspectorUtils.PlotUtils: get_plot_function, get_plot_data, render!

const MK = CairoMakie.Makie

# =============================================================================
# Plot function per type
# =============================================================================

get_plot_function(::Union{ObsSpecPlot,BFSpecPlot,BFSedPlot,MaskedSpecPlot,BinnedSFH}) = CairoMakie.stairs!
get_plot_function(::Union{ObsPhotoPlot,BFPhotoPlot,MaskedPhotoPlot}) = CairoMakie.scatter!
get_plot_function(::DensityPlot) = CairoMakie.density!
get_plot_function(::LineIndicators) = CairoMakie.linesegments!
get_plot_function(::SedErrors) = CairoMakie.band!
get_plot_function(::PhotoErrors) = CairoMakie.errorbars!
get_plot_function(::SFHErrors) = CairoMakie.band!

# =============================================================================
# Calibration
# =============================================================================

function _apply_calibration!(res::ProspectResults, ::String, p, out::Vector{<:Real})
    calib = get_bf_calib(res)
    pop!(p.kwargs, :calibration)
    return out ./ calib
end
function _apply_calibration!(::ProspectResults, cal::Vector{<:Real}, p, out::Vector{<:Real})
    pop!(p.kwargs, :calibration)
    return out ./ cal
end

# =============================================================================
# get_plot_data
# =============================================================================

const _XY = Union{typeof(CairoMakie.stairs!),typeof(CairoMakie.lines!),typeof(CairoMakie.scatter!)}

get_plot_data(r::ProspectResults, ::ObsPhotoPlot, ::_XY) =
    (copy(get_obs_pwave(r)), maggies2μJy.(copy(get_obs_pflux(r))))
get_plot_data(r::ProspectResults, ::BFPhotoPlot, ::_XY) =
    (copy(get_bf_pwave(r)), maggies2μJy.(copy(get_bf_pflux(r))))
get_plot_data(r::ProspectResults, ::BFSedPlot, ::_XY) =
    (copy(get_full_grid(r)), maggies2μJy.(copy(get_bf_sed(r))))

function get_plot_data(r::ProspectResults, p::ObsSpecPlot, ::_XY)
    y = maggies2μJy.(copy(get_obs_sflux(r)))
    x = copy(get_obs_swave(r))
    haskey(p.kwargs, :calibration) && (y .= _apply_calibration!(r, p.kwargs[:calibration], p, y))
    return x, y
end

function get_plot_data(r::ProspectResults, p::BFSpecPlot, ::_XY)
    y = maggies2μJy.(copy(get_bf_sflux(r)))
    x = copy(get_bf_swave(r))
    haskey(p.kwargs, :calibration) && (y .= _apply_calibration!(r, p.kwargs[:calibration], p, y))
    return x, y
end

function get_plot_data(r::ProspectResults, ::SedErrors, ::typeof(CairoMakie.band!))
    return copy(get_obs_swave(r)), maggies2μJy.(copy(get_obs_sflux(r))), maggies2μJy.(copy(get_obs_serr(r)))
end

function get_plot_data(r::ProspectResults, ::PhotoErrors, ::typeof(CairoMakie.errorbars!))
    return copy(get_obs_pwave(r)), maggies2μJy.(copy(get_obs_pflux(r))), maggies2μJy.(copy(get_obs_perr(r)))
end

function get_plot_data(r::ProspectResults, p::MaskedSpecPlot, ::_XY)
    y = maggies2μJy.(copy(get_obs_sflux(r)))
    x = copy(get_obs_swave(r))
    z = Bool.(get(r.obs, "mask", nothing))
    haskey(p.kwargs, :calibration) && (y .= _apply_calibration!(r, p.kwargs[:calibration], p, y))
    return x, y, z
end

function get_plot_data(r::ProspectResults, ::MaskedPhotoPlot, ::typeof(CairoMakie.scatter!))
    return copy(get_obs_pwave(r)), maggies2μJy.(copy(get_obs_pflux(r))), Bool.(get(r.obs, "phot_mask", nothing))
end

# Stairs rendering needs n+1 x-edges and a matching y; pad y by repeating the
# youngest bin so the staircase draws flat from lookback≈0. Plot-only concern —
# get_sfh returns the n physical SFR values.
function get_plot_data(r::ProspectResults, ::BinnedSFH, ::_XY)
    _, sfh = get_sfh(r)
    return sfh_lookback(get_agebins(r)), vcat(first(sfh), sfh)
end

# =============================================================================
# Masked-segment helper
# =============================================================================

function _segments(x, y, mask)
    segs = Tuple{Vector,Vector,Bool}[]
    start = 1
    current = mask[1]
    for i in 2:length(x)
        if mask[i] != current
            push!(segs, (x[start:i], y[start:i], current))
            start = i
            current = mask[i]
        end
    end
    push!(segs, (x[start:end], y[start:end], current))
    return segs
end

# =============================================================================
# render!
# =============================================================================

function render!(ax::CairoMakie.Axis, result::ProspectResults, plot::Union{AbstractProspectPlot,AbstractErrors})
    f = get_plot_function(plot)
    f(ax, get_plot_data(result, plot, f)...; plot.kwargs...)
    return ax
end

function render!(ax::CairoMakie.Axis, result::ProspectResults, cmd::PlotCommand)
    items = Any[cmd.multipliers..., cmd.base_plot]
    sort!(items, by=get_zorder)
    foreach(item -> render!(ax, result, item), items)
    return ax
end

function render!(ax::CairoMakie.Axis, result::ProspectResults, plot::MaskedPhotoPlot)
    f = get_plot_function(plot)
    x, y, mask = get_plot_data(result, plot, f)
    for i in eachindex(x)
        mask[i] && f(ax, x[i], y[i]; plot.kwargs...)
    end
    return ax
end

function render!(ax::CairoMakie.Axis, result::ProspectResults, plot::MaskedSpecPlot)
    f = get_plot_function(plot)
    wave, flux, mask = get_plot_data(result, plot, f)
    if isnothing(wave) || isnothing(flux) || isnothing(mask)
        @warn "render!(MaskedSpecPlot): wave/flux/mask is nothing"
        return nothing
    end
    for (xs, ys, is_true) in _segments(wave, flux, mask)
        is_true ? f(ax, xs, ys; plot.kwargs...) : f(ax, xs, ys; plot.mask...)
    end
    return ax
end

function render!(layout, result::ProspectResults, plot::DensityPlot)
    f = get_plot_function(plot)
    params = filter(x -> x ∉ plot.ignored, names(result.chain))
    isempty(params) && (@warn "DensityPlot: all parameters excluded"; return layout)

    n_rows = ceil(Int, sqrt(length(params)))
    n_cols = ceil(Int, length(params) / n_rows)
    idx = 1
    for i in 1:n_rows, j in 1:n_cols
        idx > length(params) && break
        data = result.chain[!, params[idx]]
        idx += 1
        isempty(data) && (@warn "DensityPlot: no data for $(params[idx-1])"; continue)

        use_sci = maximum(data) > 1e3
        xspan = -(extrema(data)...)
        q16, q50, q84 = quantile(data, [0.16, 0.5, 0.84])
        ax = CairoMakie.Axis(layout[i, j];
            title=_confidence_title(q50 - q16, q50, q84 - q50, "$(params[idx-1]) | "),
            xgridvisible=false, ygridvisible=false,
            xtickformat=use_sci ? "{:.1e}" : "{:.1f}",
            xticklabelrotation=use_sci ? π / 4 : 0.0,
            xscale=xspan > 1e6 ? log10 : identity)
        CairoMakie.hidespines!(ax, :t, :r, :l)
        CairoMakie.hideydecorations!(ax)
        f(ax, data; plot.kwargs...)
    end
    return layout
end

function render!(ax::CairoMakie.Axis, result::ProspectResults, plot::BinnedSFH)
    f = get_plot_function(plot)
    f(ax, get_plot_data(result, plot, f)...; plot.kwargs...)
    return ax
end

function render!(ax::CairoMakie.Axis, result::ProspectResults, plot::SFHErrors; step=:center)
    z = get_z(result)
    ratio_labels = filter(n -> occursin("logsfr_ratios", n), names(result.chain))
    nbins = length(ratio_labels) + 1
    agebins = zred_to_agebins(z, nbins)

    err_low = Float64[]
    err_upp = Float64[]
    for name in ratio_labels
        low, med, upp = quantile(result.chain[!, name], [0.16, 0.5, 0.84])
        push!(err_low, med - low)
        push!(err_upp, upp - med)
    end
    pushfirst!(err_low, err_low[begin]); pushfirst!(err_upp, err_upp[begin])
    push!(err_low, err_low[end]); push!(err_upp, err_upp[end])
    reverse!(err_low); reverse!(err_upp)

    logmass = bestfit(result, "logmass")
    ratios = [median(result.chain[!, name]) for name in ratio_labels]
    mass = logmass_to_masses(logmass, ratios, agebins)

    initial_bins = first.(agebins)
    final_bins = last.(agebins)
    sfr = mass ./ (10.0 .^ final_bins .- 10.0 .^ initial_bins)
    xx = 10.0^initial_bins[1] / 1e9 .+ 10.0 .^ final_bins ./ 1e9
    pushfirst!(xx, 1e-9)
    yy = log10.(vcat(sfr[1], sfr))

    s = CairoMakie.stairs!(ax, xx, yy; color=:transparent, step=step)
    _band_error(s, err_upp, err_low, ax; plot.kwargs...)
    return ax
end

# =============================================================================
# Formatting helpers
# =============================================================================

function _confidence_title(low, mid, high, label::String)
    largest = max(abs(high), abs(low))
    if largest == 0
        digits = mid == 0 ? 0 : max(0, 1 - round(Int, log10(abs(mid))))
        return @sprintf("%.*f", digits, mid)
    end

    digits = max(0, 1 - round(Int, log10(largest)))
    if digits > 4  # scientific notation
        scale = 10.0^(digits - 1)
        sup = MK.superscript(@sprintf("-%d", digits - 1))
        if round(low, digits=digits) == round(high, digits=digits)
            return MK.rich(label, "(", @sprintf("%.1f", mid * scale), " ± ",
                @sprintf("%.1f", high * scale), ") × 10", sup)
        end
        return MK.rich(label, "( ", @sprintf("%.1f", mid * scale),
            MK.left_subsup(@sprintf("-%.1f", low * scale), @sprintf("+%.1f", high * scale)),
            " ) × 10", sup)
    end

    if round(low, digits=digits) == round(high, digits=digits)
        return MK.rich(label, @sprintf("%.*f", digits, mid), " ± ", @sprintf("%.*f", digits, high))
    end
    return MK.rich(label, @sprintf("%.*f", digits, mid),
        MK.left_subsup(@sprintf("-%.*f", digits, low), @sprintf("+%.*f", digits, high)))
end

_stairpts(s) = let pts = CairoMakie.to_value(s.plots[1].converted[1])
    ([p[1] for p in pts], [p[2] for p in pts])
end

function _band_error(s, errupp, errlow, ax; kwargs...)
    bins, flux_bins = _stairpts(s)
    k = 1
    err_upper, err_lower = errupp[1], errlow[1]
    for i in 1:length(bins)-1
        if !iseven(i)
            err_upper, err_lower = errupp[k], errlow[k]
            k += 1
        end
        lw = fill(flux_bins[i] - err_lower, 2)
        up = fill(flux_bins[i] + err_upper, 2)
        CairoMakie.band!(ax, [bins[i], bins[i+1]], lw, up; kwargs...)
    end
end

end # module
