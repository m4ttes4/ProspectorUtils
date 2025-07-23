module PlotUtils
using ..DataUtils
using ..SFHUtils


using CairoMakie
using Statistics
using Printf


# import Base.|>
export AbstractProspectPlot, AbstractSedPlot, AbstractPhotoPlot, AbstractErrors
export ObsSpecPlot, BFSpecPlot, DensityPlot, ObsPhotoPlot, BFPhotoPlot, BinnedSFH, BFSedPlot
export Errors, MaskedSpecPlot, Calibration, SFHErrors
export PlotCommand, Mask#, PhotoMask, SpecMask
export get_plot_function, get_error_plot_function
export get_plot_data
export render!



abstract type AbstractProspectPlot end
abstract type AbstractSedPlot <: AbstractProspectPlot end
abstract type AbstractPhotoPlot <: AbstractProspectPlot end

abstract type AbstractMultiplier end
abstract type AbstractErrors <: AbstractMultiplier end
# =====================================
# === SED Plot Structs ===
# ====================================="ObsSpecPlot(; kwargs...)\n\nSED plot based on observed data."

mutable struct MaskedSpecPlot <: AbstractSedPlot
    kwargs::Dict{Symbol,Any}
    mask::Dict{Symbol,Any}
end

# --- Plot mutable Structs Senza plotfunc ---
"ObsSpecPlot(; kwargs...)\n\nSED plot based on observed data."
mutable struct ObsSpecPlot <: AbstractSedPlot
    kwargs::Dict{Symbol,Any}
end
ObsSpecPlot(; kwargs...) = ObsSpecPlot(Dict{Symbol,Any}(kwargs))

"BFSpecPlot(; kwargs...)\nSED plot based on best-fit model output."
mutable struct BFSpecPlot <: AbstractSedPlot
    kwargs::Dict{Symbol,Any}
end
BFSpecPlot(; kwargs...) = BFSpecPlot(Dict{Symbol,Any}(kwargs))


"BFSedPlot(; kwargs...)\nSED plot based on best-fit model output."
mutable struct BFSedPlot <: AbstractSedPlot
    kwargs::Dict{Symbol,Any}
end
BFSedPlot(; kwargs...) = BFSedPlot(Dict{Symbol,Any}(kwargs))



"DensityPlot(ignored; kwargs...)\n\nDensity plot of best-fit parameter distributions."
mutable struct DensityPlot <: AbstractProspectPlot
    ignored::Vector{String}
    kwargs::Dict{Symbol,Any}
end
DensityPlot(ignored::Vector{String}; kwargs...) = DensityPlot(ignored, Dict{Symbol,Any}(kwargs))
DensityPlot() = DensityPlot(String[], Dict{Symbol,Any}())

"ObsPhotoPlot(; kwargs...)\n\nPhotometry plot for observed data."
mutable struct ObsPhotoPlot <: AbstractPhotoPlot
    kwargs::Dict{Symbol,Any}
end
ObsPhotoPlot(; kwargs...) = ObsPhotoPlot(Dict{Symbol,Any}(kwargs))

"BFPhotoPlot(; kwargs...)\n\nPhotometry plot for best-fit model predictions."
mutable struct BFPhotoPlot <: AbstractPhotoPlot
    kwargs::Dict{Symbol,Any}
end
BFPhotoPlot(; kwargs...) = BFPhotoPlot(Dict{Symbol,Any}(kwargs))

mutable struct BinnedSFH <: AbstractProspectPlot
    kwargs::Dict{Symbol, Any}
end
BinnedSFH(;kwargs...) = BinnedSFH(Dict(kwargs))


mutable struct SedErrors <: AbstractErrors
    kwargs::Dict{Symbol,Any}
end
SedErrors(;kwargs...) = SedErrors(Dict(kwargs))

mutable struct PhotoErrors <: AbstractErrors
    kwargs::Dict{Symbol,Any}
end
PhotoErrors(; kwargs...) = PhotoErrors(Dict(kwargs))

mutable struct SFHErrors <: AbstractErrors
    kwargs::Dict{Symbol,Any}
end
SFHErrors(;kwargs...) = SFHErrors(Dict(kwargs))
mutable struct Errors <: AbstractMultiplier
    kwargs::Dict{Symbol,Any}
end

Errors(; kwargs...) = Errors(Dict(kwargs))

# Implementazione dispatch
_errors(::AbstractSedPlot, kwargs) = SedErrors(kwargs)
_errors(::BFPhotoPlot) = error("Erros can be plotted only for OBS currently")

_errors(::AbstractPhotoPlot, kwargs) = PhotoErrors(kwargs)
_errors(::BFSpecPlot) = error("Erros can be plotted only for OBS currently")

_errors(::BinnedSFH, kwargs) = SFHErrors(kwargs) 
# Call method
function (self::Errors)(p::AbstractProspectPlot)
    _errors(p, self.kwargs)
end


struct Calibration
    # key::Union{Nothing,String}
    # poly::Union{Nothing,Vector{T}}
    calibration::Union{String, Vector{<: Real}}
end

Calibration() = Calibration("speccal")

function _apply_calibration!(res::ProspectResults, ::String, p::AbstractProspectPlot, out::Vector{<:Real})
    calib = get_bf_calib(res)
    pop!(p.kwargs, :calibration)
    return out ./ calib
end
function _apply_calibration!(::ProspectResults, cal::Vector{<:Real}, p::AbstractProspectPlot, out::Vector{<:Real})
    pop!(p.kwargs, :calibration)
    return out ./ cal
end


struct Mask <: AbstractMultiplier
    kwargs::Dict{Symbol, Any}
end
Mask(; kwargs...) = Mask(Dict(kwargs))
# struct PhotoMask <: AbstractMultiplier
#     kwargs::Dict{Symbol,Any}
# end
# struct SpecMask <: AbstractMultiplier
#     kwargs::Dict{Symbol,Any}
# end

# Mask(;kwargs...) = Mask(Dict(kwargs))

# _mask(::AbstractSedPlot, kwargs) = SpecMask(kwargs)
# _mask(::AbstractPhotoPlot, kwargs) = PhotoMask(kwargs)

# function (self::Mask)(p::AbstractProspectPlot)
#     _mask(p, self.kwargs)
# end

# =====================================
# ===  PlotCommand ===
# =====================================

struct PlotCommand
    base_plot::AbstractProspectPlot
    multipliers::Vector{<:AbstractMultiplier}
end


function Base.:|>(left::AbstractProspectPlot, right::Calibration)
    kwargs = copy(left.kwargs)
    kwargs[:calibration] = right.calibration
    return typeof(left)(kwargs)
end

function Base.:|>(left::MaskedSpecPlot, right::Calibration)
    kwargs = copy(left.kwargs)
    kwargs[:calibration] = right.calibration
    return MaskedSpecPlot(kwargs, left.mask)
end



function Base.:|>(left::AbstractSedPlot, right::Mask)
    return MaskedSpecPlot(left.kwargs, right.kwargs)
end

# 1. Plot base diventa un PlotCommand con il primo modificatore
function Base.:|>(left::AbstractProspectPlot, right::AbstractMultiplier)
    # Errors(left) = Sed/Photo Error
    multiplier_instance = right(left)
    return PlotCommand(left, [multiplier_instance])
end

function Base.:|>(left::PlotCommand, right::AbstractMultiplier)
    # Errors(left) = Sed/Photo Error
    multiplier_instance = right(left.base_plot)

    return PlotCommand(left.base_plot, [left.multipliers..., multiplier_instance])
end



# 2. Funzioni di Dispatch per la Funzione di Plotting
# ----------------------------------------------------

"Restituisce la funzione di plotting da usare per un dato AbstractProspectPlot."
function get_plot_function end # Definiamo la funzione generica

# Metodi di default (questi sono i metodi che l'utente può estendere)
get_plot_function(::ObsSpecPlot) = CairoMakie.stairs!
get_plot_function(::BFSpecPlot) = CairoMakie.stairs!
get_plot_function(::ObsPhotoPlot) = CairoMakie.scatter!
get_plot_function(::BFPhotoPlot) = CairoMakie.scatter!

get_plot_function(::BFSedPlot) = CairoMakie.stairs!

get_plot_function(::DensityPlot) = CairoMakie.density!

get_plot_function(::PhotoErrors) = CairoMakie.errorbars!
get_plot_function(::SedErrors) = CairoMakie.band!
# get_plot_function(::SpecMask) = CairoMakie.stairs!
# get_plot_function(::PhotoMask) = CairoMakie.scatter!
get_plot_function(::MaskedSpecPlot) = CairoMakie.stairs!


get_plot_function(::BinnedSFH) = CairoMakie.stairs!
get_plot_function(::SFHErrors) = CairoMakie.Band!

# =====================================
# === Plot data vaidations ===
# =====================================
function get_plot_data end

function get_segments(x, y, mask)
    segments = []
    start = 1
    current_mask = mask[1]
    for i in 2:length(x)
        if mask[i] != current_mask
            push!(segments, (x[start:i], y[start:i], current_mask))
            start = i
            current_mask = mask[i]
        end
    end
    # Aggiungi l'ultimo segmento
    push!(segments, (x[start:end], y[start:end], current_mask))
    return segments
end


# Caso generale: Dati per funzioni di plot (x, y) come scatter!, lines!, stairs!
function get_plot_data(result::ProspectResults, ::ObsPhotoPlot, ::Union{typeof(CairoMakie.scatter!),typeof(CairoMakie.lines!)})

    y = get_obs_pflux(result) |> copy .|> maggies2μJy
    x = get_obs_pwave(result) |> copy
    return x, y
end

function get_plot_data(result::ProspectResults, p::ObsSpecPlot, ::Union{typeof(CairoMakie.stairs!),typeof(CairoMakie.lines!)})

    y = get_obs_sflux(result) |> copy .|> maggies2μJy
    x = get_obs_swave(result) |> copy
    if haskey(p.kwargs, :calibration)
        y .= _apply_calibration!(result, p.kwargs[:calibration], p, y)
    end
    return x, y
end

# Caso generale: Dati per funzioni di plot (x, y) come scatter!, lines!, stairs!
function get_plot_data(result::ProspectResults, ::BFPhotoPlot, ::Union{typeof(CairoMakie.scatter!),typeof(CairoMakie.lines!)})

    y = get_bf_pflux(result) |> copy .|> maggies2μJy
    x = get_bf_pwave(result) |> copy 
    return x, y
end

function get_plot_data(result::ProspectResults, p::BFSpecPlot, ::Union{typeof(CairoMakie.stairs!),typeof(CairoMakie.lines!)})

    y = get_bf_sflux(result) |> copy .|> maggies2μJy
    x = get_bf_swave(result) |> copy 
    if haskey(p.kwargs, :calibration)
        y .= _apply_calibration!(result, p.kwargs[:calibration], p, y)
    end
    return x, y
end



function get_plot_data(result::ProspectResults, ::SedErrors, ::typeof(CairoMakie.band!))

    y = get_obs_sflux(result) |> copy .|> maggies2μJy
    x = get_obs_swave(result) |> copy 
    z = get_obs_serr(result) |> copy .|> maggies2μJy

    return x, y, z
end


function get_plot_data(result::ProspectResults, ::PhotoErrors, ::typeof(CairoMakie.errorbars!))

    y = get_obs_pflux(result) |> copy .|> maggies2μJy
    x = get_obs_pwave(result) |> copy 
    z = get_obs_perr(result) |> copy .|> maggies2μJy

    return x, y, z
end

# Caso generale: Dati per funzioni di plot (x, y) come scatter!, lines!, stairs!
function get_plot_data(result::ProspectResults, ::BFSedPlot, ::Union{typeof(CairoMakie.stairs!),typeof(CairoMakie.lines!)})

    y = get_bf_sed(result) |> copy .|> maggies2μJy
    x = get_full_grid(result) |> copy 
    return x, y
end

# function get_plot_data(result::ProspectResults, ::PhotoMask,::Union{typeof(CairoMakie.stairs!),typeof(CairoMakie.lines!)})

#     y = get_obs_pflux(result) |> copy .|> maggies2μJy
#     x = get_obs_pwave(result) |> copy .|> maggies2μJy
#     z = get(result.obs,"photo_mask", nothing)

#     return x, y, z
# end

# function get_plot_data(result::ProspectResults, ::SpceMask,::Union{typeof(CairoMakie.stairs!),typeof(CairoMakie.lines!)})

#     y = get_obs_pflux(result) |> copy .|> maggies2μJy
#     x = get_obs_pwave(result) |> copy .|> maggies2μJy
#     z = get(result.bestfit, "mask", nothing) 

#     return x, y, z
# end

function get_plot_data(result::ProspectResults, p::MaskedSpecPlot, ::Union{typeof(CairoMakie.stairs!),typeof(CairoMakie.lines!)})
    y = get_obs_sflux(result) |> copy .|> maggies2μJy
    x = get_obs_swave(result) |> copy 
    z = get(result.obs, "mask", nothing) .|> Bool

    if haskey(p.kwargs, :calibration)
        y .= _apply_calibration!(result, p.kwargs[:calibration], p, y)
    end

    return x,y,z
end

get_plot_data(result::ProspectResults, ::BinnedSFH, ::Union{typeof(CairoMakie.lines!),typeof(CairoMakie.stairs!)}) = get_sfh(result)

#plot_data = get_plot_data(result, plot_spec, plot_func)

# =====================================
# === Rendering functions ===
# =====================================

"Renders a specific plot on an axis using data from the result."
function render!(ax::CairoMakie.Axis, result::ProspectResults, plot::Union{AbstractProspectPlot, AbstractErrors}) 
    # 1. Get the actual plotting function (even if customized)
    plot_func = get_plot_function(plot)

    # 2. assume no want to plot errors
    data = get_plot_data(result, plot, plot_func)

    # 4. Perform the plotting
    plot_func(ax, data... ; plot.kwargs...)

    return ax
end

function render!(ax::CairoMakie.Axis, result::ProspectResults, cmd::PlotCommand)
    # 1- renderizzo i multipliers
    if !isempty(cmd.multipliers)
        for m in cmd.multipliers
            render!(ax, result, m)
        end
    end
    render!(ax, result, cmd.base_plot)
end

"Renders a specific plot on an axis using data from the result."
function render!(ax::CairoMakie.Axis, result::ProspectResults, plot::MaskedSpecPlot)
    # 1. Get the actual plotting function (even if customized)
    plot_func = get_plot_function(plot)

    # 2. assume no want to plot errors
    wave, flux, mask = get_plot_data(result, plot, plot_func)

    if isnothing(wave) || isnothing(flux) || isnothing(mask)
        @warn "One between wave, flux, mask is nothing"
        return nothing
    end


    seg = get_segments(wave, flux, mask)


    # # lines!(ax, wave, flux, ;plot.kwargs)
    for (xs, ys, is_true) in seg

        if is_true
            plot_func(ax, xs, ys; plot.kwargs...)
        else
            plot_func(ax, xs, ys; plot.mask...)
        end
    end
   
    return ax
end


function render!(layout, result::ProspectResults, plot::DensityPlot)
    # Estrai i nomi delle colonne da plottare, escludendo quelli ignorati
    plot_func = get_plot_function(plot)
    all_labels = names(result.chain)
    labels = filter(x -> x ∉ plot.ignored, all_labels)
    n = length(labels)

    if n == 0
        @warn "No parameters to plot. All excluded by `plot.ignored`."
        return
    end

    # Layout dinamico griglia
    n_rows = ceil(Int, sqrt(n))
    n_cols = ceil(Int, n / n_rows)

    idx = 1
    for i in 1:n_rows, j in 1:n_cols
        idx > n && break
        label = labels[idx]
        data = result.chain[!, label]

        if isempty(data)
            @warn "No data for parameter $label"
            idx += 1
            continue
        end

        # Formattazione asse X
        use_sci = maximum(data) > 1e3
        form = use_sci ? "{:.1e}" : "{:.1f}"
        xticklabelrotation = use_sci ? π / 4 : 0

        # Scala logaritmica opzionale (basata su range)
        range = extrema(data) |> x -> (x[1] - x[2])
        xscale = range > 1e6 ? log10 : identity

        # Calcolo quantili
        q16, q50, q84 = quantile(data, [0.16, 0.5, 0.84])
        title = margin_confidence_default_formatter(q50 - q16, q50, q84 - q50, "$label | ")


        ax = Axis(layout[i, j];
            title=title,
            xgridvisible=false,
            ygridvisible=false,
            xtickformat=form,
            xticklabelrotation=xticklabelrotation,
            xscale=xscale
        )

        hidespines!(ax, :t, :r, :l)
        hideydecorations!(ax)

        
        plot_func(ax, data;
            plot.kwargs...  # scatter, linewidth, colormap, etc.
        )


        idx += 1
    end
    return layout
end



function render!(ax::CairoMakie.Axis, result::ProspectResults, plot::BinnedSFH)
    plot_func = get_plot_function(plot)
    plot_data = get_plot_data(result, plot, plot_func)

    plot_func(ax, plot_data...; plot.kwargs...)
    
    return ax
end


function render!(ax::CairoMakie.Axis, result::ProspectResults, plot::SFHErrors)
    z = get_z(result)

    idx = findall(x -> occursin.("logsfr_ratios", x), names(result.chain))
    nbins = length(idx) + 1 
    agebins = zred_to_agebins(z, nbins)

    # calcol errori
    labels = names(result.chain)[findall(x -> occursin.("logsfr_ratios", x), names(result.chain))]
    err_low = []
    err_upp = []
    for name in labels
        low, med, upp = quantile(result.chain[!, name], [0.16, 0.5, 0.84])
        push!(err_low, med - low)
        push!(err_upp, upp - med)
    end
    # force the bins
    pushfirst!(err_low, copy(err_low[begin]))
    pushfirst!(err_upp, copy(err_upp[begin]))

    push!(err_low, copy(err_low[end]))
    push!(err_upp, copy(err_upp[end]))

    reverse!(err_low)
    reverse!(err_upp)

    idx_logmass = findall(x -> x == "logmass", result.sampling["theta_labels"])
    logmass = result.bestfit["parameter"][idx_logmass][1]

    sfr_ratios = [median(result.chain[!, name]) for name in labels]

    # Calcola la massa per ogni bin di età
    mass = SFHUtils.logmass_to_masses(logmass, sfr_ratios, agebins)

    # Estrai i valori di età iniziali e finali dai bin
    initial_bins = [b[1] for b in agebins]
    final_bins = [b[2] for b in agebins]

    # Calcola la SFR usando la massa e i bin di età
    sfr = mass ./ (10.0 .^ final_bins .- 10.0 .^ initial_bins)

    # Crea l'array di età in giga anni
    xx = 10 .^ initial_bins[1] ./ 1e9 .+ 10 .^ final_bins ./ 1e9
    pushfirst!(xx, 1e-9)

    # Calcola il logaritmo della SFR
    yy = log10.(vcat(sfr[1], sfr))

    # little hack
    s = stairs!(ax, xx, yy; color=:transparent)


    _plot_band_error(s, err_upp, err_low, ax; plot.kwargs...)

end




function margin_confidence_default_formatter(low, mid, high, label::String)
    largest_error = max(abs(high), abs(low))
    # Fallback for series with no variance
    if largest_error == 0
        if mid == 0
            digits_after_dot = 0
        else
            digits_after_dot = max(0, 1 - round(Int, log10(abs(mid))))
        end
        title = @sprintf("%.*f", digits_after_dot, mid,)
        return title
    end

    digits_after_dot = max(0, 1 - round(Int, log10(largest_error)))
    use_scientific = digits_after_dot > 4

    if use_scientific
        if round(low, digits=digits_after_dot) == round(high, digits=digits_after_dot)
            title = Makie.rich(label,
                "(",
                @sprintf("%.1f", mid * 10^(digits_after_dot - 1)),
                " ± ",
                @sprintf("%.1f", high * 10^(digits_after_dot - 1)),
                ") × 10",
                Makie.superscript(@sprintf("-%d", digits_after_dot - 1))
            )
        else
            title = Makie.rich(label,
                "( ",
                @sprintf("%.1f", mid * 10^(digits_after_dot - 1)),
                Makie.left_subsup(
                    @sprintf("-%.1f", low * 10^(digits_after_dot - 1)),
                    @sprintf("+%.1f", high * 10^(digits_after_dot - 1))
                ),
                " ) × 10",
                Makie.superscript(@sprintf("-%d", digits_after_dot - 1))
            )
        end
    else
        if round(low, digits=digits_after_dot) == round(high, digits=digits_after_dot)
            title = Makie.rich(label,
                @sprintf("%.*f", digits_after_dot, mid,),
                " ± ",
                @sprintf("%.*f", digits_after_dot, high),
            )
        else
            title = Makie.rich(label,
                @sprintf("%.*f", digits_after_dot, mid),
                Makie.left_subsup(
                    @sprintf("-%.*f", digits_after_dot, low),
                    @sprintf("+%.*f", digits_after_dot, high),
                )
            )
        end
    end

    return title
end


function stairpts(s)#s = stariplot
    pts = s.plots[1].converted[1][]
    [p[1] for p in pts], [p[2] for p in pts]
end

function _plot_band_error(s, errupp, errlow, ax; kwargs...)
    bins, flux_bins = stairpts(s)
    k = 1
    err_upper = errupp[1]
    err_lower = errlow[1]
    for i in 1:length(bins)-1
        if !iseven(i)
            err_upper = errupp[k]
            err_lower = errlow[k]
            k += 1
        end
        lw = [flux_bins[i] - err_lower, flux_bins[i] - err_lower]
        up = [flux_bins[i] + err_upper, flux_bins[i] + err_upper]

        band!(ax, [bins[i], bins[i+1]], lw, up; kwargs...)
    end
end


end #MODULE


"
obs => (asse , SedPlot() |> Errors() |> Mask())

Errors() -
    - AbstractSedPlot   =   SedError()
    - AbstractPhotoPlot =   PhotError()


"

"""
    render!(layout::GridLayout, result::ProspectResults, plot::DensityPlot)

Render density plots for selected parameters from the posterior chain into a GridLayout.
"""