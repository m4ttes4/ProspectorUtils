module PlotUtils
using ..DataUtils
using ..SFHUtils

using CairoMakie
# import Base.|>
export AbstractProspectPlot, AbstractSedPlot, AbstractPhotoPlot, AbstractMultiplier
export ObsSpecPlot, BFSedPlot, DensityPlot, ObsPhotoPlot, BFPhotoPlot
export Errors, MaskPlot, Calibration
export PlotCommand
export get_plot_function, get_error_plot_function
export get_plot_data#, validate_plot_data
export render!



abstract type AbstractProspectPlot end
abstract type AbstractSedPlot <: AbstractProspectPlot end
abstract type AbstractPhotoPlot <: AbstractProspectPlot end

abstract type AbstractMultiplier{T<:Union{AbstractSedPlot,AbstractPhotoPlot}} end
# =====================================
# === SED Plot Structs ===
# ====================================="ObsSpecPlot(; kwargs...)\n\nSED plot based on observed data."


# --- Plot Structs Senza plotfunc ---
"ObsSpecPlot(; kwargs...)\n\nSED plot based on observed data."
struct ObsSpecPlot <: AbstractSedPlot
    kwargs::Dict{Symbol,Any}
end
ObsSpecPlot(; kwargs...) = ObsSpecPlot(Dict{Symbol,Any}(kwargs))

"BFSedPlot(; kwargs...)\nSED plot based on best-fit model output."
struct BFSedPlot <: AbstractSedPlot
    kwargs::Dict{Symbol,Any}
end
BFSedPlot(; kwargs...) = BFSedPlot(Dict{Symbol,Any}(kwargs))

"DensityPlot(ignored; kwargs...)\n\nDensity plot of best-fit parameter distributions."
struct DensityPlot <: AbstractProspectPlot
    ignored::Vector{String}
    kwargs::Dict{Symbol,Any}
end
DensityPlot(ignored::Vector{String}; kwargs...) = DensityPlot(ignored, Dict{Symbol,Any}(kwargs))

"ObsPhotoPlot(; kwargs...)\n\nPhotometry plot for observed data."
struct ObsPhotoPlot <: AbstractPhotoPlot
    kwargs::Dict{Symbol,Any}
end
ObsPhotoPlot(; kwargs...) = ObsPhotoPlot(Dict{Symbol,Any}(kwargs))

"BFPhotoPlot(; kwargs...)\n\nPhotometry plot for best-fit model predictions."
struct BFPhotoPlot <: AbstractPhotoPlot
    kwargs::Dict{Symbol,Any}
end
BFPhotoPlot(; kwargs...) = BFPhotoPlot(Dict{Symbol,Any}(kwargs))

#

# Multiplier Structs (Ora con un singolo costruttore interno se non hai bisogno di logica complessa)
struct Errors{T<:Union{AbstractSedPlot,AbstractPhotoPlot}} <: AbstractMultiplier{T}
    plotfunc::Function
    kwargs::Dict{Symbol,Any}
    # Costruttore interno per uso del framework, non per l'utente finale
    Errors{T}(plotfunc::Function, kwargs::Dict{Symbol,Any}) where {T<:Union{AbstractSedPlot,AbstractPhotoPlot}} = new(plotfunc, kwargs)
end

struct MaskPlot{T<:Union{AbstractSedPlot,AbstractPhotoPlot}} <: AbstractMultiplier{T}
    plotfunc::Function
    kwargs::Dict{Symbol,Any}
    MaskPlot{T}(plotfunc::Function, kwargs::Dict{Symbol,Any}) where {T<:Union{AbstractSedPlot,AbstractPhotoPlot}} = new(plotfunc, kwargs)
end


# Definisci Errors e MaskPlot come FUNZIONI che creano le istanze parametrizzate
# Questi saranno gli "oggetti" che l'utente userà con il pipe.
function Errors(; kwargs...)
    # Quando questa funzione viene chiamata, non sa ancora a quale tipo di plot si applicherà.
    # Restituisce una "funzione fabbrica" che verrà chiamata dal pipe per creare l'istanza tipizzata.
    return (plot_type) -> begin
        if plot_type <: AbstractSedPlot
            Errors{AbstractSedPlot}(CairoMakie.band!, Dict{Symbol,Any}(kwargs))
        elseif plot_type <: AbstractPhotoPlot
            Errors{AbstractPhotoPlot}(CairoMakie.errorbars!, Dict{Symbol,Any}(kwargs))
        else
            error("Errors multiplier is only applicable to AbstractSedPlot or AbstractPhotoPlot, not $plot_type.")
        end
    end
end

function MaskPlot(; kwargs...)
    return (plot_type) -> begin
        if plot_type <: AbstractSedPlot
            MaskPlot{AbstractSedPlot}(CairoMakie.stairs!, Dict{Symbol,Any}(kwargs))
        elseif plot_type <: AbstractPhotoPlot
            MaskPlot{AbstractPhotoPlot}(CairoMakie.stairs!, Dict{Symbol,Any}(kwargs))
        else
            error("MaskPlot multiplier is only applicable to AbstractSedPlot or AbstractPhotoPlot, not $plot_type.")
        end
    end
end

struct Calibration{T}
    key::Union{Nothing,String}
    poly::Union{Nothing,Vector{T}}
end
Calibration() = Calibration("speccal", nothing)




# =====================================
# ===  PlotCommand ===
# =====================================

struct PlotCommand
    base_plot::AbstractProspectPlot
    multipliers::AbstractMultiplier
end



# 1. Plot base diventa un PlotCommand con il primo modificatore
function Base.:|>(left::AbstractProspectPlot, right_factory::Function)
    # right_factory è la funzione anonima restituita da Errors() o MaskPlot()
    # Esegui la funzione fabbrica per ottenere l'istanza del moltiplicatore tipizzata
    multiplier_instance = right_factory(typeof(left))
    return PlotCommand(left, [multiplier_instance])
end

# 2. PlotCommand esistente aggiunge un altro modificatore
function Base.:|>(left::PlotCommand, right_factory::Function)
    multiplier_instance = right_factory(typeof(left.base_plot))
    push!(left.multipliers, multiplier_instance)
    return left
end


# 2. Funzioni di Dispatch per la Funzione di Plotting
# ----------------------------------------------------

"Restituisce la funzione di plotting da usare per un dato AbstractProspectPlot."
function get_plot_function end # Definiamo la funzione generica

# Metodi di default (questi sono i metodi che l'utente può estendere)
get_plot_function(::ObsSpecPlot) = CairoMakie.stairs!
get_plot_function(::BFSedPlot) = CairoMakie.stairs!
get_plot_function(::ObsPhotoPlot) = CairoMakie.scatter!
get_plot_function(::BFPhotoPlot) = CairoMakie.scatter!
get_plot_function(::DensityPlot) = CairoMakie.density!
get_plot_function(::AbstractMultiplier{AbstractPhotoPlot}) = CairoMakie.errorbars!
get_plot_function(::AbstractMultiplier{AbstractSedPlot}) = CairoMakie.band!


# Funzioni per ottenere la funzione di plot degli errori (come prima)
function get_error_plot_function(p::T) where {T<:AbstractSedPlot}
    return CairoMakie.band!
end
function get_error_plot_function(p::T) where {T<:AbstractPhotoPlot}
    return CairoMakie.errorbars!
end


# =====================================
# === Plot data vaidations ===
# =====================================
function get_plot_data end


# Caso generale: Dati per funzioni di plot (x, y) come scatter!, lines!, stairs!
function get_plot_data(result::ProspectResults, ::ObsPhotoPlot, ::Union{typeof(CairoMakie.scatter!),typeof(CairoMakie.lines!)})

    y = get_obs_pflux(result) |> copy .|> maggies2μJy
    x = get_obs_pwave(result) |> copy .|> maggies2μJy
    return x, y
end

function get_plot_data(result::ProspectResults, ::ObsSpecPlot, ::Union{typeof(CairoMakie.stairs!),typeof(CairoMakie.lines!)})

    y = get_obs_sflux(result) |> copy .|> maggies2μJy
    x = get_obs_swave(result) |> copy .|> maggies2μJy
    return x, y
end


function get_plot_data(result::ProspectResults, ::AbstractMultiplier{AbstractSedPlot}, ::typeof(CairoMakie.band!))

    y = get_obs_sflux(result) |> copy .|> maggies2μJy
    x = get_obs_swave(result) |> copy .|> maggies2μJy
    z = get_obs_serr(result) |> copy .|> maggies2μJy

    return x, y, z
end


function get_plot_data(result::ProspectResults, ::AbstractMultiplier{AbstractPhotoPlot}, ::typeof(CairoMakie.errorbars!))

    y = get_obs_pflux(result) |> copy .|> maggies2μJy
    x = get_obs_pwave(result) |> copy .|> maggies2μJy
    z = get_obs_perr(result) |> copy .|> maggies2μJy

    return x, y, z
end

#plot_data = get_plot_data(result, plot_spec, plot_func)

# =====================================
# === Rendering functions ===
# =====================================

"Renders a specific plot on an axis using data from the result."
function render!(ax::CairoMakie.Axis, result::ProspectResults, plot::Union{AbstractProspectPlot, AbstractMultiplier{T}}) where T
    # 1. Get the actual plotting function (even if customized)
    plot_func = get_plot_function(plot)

    # 2. assume no want to plot errors
    data = get_plot_data(result, plot, plot_func)

    # 3. Validate the data, passing the plot function for specific validation rules
    # if !validate_plot_data(plot_data, plot_func)
    #     @warn "Invalid data for $(typeof(plot)) with function $(plot_func). Plot will not be drawn."
    #     return
    # end

    # 4. Perform the plotting
    plot_func(ax, data... ; plot.kwargs...)

    return nothing
end




end #MODULE