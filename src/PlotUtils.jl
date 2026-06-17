module PlotUtils

# Plot *data model* only. All CairoMakie-dependent behaviour (get_plot_function,
# get_plot_data, render!) lives in ext/ProspectorUtilsCairoMakieExt.jl and is
# loaded automatically once the user does `using CairoMakie`.

using ..DataUtils

export AbstractProspectPlot, AbstractSedPlot, AbstractPhotoPlot, AbstractMultiplier, AbstractErrors
export ObsSpecPlot, BFSpecPlot, BFSedPlot, DensityPlot, ObsPhotoPlot, BFPhotoPlot, BinnedSFH, LineIndicators
export MaskedSpecPlot, MaskedPhotoPlot
export SedErrors, PhotoErrors, SFHErrors, Errors, Mask, Calibration, PlotCommand
export get_zorder, get_plot_function, get_error_plot_function, get_plot_data, render!

# Functions whose methods are provided by the CairoMakie extension.
function get_plot_function end
function get_error_plot_function end
function get_plot_data end
function render! end

# =============================================================================
# Plot types
# =============================================================================

abstract type AbstractProspectPlot end
abstract type AbstractSedPlot <: AbstractProspectPlot end
abstract type AbstractPhotoPlot <: AbstractProspectPlot end
abstract type AbstractMultiplier end
abstract type AbstractErrors <: AbstractMultiplier end

# Simple kwargs-carrying plots generated with identical boilerplate.
for (T, super) in ((:ObsSpecPlot, :AbstractSedPlot), (:BFSpecPlot, :AbstractSedPlot),
                   (:BFSedPlot, :AbstractSedPlot), (:ObsPhotoPlot, :AbstractPhotoPlot),
                   (:BFPhotoPlot, :AbstractPhotoPlot), (:BinnedSFH, :AbstractProspectPlot),
                   (:SedErrors, :AbstractErrors), (:PhotoErrors, :AbstractErrors),
                   (:SFHErrors, :AbstractErrors), (:Errors, :AbstractMultiplier),
                   (:Mask, :AbstractMultiplier))
    @eval begin
        mutable struct $T <: $super
            kwargs::Dict{Symbol,Any}
        end
        $T(; kwargs...) = $T(Dict{Symbol,Any}(kwargs))
    end
end

mutable struct MaskedSpecPlot <: AbstractSedPlot
    kwargs::Dict{Symbol,Any}
    mask::Dict{Symbol,Any}
end

mutable struct MaskedPhotoPlot <: AbstractPhotoPlot
    kwargs::Dict{Symbol,Any}
    mask::Dict{Symbol,Any}
end

mutable struct DensityPlot <: AbstractProspectPlot
    ignored::Vector{String}
    kwargs::Dict{Symbol,Any}
end
DensityPlot(ignored::Vector{String}=String[]; kwargs...) = DensityPlot(ignored, Dict{Symbol,Any}(kwargs))

mutable struct LineIndicators <: AbstractProspectPlot
    lines::Dict{String,Float64}
    kwargs::Dict{Symbol,Any}
end
LineIndicators(lines; kwargs...) = LineIndicators(lines, Dict{Symbol,Any}(kwargs))

# Render order (lower drawn first). No CairoMakie dependency.
get_zorder(::AbstractProspectPlot) = 10
get_zorder(::AbstractErrors) = 5
get_zorder(::LineIndicators) = 15

# =============================================================================
# Multipliers: errors, masks, calibration
# =============================================================================

_errors(::AbstractSedPlot, kwargs) = SedErrors(kwargs)
_errors(::AbstractPhotoPlot, kwargs) = PhotoErrors(kwargs)
_errors(::BinnedSFH, kwargs) = SFHErrors(kwargs)
(self::Errors)(p::AbstractProspectPlot) = _errors(p, self.kwargs)

struct Calibration
    calibration::Union{String,Vector{<:Real}}
end
Calibration() = Calibration("speccal")

struct PlotCommand
    base_plot::AbstractProspectPlot
    multipliers::Vector{<:AbstractMultiplier}
end

# Pipe DSL: attach calibration / mask / errors to a base plot.
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

Base.:|>(left::AbstractSedPlot, right::Mask) = MaskedSpecPlot(left.kwargs, right.kwargs)
Base.:|>(left::AbstractPhotoPlot, right::Mask) = MaskedPhotoPlot(left.kwargs, right.kwargs)

function Base.:|>(left::AbstractProspectPlot, right::AbstractMultiplier)
    return PlotCommand(left, [right(left)])
end

function Base.:|>(left::PlotCommand, right::AbstractMultiplier)
    return PlotCommand(left.base_plot, [left.multipliers..., right(left.base_plot)])
end

end # module
