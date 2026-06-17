using CairoMakie
using DataFrames
using Printf
using Statistics
using StatsBase

using ProspectorUtils

const DEFAULT_INPUT = "/home/matteo/Dottorato/DATA/DAWN_Archive/prospector/results/abell2744-castellano1-v4_prism-clear_3073_18819/abell2744-castellano1-v4_prism-clear_3073_18819.h5"
const DEFAULT_OUTPUT = joinpath(@__DIR__, "sfh_branch_comparison.png")
const DEFAULT_REPORT = joinpath(@__DIR__, "sfh_branch_comparison.md")

bin_width(bin) = 10.0^bin[2] - 10.0^bin[1]

function zred_agebins_main(p)
    return zred_to_agebins(get_z(p), n_bins(p))
end

function lookback_from_agebins(agebins)
    offset = 10.0^first(agebins[1]) / 1e9
    return vcat(1e-9, offset .+ 10.0 .^ last.(agebins) ./ 1e9)
end

function sfh_from_components(logmass, ratios, agebins)
    mass = logmass_to_masses(logmass, ratios, agebins)
    sfr = mass ./ map(bin_width, agebins)
    map!(x -> (isfinite(x) && x > 0) ? x : 1e-30, sfr, sfr)
    return lookback_from_agebins(agebins), log10.(vcat(sfr[1], sfr)), mass, sfr
end

function main_canonical_sfh(p)
    ratio_labels = filter(n -> occursin("logsfr_ratios", n), names(p.chain))
    agebins = zred_agebins_main(p)
    logmass = median(skipmissing(p.chain[!, "logmass"]))
    ratios = [median(skipmissing(p.chain[!, name])) for name in ratio_labels]
    lookback, sfh, mass, sfr = sfh_from_components(logmass, ratios, agebins)
    return (; label="main canonical: unweighted median + zred agebins",
            lookback, sfh, mass, sfr, agebins, logmass, ratios)
end

function branch_sfh(p, est, label)
    lookback, sfh, mass = get_sfh(p, est; return_mass=true)
    agebins = get_agebins(p)
    sfr = mass ./ map(bin_width, agebins)
    ratio_labels = filter(n -> occursin("logsfr_ratios", n), labels(p))
    return (; label, lookback, sfh, mass, sfr, agebins,
            logmass=estimate(p, "logmass", est),
            ratios=[estimate(p, name, est) for name in ratio_labels])
end

function max_abs_diff(a, b)
    length(a) == length(b) || return NaN
    return maximum(abs.(a .- b))
end

function format_vec(v)
    return "[" * join((@sprintf("%.6e", x) for x in v), ", ") * "]"
end

function write_report(path, input, p, curves)
    stored = get(p.bestfit, "agebins", nothing)
    agebins_current = get_agebins(p)
    agebins_main = zred_agebins_main(p)
    open(path, "w") do io
        println(io, "# SFH branch comparison")
        println(io)
        println(io, "- input: `$(input)`")
        println(io, "- chain rows: `$(nrow(p.chain))`")
        println(io, "- chain parameters: `$(length(labels(p)))`")
        println(io, "- has weights: `$(has_weights(p))`")
        println(io, "- default current estimator: `$(nameof(typeof(default_estimator(p))))`")
        println(io, "- redshift: `$(get_z(p))`")
        println(io, "- stored bestfit agebins available but ignored: `$(stored isa AbstractMatrix)`")
        println(io, "- max abs agebin diff current-vs-main zred: `$(max_abs_diff(reduce(vcat, collect.(agebins_current)), reduce(vcat, collect.(agebins_main))))`")
        println(io)
        println(io, "| Method | logmass input | max abs SFH diff vs current default | max abs SFH diff vs main canonical | total formed mass | recent log10 SFR |")
        println(io, "|---|---:|---:|---:|---:|---:|")
        default_curve = curves[1]
        main_curve = curves[end]
        for curve in curves
            println(io, "| $(curve.label) | $(@sprintf("%.8f", curve.logmass)) | $(@sprintf("%.6e", max_abs_diff(curve.sfh, default_curve.sfh))) | $(@sprintf("%.6e", max_abs_diff(curve.sfh, main_curve.sfh))) | $(@sprintf("%.6e", sum(curve.mass))) | $(@sprintf("%.6f", curve.sfh[1])) |")
        end
        println(io)
        println(io, "## Ratio Inputs")
        println(io)
        for curve in curves
            println(io, "- $(curve.label): `$(format_vec(curve.ratios))`")
        end
    end
end

function make_plot(output, curves)
    fig = Figure(size=(1180, 760), fontsize=20)
    ax = Axis(fig[1, 1],
        xlabel="Lookback time [Gyr]",
        ylabel="log10 SFR [Msun yr^-1]",
        title="SFH method comparison",
        xscale=log10,
        xgridvisible=true,
        ygridvisible=true)

    palette = (:black, :dodgerblue3, :orange2, :seagreen4)
    styles = (:solid, :dash, :dot, :dashdot)
    widths = (4, 3, 3, 3)
    for (i, curve) in enumerate(curves)
        stairs!(ax, curve.lookback, curve.sfh;
            label=curve.label,
            color=palette[i],
            linestyle=styles[i],
            linewidth=widths[i],
            step=:post)
    end
    axislegend(ax, position=:lb, framevisible=true)
    save(output, fig, px_per_unit=2)
    return output
end

function main()
    input = get(ARGS, 1, DEFAULT_INPUT)
    output = get(ARGS, 2, DEFAULT_OUTPUT)
    report = get(ARGS, 3, DEFAULT_REPORT)

    p = ProspectResults(input; verbose=false)
    default_est = default_estimator(p)
    default_name = nameof(typeof(default_est))
    curves = [
        branch_sfh(p, default_est, "current default: $(default_name) + zred agebins"),
        branch_sfh(p, Median(), "current Median + zred agebins"),
        branch_sfh(p, BestFit(), "current BestFit + zred agebins"),
        main_canonical_sfh(p),
    ]

    make_plot(output, curves)
    write_report(report, input, p, curves)
    println(output)
    println(report)
end

main()
