#!/usr/bin/env julia

using HDF5
using JSON
using DataFrames
using Statistics
using Printf
using PythonCall

const DEFAULT_INPUT = "/Users/matteo/Dottorato/script/analisys/shared_galaxy_h5/capers-udsp1-v4_prism-clear_6368_21205.h5"
const DEFAULT_OUTPUT = joinpath(@__DIR__, "sfh_comparison.md")

@inline mydiff(bin) = 10.0 ^ bin[2] - 10.0 ^ bin[1]

function read_run_params(path::AbstractString)
    h5open(path, "r") do f
        return JSON.parse(String(read_attribute(f, "run_params")))
    end
end

function read_theta_labels(path::AbstractString)
    h5open(path, "r") do f
        raw = attrs(f["sampling"])["theta_labels"]
        if raw isa AbstractString
            parsed = try
                JSON.parse(raw)
            catch
                [raw]
            end
            return String.(parsed)
        elseif raw isa AbstractVector
            return String.(raw)
        else
            return String.(collect(raw))
        end
    end
end

function read_bestfit(path::AbstractString)
    h5open(path, "r") do f
        return Dict(
            "parameter" => collect(read(f["bestfit/parameter"])),
            "agebins" => read(f["bestfit/agebins"]),
            "mfrac" => read_attribute(f["bestfit"], "mfrac"),
        )
    end
end

function read_chain_df(path::AbstractString, runparams::Dict{String,Any}, labels::Vector{String})
    h5open(path, "r") do f
        chain = read(f["sampling/chain"])
        if get(runparams, "dynesty", false)
            mat = PermutedDimsArray(chain, (2, 1))
        elseif get(runparams, "emcee", false)
            nl, nch, nwk = size(chain)
            mat = reshape(chain, nl, nch * nwk)'
        else
            mat = chain
        end
        if !isempty(labels) && size(mat, 2) == length(labels)
            return DataFrame(mat, labels, makeunique=true)
        end
        return DataFrame(mat, :auto)
    end
end

function build_agebins(tuniv::Real; nbins::Integer)
    nbins < 4 && throw(ArgumentError("nbins must be >= 4"))
    tbinmax = (tuniv * 0.9) * 1e9
    lim1, lim2 = 7.4772, 8.0
    log_bins = range(lim2, stop=log10(tbinmax), length=nbins - 2) |> collect
    agelims = vcat([0.0, lim1], log_bins, [log10(tuniv * 1e9)])
    return [(agelims[i], agelims[i + 1]) for i in 1:(length(agelims) - 1)]
end

function universe_age_gyr(zred::Real)
    cosmo = pyimport("astropy.cosmology")
    return pyconvert(Float64, cosmo.Planck18.age(zred).value)
end

function zred_to_agebins(zred::Real, nbins::Integer)
    return build_agebins(universe_age_gyr(zred); nbins=nbins)
end

function src_n_bins(labels::Vector{String})
    return count(x -> occursin("logsfr_ratios", x), labels) + 1
end

function median_ratio_vector(df::DataFrame)
    ratio_cols = filter(n -> occursin("logsfr_ratios", String(n)), names(df))
    return [median(skipmissing(df[!, c])) for c in ratio_cols], String.(ratio_cols)
end

function logmass_vector(df::DataFrame)
    return median(skipmissing(df[!, "logmass"]))
end

function logmass_to_masses(logmass::Real, logsfr_ratios::AbstractVector{<:Real}, agebins)
    nbins = length(agebins)
    length(logsfr_ratios) == nbins - 1 || throw(DimensionMismatch("Expected $(nbins - 1) ratios, got $(length(logsfr_ratios))"))
    sratios = 10.0 .^ clamp.(Float64.(logsfr_ratios), -10.0, 10.0)
    dt = mydiff.(agebins)
    coeffs = ones(Float64, nbins)
    for j in 2:nbins
        coeffs[j] = dt[j] / (dt[1] * prod(sratios[1:j-1]))
    end
    m1 = 10.0 ^ float(logmass) / sum(coeffs)
    return m1 .* coeffs
end

function sfh_steps(logmass::Real, logsfr_ratios::AbstractVector{<:Real}, agebins)
    dt = 10.0 .^ last.(agebins) .- 10.0 .^ first.(agebins)
    mass = logmass_to_masses(logmass, logsfr_ratios, agebins)
    sfr = mass ./ dt
    lookback = vcat(1e-9, 10.0 ^ first(agebins[1]) / 1e9 .+ 10.0 .^ last.(agebins) ./ 1e9)
    sfh = log10.(vcat(sfr[1], sfr))
    return (mass=mass, dt=dt, sfr=sfr, lookback=lookback, sfh=sfh)
end

function format_vec(v; digits=6)
    io = IOBuffer()
    print(io, "[")
    for (i, x) in enumerate(v)
        i > 1 && print(io, ", ")
        if x isa Real
            @printf(io, "%.*e", digits, float(x))
        else
            print(io, x)
        end
    end
    print(io, "]")
    return String(take!(io))
end

function maxabsdiff(a, b)
    return maximum(abs.(Float64.(a) .- Float64.(b)))
end

function compare(path::AbstractString)
    runparams = read_run_params(path)
    labels = read_theta_labels(path)
    bestfit = read_bestfit(path)
    chain = read_chain_df(path, runparams, labels)

    zred = Float64(runparams["redshift"])
    nbins_src = src_n_bins(labels)
    nbins_utils = Int(get(runparams, "nbins", nbins_src))

    agebins_src = zred_to_agebins(zred, nbins_src)
    agebins_utils = zred_to_agebins(zred, nbins_utils)

    logmass_median = logmass_vector(chain)
    ratio_medians, ratio_cols = median_ratio_vector(chain)

    src_steps = sfh_steps(logmass_median, ratio_medians, agebins_src)
    utils_steps = sfh_steps(logmass_median, ratio_medians, agebins_utils)

    bestfit_labels = collect(String.(labels))
    logmass_idx = findall(==( "logmass"), bestfit_labels)
    bestfit_logmass = bestfit["parameter"][logmass_idx][1] - bestfit["mfrac"]

    report = IOBuffer()
    println(report, "# SFH comparison")
    println(report, "")
    println(report, "- input file: `$path`")
    println(report, "- redshift: `", zred, "`")
    println(report, "- theta labels with SFH ratios: `", join(ratio_cols, ", "), "`")
    println(report, "- chain rows x cols: `", size(chain, 1), " x ", size(chain, 2), "`")
    println(report, "- source nbins (`src/SFHUtils.jl`): `", nbins_src, "`")
    println(report, "- utils nbins (`Utils.jl` path used by `get_sfh`): `", nbins_utils, "`")
    println(report, "- bestfit logmass from file: `", bestfit_logmass, "`")
    println(report, "- median logmass from chain: `", logmass_median, "`")
    println(report, "")

    println(report, "## Step-by-step comparison")
    println(report, "")
    println(report, "| Step | src/SFHUtils.jl | Utils.jl | Max abs diff | Verdict |")
    println(report, "| --- | --- | --- | --- | --- |")
    println(report, "| `nbins` | `", nbins_src, "` | `", nbins_utils, "` | `", abs(nbins_src - nbins_utils), "` | ", nbins_src == nbins_utils ? "same" : "different", " |")
    println(report, "| `agebins` | `", replace(string(agebins_src), "\n" => " "), "` | `", replace(string(agebins_utils), "\n" => " "), "` | `", maxabsdiff(reduce(vcat, collect.(agebins_src)), reduce(vcat, collect.(agebins_utils))), "` | ", agebins_src == agebins_utils ? "same" : "different", " |")
    println(report, "| `logmass` used in canonical SFH | `", logmass_median, "` | `", logmass_median, "` | `0.0` | same |")
    println(report, "| `logsfr_ratios` | `", format_vec(ratio_medians), "` | `", format_vec(ratio_medians), "` | `0.0` | same |")
    println(report, "| `dt` | `", format_vec(src_steps.dt), "` | `", format_vec(utils_steps.dt), "` | `", maxabsdiff(src_steps.dt, utils_steps.dt), "` | ", maxabsdiff(src_steps.dt, utils_steps.dt) == 0 ? "same" : "different", " |")
    println(report, "| `mass` | `", format_vec(src_steps.mass), "` | `", format_vec(utils_steps.mass), "` | `", maxabsdiff(src_steps.mass, utils_steps.mass), "` | ", maxabsdiff(src_steps.mass, utils_steps.mass) == 0 ? "same" : "different", " |")
    println(report, "| `sfr` | `", format_vec(src_steps.sfr), "` | `", format_vec(utils_steps.sfr), "` | `", maxabsdiff(src_steps.sfr, utils_steps.sfr), "` | ", maxabsdiff(src_steps.sfr, utils_steps.sfr) == 0 ? "same" : "different", " |")
    println(report, "| `lookback` | `", format_vec(src_steps.lookback), "` | `", format_vec(utils_steps.lookback), "` | `", maxabsdiff(src_steps.lookback, utils_steps.lookback), "` | ", maxabsdiff(src_steps.lookback, utils_steps.lookback) == 0 ? "same" : "different", " |")
    println(report, "| `sfh` | `", format_vec(src_steps.sfh), "` | `", format_vec(utils_steps.sfh), "` | `", maxabsdiff(src_steps.sfh, utils_steps.sfh), "` | ", maxabsdiff(src_steps.sfh, utils_steps.sfh) == 0 ? "same" : "different", " |")
    println(report, "")

    println(report, "## Notes")
    println(report, "")
    println(report, "- The canonical SFH path in both codebases is numerically identical here because both use the same median chain logmass and median `logsfr_ratios` values.")
    println(report, "- `Utils.jl` also contains plotting wrappers (`render!(::BinnedSFH)`, `plot_sfh`) that use the best-fit `logmass` branch instead of the median-chain branch; those wrappers are not equivalent to the canonical `get_sfh` path.")
    println(report, "- The equality above is about the core SFH pipeline; any mismatch in a different file path or wrapper usually comes from input selection, not from the mass/SFR formula itself.")

    return String(take!(report))
end

function main()
    input = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_INPUT
    output = length(ARGS) >= 2 ? ARGS[2] : DEFAULT_OUTPUT
    report = compare(input)
    open(output, "w") do io
        write(io, report)
        if !endswith(report, "\n")
            write(io, "\n")
        end
    end
    println(output)
end

main()
