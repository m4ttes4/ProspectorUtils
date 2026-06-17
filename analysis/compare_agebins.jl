using Printf

using ProspectorUtils

const DEFAULT_INPUT = "/home/matteo/Dottorato/DATA/DAWN_Archive/prospector/results/abell2744-castellano1-v4_prism-clear_3073_18819/abell2744-castellano1-v4_prism-clear_3073_18819.h5"
const DEFAULT_OUTPUT = joinpath(@__DIR__, "agebins_comparison.md")

bin_width_yr(bin) = 10.0^bin[2] - 10.0^bin[1]
bin_edge_gyr(bin) = 10.0^bin[2] / 1e9

function agebins_to_tuples(m::AbstractMatrix)
    if size(m, 2) == 2
        return [(Float64(m[i, 1]), Float64(m[i, 2])) for i in axes(m, 1)]
    elseif size(m, 1) == 2
        return [(Float64(m[1, j]), Float64(m[2, j])) for j in axes(m, 2)]
    end
    throw(ArgumentError("unexpected agebins shape $(size(m))"))
end

function flat_agebins(agebins)
    return reduce(vcat, collect.(agebins))
end

function max_abs_diff(a, b)
    length(a) == length(b) || return NaN
    return maximum(abs.(a .- b))
end

function write_report(path, input)
    p = ProspectResults(input; verbose=false)
    main_agebins = zred_to_agebins(get_z(p), n_bins(p))
    current_agebins = get_agebins(p)
    stored = get(p.bestfit, "agebins", nothing)
    stored_agebins = stored isa AbstractMatrix ? agebins_to_tuples(stored) : nothing

    open(path, "w") do io
        println(io, "# Agebins comparison")
        println(io)
        println(io, "- input: `$(input)`")
        println(io, "- redshift: `$(get_z(p))`")
        println(io, "- nbins: `$(n_bins(p))`")
        println(io, "- stored `bestfit/agebins` present but ignored: `$(stored isa AbstractMatrix)`")
        println(io, "- max abs diff current-vs-main in log10(years): `$(max_abs_diff(flat_agebins(main_agebins), flat_agebins(current_agebins)))`")
        if stored_agebins !== nothing
            println(io, "- max abs diff stored-vs-main in log10(years): `$(max_abs_diff(flat_agebins(main_agebins), flat_agebins(stored_agebins)))`")
        end
        println(io)
        println(io, "| bin | main low | main high | current low | current high | delta low | delta high | main width [Gyr] | current width [Gyr] | width ratio current/main |")
        println(io, "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")

        for i in eachindex(main_agebins)
            m = main_agebins[i]
            c = current_agebins[i]
            m_width = bin_width_yr(m) / 1e9
            c_width = bin_width_yr(c) / 1e9
            println(io,
                "| ", i,
                " | ", @sprintf("%.8f", m[1]),
                " | ", @sprintf("%.8f", m[2]),
                " | ", @sprintf("%.8f", c[1]),
                " | ", @sprintf("%.8f", c[2]),
                " | ", @sprintf("%+.8f", c[1] - m[1]),
                " | ", @sprintf("%+.8f", c[2] - m[2]),
                " | ", @sprintf("%.8e", m_width),
                " | ", @sprintf("%.8e", c_width),
                " | ", @sprintf("%.6f", c_width / m_width),
                " |")
        end

        if stored_agebins !== nothing
            println(io)
            println(io, "## Stored bestfit agebins ignored by current code")
            println(io)
            println(io, "| bin | main high [Gyr] | current high [Gyr] | stored high [Gyr] | stored-main delta [Gyr] |")
            println(io, "|---:|---:|---:|---:|---:|")
            for i in eachindex(main_agebins)
                m_edge = bin_edge_gyr(main_agebins[i])
                c_edge = bin_edge_gyr(current_agebins[i])
                s_edge = bin_edge_gyr(stored_agebins[i])
                println(io,
                    "| ", i,
                    " | ", @sprintf("%.8e", m_edge),
                    " | ", @sprintf("%.8e", c_edge),
                    " | ", @sprintf("%.8e", s_edge),
                    " | ", @sprintf("%+.8e", s_edge - m_edge),
                    " |")
            end
        end

        println(io)
        println(io, "## Upper edges")
        println(io)
        println(io, "| edge after bin | main [Gyr] | current [Gyr] | delta [Gyr] |")
        println(io, "|---:|---:|---:|---:|")
        for i in eachindex(main_agebins)
            m_edge = bin_edge_gyr(main_agebins[i])
            c_edge = bin_edge_gyr(current_agebins[i])
            println(io,
                "| ", i,
                " | ", @sprintf("%.8e", m_edge),
                " | ", @sprintf("%.8e", c_edge),
                " | ", @sprintf("%+.8e", c_edge - m_edge),
                " |")
        end
    end
    return path
end

input = get(ARGS, 1, DEFAULT_INPUT)
output = get(ARGS, 2, DEFAULT_OUTPUT)
println(write_report(output, input))
