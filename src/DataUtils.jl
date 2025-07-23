module DataUtils

export ProspectResults, ProspectorBestFit, ProspectorObs, ProspectorObs
export maggies2μJy, get_z,labels,bestfit,get_obs_sflux,get_obs_swave,get_obs_serr
export get_obs_pflux,get_obs_pwave,get_obs_perr,get_bf_sflux,get_bf_swave,get_bf_pflux,get_bf_pwave
export get_bf_cont, get_bf_calib, logmass_to_masses, n_bins



using HDF5
using DataFrames
using JSON
using PythonCall

# using LinearAlgebra
# -------------------------------
#--- Constants for HDF5 groups ---
# -------------------------------
const H5_GROUPS = ("obs", "bestfit", "sampling")

# -------------------------------
#--- Abstract Types & Structs ---
# -------------------------------

abstract type AbstractProspectResult end

"""
Container for full Prospector output
"""
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


# -------------------------------
#--- Low-level HDF5 Reader -------
# -------------------------------

function _read_h5_internal(f::HDF5.File, path::String, ::Val{true})
    # Controlla prima se il path esiste per evitare errori
    if !haskey(f, path)
        return Dict{String,Any}() # Restituisci un dizionario vuoto se il path non esiste
    end

    attributes = attrs(f[path])
    res = Dict{String,Any}()
    for key in keys(attributes)
        tmp = read_attribute(f[path], key)
        # Semplifichiamo il parsing del JSON
        if isa(tmp, AbstractString) || isa(tmp, Vector{UInt8})
            str_val = isa(tmp, AbstractString) ? tmp : String(tmp)
            try
                res[key] = JSON.parse(str_val)
            catch
                res[key] = str_val # fallback se non è JSON valido
            end
        else
            res[key] = tmp
        end
    end
    return res
end

# Metodo per at=false: restituisce SEMPRE il tipo di read()
function _read_h5_internal(f::HDF5.File, path::String, ::Val{false})

    !haskey(f, path) && return nothing

    return read(f[path])
end

function _read_h5(f::HDF5.File, path::String; at::Bool=false)
    return _read_h5_internal(f, path, Val(at))
end


"""
Generic reader for all HDF5 groups in one go
Returns NamedTuple(:obs, :bestfit, :sampling) of Dicts
"""
function _read_all_groups(f::HDF5.File; at::Bool=false)
    data = ntuple(i -> _read_h5(f, H5_GROUPS[i]; at=at), length(H5_GROUPS))
    return (; obs=data[1], bestfit=data[2], sampling=data[3])
end

# --------------------------------
#--- Public API Constructors -----
# --------------------------------
"""
Main constructor: REFACTORED to perform all I/O once and pass pre-loaded data 
to helper functions, ensuring high performance and maintainability.
"""
function ProspectResults(filename::String; verbose::Bool=true)
    # Apri il file una sola volta per tutte le operazioni di lettura
    h5open(filename, "r") do f
        # --- 1. I/O OTTIMIZZATO: Leggi tutti i dati e gli attributi una sola volta ---
        attrs_nt = _read_all_groups(f; at=true)
        data_nt = _read_all_groups(f; at=false)

        # Leggi i parametri di run una sola volta
        run_params_dict = _get_run_params(f)

        # Unisci gli attributi nei rispettivi dizionari di dati
        merge!(data_nt.obs, attrs_nt.obs)
        merge!(data_nt.bestfit, attrs_nt.bestfit)
        merge!(data_nt.sampling, attrs_nt.sampling)

        # Aggiungi riferimenti ai dati di lunghezza d'onda per convenienza
        data_nt.bestfit["wavelength"] = data_nt.obs["wavelength"]
        data_nt.bestfit["phot_wave"] = data_nt.obs["phot_wave"]

        # --- 2. CHIAMATA CORRETTA: Passa solo i dati necessari a `_build_chain_df` ---
        # La funzione riceve i dizionari `sampling`, non l'intero `NamedTuple`.
        df_chain = _build_chain_df(data_nt.sampling, attrs_nt.sampling, run_params_dict; verbose=verbose)


        # Costruisci l'oggetto finale con tutti i dati processati
        return ProspectResults(df_chain,
            run_params_dict, # Usa la struct type-stable
            data_nt.bestfit,
            data_nt.sampling,
            data_nt.obs)
    end
end

# --------------------------------
#--- Accessor Functions ---------
# --------------------------------
get_z(p::ProspectResults) = get(p.runparams, "redshift", nothing)
labels(p::ProspectResults) = names(p.chain)
bestfit(p::ProspectResults) = get(p.bestfit, "parameter", nothing)
n_bins(p::ProspectResults) = sum(occursin.("logsfr_ratios", names(p.chain))) + 1


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


get_mass(p::ProspectResults) = first(p.bestfit["parameter"][findall(x -> x == "logmass", p.sampling["theta_labels"])]) - p.bestfit["mfrac"]
maggies2μJy(mag::Real) = mag * 1e6 * 3631
maggies2μJy(::Nothing) = nothing

function bestfit(p::ProspectResults, param::String)
    return first(p.bestfit["parameter"][findall(x -> x == param, p.sampling["theta_labels"])]) 
end


# --------------------------------
#--- Run Params Reader -----------
# --------------------------------

"""
Read `run_params` attribute and parse JSON
Takes open file handle f
"""
function _get_run_params(f::HDF5.File)::Dict{String,Any}
    try
        raw = read_attribute(f, "run_params")
        j = isa(raw, AbstractString) ? raw : String(raw)
        return JSON.parse(j)
    catch e
        @warn "Failed to parse run_params: $e"
        return Dict{String,Any}()
    end
end

# --------------------------------
#--- Chain DataFrame Builder -----
# --------------------------------


"""
Build DataFrame for MCMC chain from pre-loaded data.
This version avoids redundant I/O and minimizes memory allocations by using views.
"""
function _build_chain_df(sampling_data::Dict,
    sampling_attrs::Dict,
    run_params::Dict;
    verbose::Bool=true)::DataFrame

    labels = get(sampling_attrs, "theta_labels", String[])
    chain = get(sampling_data, "chain", Float64[])

    # Se la chain non esiste o è vuota, restituisci un DataFrame vuoto.
    isempty(chain) && return DataFrame()

    # Anche qui, l'accesso a run_params è più efficiente se si usa una struct tipizzata
    has_dyn = get(run_params, "dynesty", false)
    has_e = get(run_params, "emcee", false)

    if has_dyn && has_e
        @warn "Both dynesty and emcee flags true—defaulting to dynesty behavior"
    end

    local mat
    if has_dyn
        verbose && @info "Building chain from Dynesty sampler"
        # 2. NO ALLOCATION: PermutedDimsArray è una vista, non una copia.
        mat = PermutedDimsArray(chain, (2, 1))
    elseif has_e
        verbose && @info "Building chain from Emcee sampler"
        nl, nch, nwk = size(chain)
        # 2. LOW ALLOCATION: reshape e ' creano viste. La materializzazione 
        # avviene al momento della creazione del DataFrame, ma abbiamo evitato la prima copia.
        mat = reshape(chain, nl, nch * nwk)'
    else
        mat = chain
    end

    # Controllo di sicurezza per le etichette
    if !isempty(labels) && size(mat, 2) != length(labels)
        @warn "Mismatch between number of columns ($(size(mat, 2))) and labels ($(length(labels))). Using generic labels."
        return DataFrame(mat, :auto)
    end

    # L'opzione makeunique=true è una buona pratica per evitare errori con etichette duplicate
    return DataFrame(mat, labels, makeunique=true)
end

mydiff(bins) = 10^bins[2] - 10^bins[1]

"""
    logmass_to_masses(logmass, logsfr_ratios, agebins) -> Vector{Float64}

Converts a value of log₁₀(∑ᵢ Mᵢ) and an array of log₁₀(SFR_j / SFR₍ⱼ₊₁₎) into Mᵢ values.

## Arguments
- `logmass::Real`: The log₁₀ value of the total mass ∑ Mᵢ.
- `logsfr_ratios::Vector{Real}`: Vector of size (nbins-1) containing the log₁₀ of SFR ratios.
- `agebins::Vector{Vector{Real}}`: Matrix of size (nbins, 2) with the age bin limits in log₁₀(years).

## Returns
- `Vector{Float64}`: An array containing the Mᵢ values.

## Notes
- Assumes that j=0 (in Python) corresponds to the most recent bin.
- This function follows the behavior of `prospector`.
- Assumes that `logmass` is the median value from the sampling chain rather than the best-fit value.
"""
function logmass_to_masses(logmass::T, logsfr_ratios::Vector{T}, agebins::Vector{Tuple{T,T}}) where {T<:Real}
    nbins = size(agebins, 1)
    sratios = 10 .^ clamp.(logsfr_ratios, -10, 10)
    dt = mydiff.(agebins)
    coeffs = ones(nbins)
    for j in 2:nbins
        coeffs[j] = dt[j] / (dt[1] * prod(sratios[1:j-1]))
    end
    m1 = 10^logmass ./ sum(coeffs)
    return m1 .* coeffs
end


function Base.show(io::IO, pr::ProspectResults)
    println(io, "ProspectResults Summary")
    println(io, "────────────────────────")
    println(io, "Chain dataframe:        ", size(pr.chain, 1), " rows × ", size(pr.chain, 2), " columns")
    println(io, "Run parameters:         ", length(pr.runparams), " entries")
    println(io, "Bestfit parameters:     ", length(pr.bestfit), " entries")
    println(io, "Sampling info:          ", length(pr.sampling), " entries")
    println(io, "Observational data:     ", length(pr.obs), " entries")
end

function Base.show(io::IO, po::ProspectorObs)
    println(io, "ProspectorObs")
    println(io, "──────────────")
    println(io, "Entries: ", length(po.obs))
    println(io, "Keys:    ", join(keys(po.obs), ", "))
end

function Base.show(io::IO, pb::ProspectorBestFit)
    println(io, "ProspectorBestFit")
    println(io, "──────────────────")
    println(io, "Entries: ", length(pb.bestfit))
    println(io, "Keys:    ", join(keys(pb.bestfit), ", "))
end

function Base.show(io::IO, ps::ProspectorSampling)
    println(io, "ProspectorSampling")
    println(io, "───────────────────")
    println(io, "Entries: ", length(ps.sampling))
    println(io, "Keys:    ", join(keys(ps.sampling), ", "))
end
    
end