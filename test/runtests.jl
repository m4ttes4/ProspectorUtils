using Test
using DataFrames
using ProspectorUtils

# -----------------------------------------------------------------------------
# Synthetic in-memory result (no HDF5 I/O)
# -----------------------------------------------------------------------------
function make_result(; weights=true)
    theta = ["logmass", "logsfr_ratios_1", "logsfr_ratios_2", "logsfr_ratios_3"]
    chain = DataFrame(
        logmass         = [9.0, 9.2, 9.4, 9.6, 9.8],
        logsfr_ratios_1 = [0.0, 0.1, -0.1, 0.2, -0.2],
        logsfr_ratios_2 = [0.0, -0.1, 0.1, -0.2, 0.2],
        logsfr_ratios_3 = [0.1, 0.0, -0.1, 0.1, 0.0],
    )
    weights && (chain.weights = [1.0, 1.0, 1.0, 1.0, 1.0])
    bestfit = Dict{String,Any}(
        "parameter" => [9.5, 0.05, -0.05, 0.02],   # aligned to theta order
        "mfrac"     => 0.5,
        "agebins"   => [0.0 8.0; 8.0 8.5; 8.5 9.0; 9.0 13.7],
    )
    sampling = Dict{String,Any}("theta_labels" => theta)
    runparams = Dict{String,Any}("redshift" => 5.0)
    return ProspectResults(chain, runparams, bestfit, sampling, Dict{String,Any}())
end

@testset "ProspectorUtils" begin

    @testset "estimators" begin
        p = make_result()
        @test has_weights(p)
        @test estimate(p, "logmass", BestFit()) == 9.5
        @test estimate(p, "logmass", Median()) ≈ 9.4
        @test estimate(p, "logmass", WeightedMedian()) ≈ 9.4

        pnw = make_result(weights=false)
        @test !has_weights(pnw)
        @test_throws ArgumentError estimate(pnw, "logmass", WeightedMedian())
    end

    @testset "get_mass scale fix" begin
        p = make_result()
        # log10 surviving mass = logmass + log10(mfrac), NOT logmass - mfrac
        @test get_mass(p, BestFit()) ≈ 9.5 + log10(0.5)
        @test get_mass(p, Median()) ≈ 9.4 + log10(0.5)
    end

    @testset "quantiles" begin
        p = make_result()
        q = quantiles(p, "logmass"; weighted=false)
        @test length(q) == 3
        @test issorted(q)
        @test q[2] ≈ 9.4                      # equal weights → median
        @test quantiles(p, "logmass")[2] ≈ 9.4  # weighted default
        @test labels(p) == ["logmass", "logsfr_ratios_1", "logsfr_ratios_2", "logsfr_ratios_3"]  # excludes weights
    end

    @testset "logmass_to_masses" begin
        agebins = [(0.0, 8.0), (8.0, 8.5), (8.5, 9.0)]
        # total mass invariant: sum of per-bin masses = 10^logmass
        m = logmass_to_masses(9.0, [0.0, 0.0], agebins)
        @test length(m) == 3
        @test sum(m) ≈ 10.0^9.0
        @test_throws DimensionMismatch logmass_to_masses(9.0, [0.0], agebins)
        # estimator path reconstructs agebins from redshift.
        p = make_result()
        @test sum(logmass_to_masses(p, BestFit())) ≈ 10.0^9.5
    end

    @testset "build_agebins" begin
        ab = build_agebins(tuniv=13.7, nbins=7)
        @test length(ab) == 7
        @test all(b -> b[1] < b[2], ab)
        @test issorted(first.(ab))
        @test_throws ArgumentError build_agebins(nbins=3)
    end

    @testset "sfh_lookback" begin
        agebins = [(0.0, 8.0), (8.0, 8.5), (8.5, 9.0)]
        lb = sfh_lookback(agebins)
        @test length(lb) == length(agebins) + 1
        @test lb[1] ≈ 1e-9
        @test issorted(lb)
    end

    @testset "get_agebins reconstructs from redshift" begin
        # Stored bestfit/agebins can be inconsistent with the object redshift.
        # get_agebins must ignore them and rebuild Prospector bins from zred.
        theta = ["logmass", "logsfr_ratios_1", "logsfr_ratios_2", "logsfr_ratios_3"]
        chain = DataFrame(logmass=[9.0], logsfr_ratios_1=[0.0], logsfr_ratios_2=[0.0], logsfr_ratios_3=[0.0])
        sampling = Dict{String,Any}("theta_labels" => theta)
        bf = Dict{String,Any}("agebins" => [0.0 8.0; 8.0 8.5; 8.5 9.0; 9.0 13.7])
        p = ProspectResults(chain, Dict{String,Any}("redshift" => 5.0), bf, sampling, Dict{String,Any}())
        @test get_agebins(p) == zred_to_agebins(5.0, 4)
    end

    @testset "get_sfh / get_sfr" begin
        p = make_result()
        lookback, sfh = get_sfh(p)
        @test length(lookback) == length(sfh) == n_bins(p) + 1
        @test all(isfinite, sfh)
        lookback, sfh, mass = get_sfh(p; return_mass=true)
        @test length(mass) == n_bins(p)
        @test isfinite(get_sfr(p; nbins=2))
        # get_sfr must skip sfh[1] (duplicated youngest-bin stairs point)
        _, sfh_full = get_sfh(p)
        @test get_sfr(p; nbins=2) ≈ (sfh_full[2] + sfh_full[3]) / 2
    end

    @testset "maggies2μJy" begin
        @test maggies2μJy(1.0) ≈ 3631e6
        @test maggies2μJy(nothing) === nothing
    end

    @testset "n_bins" begin
        @test n_bins(make_result()) == 4
    end
end

@testset "Aqua" begin
    using Aqua
    Aqua.test_all(ProspectorUtils; ambiguities=false)
end
