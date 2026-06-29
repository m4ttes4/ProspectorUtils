module ProspectorUtils

include("DataUtils.jl")
using .DataUtils
export AbstractProspectResult, ProspectResults, ProspectorObs, ProspectorBestFit, ProspectorSampling
export Estimator, BestFit, Median, WeightedMedian, has_weights, estimate
export quantiles
export get_mass, maggies2μJy, get_z, labels, bestfit, n_bins
export get_obs_sflux, get_obs_swave, get_obs_serr, get_obs_pflux, get_obs_pwave, get_obs_perr
export get_bf_sflux, get_bf_swave, get_bf_pflux, get_bf_pwave, get_bf_sed, get_bf_cont, get_bf_calib, get_full_grid

include("SFHUtils.jl")
using .SFHUtils
export zred_to_agebins, build_agebins, get_agebins, logmass_to_masses, get_sfh, get_sfr, sfh_lookback

include("PlotUtils.jl")
using .PlotUtils
export AbstractProspectPlot, AbstractSedPlot, AbstractPhotoPlot, AbstractMultiplier, AbstractErrors
export ObsSpecPlot, BFSpecPlot, BFSedPlot, DensityPlot, ObsPhotoPlot, BFPhotoPlot, BinnedSFH, LineIndicators
export MaskedSpecPlot, MaskedPhotoPlot
export SedErrors, PhotoErrors, SFHErrors, Errors, Mask, Calibration, PlotCommand
export get_zorder, get_plot_function, get_error_plot_function, get_plot_data, render!

end # module ProspectorUtils
