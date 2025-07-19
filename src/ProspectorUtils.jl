module ProspectorUtils

include("DataUtils.jl")
using .DataUtils
export ProspectResults, ProspectorBestFit, ProspectorObs, ProspectorObs
export get_z, get_phot_flux, get_phot_wave, get_spec_cont, get_spec_flux, get_spec_err, get_spec_wave, get_spec_calib, labels, bestfit, maggies2μJy
export get_obs_sflux, get_obs_swave, get_obs_serr, get_obs_pflux, get_obs_pwave, get_obs_perr


include("SFHUtils.jl")
using .SFHUtils
export zred_to_agebins, build_agebins, logmass_to_masses, get_sfh, pyshow

include("PlotUtils.jl")
using .PlotUtils
export AbstractProspectPlot, AbstractSedPlot, AbstractPhotoPlot, AbstractMultiplier
export ObsSpecPlot, BFSedPlot, DensityPlot, ObsPhotoPlot, BFPhotoPlot
export Errors, MaskPlot, Calibration
export PlotCommand
export get_plot_function, get_error_plot_function
export get_plot_data
export render!


end # module ProspectorUtils
