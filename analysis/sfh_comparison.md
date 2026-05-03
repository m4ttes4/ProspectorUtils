# SFH comparison

- input file: `/Users/matteo/Dottorato/script/analisys/shared_galaxy_h5/capers-udsp1-v4_prism-clear_6368_21205.h5`
- redshift: `3.2142`
- theta labels with SFH ratios: `logsfr_ratios_1, logsfr_ratios_2, logsfr_ratios_3, logsfr_ratios_4, logsfr_ratios_5, logsfr_ratios_6, logsfr_ratios_7`
- chain rows x cols: `15752 x 22`
- source nbins (`src/SFHUtils.jl`): `8`
- utils nbins (`Utils.jl` path used by `get_sfh`): `8`
- bestfit logmass from file: `10.985260630891148`
- median logmass from chain: `11.48283447773861`

## Step-by-step comparison

| Step | src/SFHUtils.jl | Utils.jl | Max abs diff | Verdict |
| --- | --- | --- | --- | --- |
| `nbins` | `8` | `8` | `0` | same |
| `agebins` | `[(0.0, 7.4772), (7.4772, 8.0), (8.0, 8.250350867338568), (8.250350867338568, 8.500701734677136), (8.500701734677136, 8.751052602015706), (8.751052602015706, 9.001403469354274), (9.001403469354274, 9.251754336692843), (9.251754336692843, 9.297511827253517)]` | `[(0.0, 7.4772), (7.4772, 8.0), (8.0, 8.250350867338568), (8.250350867338568, 8.500701734677136), (8.500701734677136, 8.751052602015706), (8.751052602015706, 9.001403469354274), (9.001403469354274, 9.251754336692843), (9.251754336692843, 9.297511827253517)]` | `0.0` | same |
| `logmass` used in canonical SFH | `11.48283447773861` | `11.48283447773861` | `0.0` | same |
| `logsfr_ratios` | `[6.791501e-01, -1.220516e-01, -2.755667e-01, 2.930666e-02, 1.734181e-01, 1.017313e-01, 8.766153e-03]` | `[6.791501e-01, -1.220516e-01, -2.755667e-01, 2.930666e-02, 1.734181e-01, 1.017313e-01, 8.766153e-03]` | `0.0` | same |
| `dt` | `[3.000544e+07, 6.999456e+07, 7.797167e+07, 1.387675e+08, 2.469668e+08, 4.395309e+08, 7.822405e+08, 1.983864e+08]` | `[3.000544e+07, 6.999456e+07, 7.797167e+07, 1.387675e+08, 2.469668e+08, 4.395309e+08, 7.822405e+08, 1.983864e+08]` | `0.0` | same |
| `mass` | `[1.388114e+10, 6.778589e+09, 1.000146e+10, 3.357230e+10, 5.585027e+10, 6.667424e+10, 9.388097e+10, 2.333367e+10]` | `[1.388114e+10, 6.778589e+09, 1.000146e+10, 3.357230e+10, 5.585027e+10, 6.667424e+10, 9.388097e+10, 2.333367e+10]` | `0.0` | same |
| `sfr` | `[4.626208e+02, 9.684452e+01, 1.282705e+02, 2.419320e+02, 2.261448e+02, 1.516941e+02, 1.200155e+02, 1.176173e+02]` | `[4.626208e+02, 9.684452e+01, 1.282705e+02, 2.419320e+02, 2.261448e+02, 1.516941e+02, 1.200155e+02, 1.176173e+02]` | `0.0` | same |
| `lookback` | `[1.000000e-09, 3.000544e-02, 1.000000e-01, 1.779717e-01, 3.167391e-01, 5.637059e-01, 1.003237e+00, 1.785477e+00, 1.983864e+00]` | `[1.000000e-09, 3.000544e-02, 1.000000e-01, 1.779717e-01, 3.167391e-01, 5.637059e-01, 1.003237e+00, 1.785477e+00, 1.983864e+00]` | `0.0` | same |
| `sfh` | `[2.665225e+00, 2.665225e+00, 1.986075e+00, 2.108127e+00, 2.383693e+00, 2.354387e+00, 2.180969e+00, 2.079237e+00, 2.070471e+00]` | `[2.665225e+00, 2.665225e+00, 1.986075e+00, 2.108127e+00, 2.383693e+00, 2.354387e+00, 2.180969e+00, 2.079237e+00, 2.070471e+00]` | `0.0` | same |

## Notes

- The canonical SFH path in both codebases is numerically identical here because both use the same median chain logmass and median `logsfr_ratios` values.
- `Utils.jl` also contains plotting wrappers (`render!(::BinnedSFH)`, `plot_sfh`) that use the best-fit `logmass` branch instead of the median-chain branch; those wrappers are not equivalent to the canonical `get_sfh` path.
- The equality above is about the core SFH pipeline; any mismatch in a different file path or wrapper usually comes from input selection, not from the mass/SFR formula itself.
