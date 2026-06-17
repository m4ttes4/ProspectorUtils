# Plotting as a package extension

CairoMakie is a weak dependency. The plot *data model* (plot types, the `|>` composition
DSL, `get_zorder`) lives in the core `PlotUtils` submodule with no CairoMakie dependency;
all CairoMakie-bound behaviour (`get_plot_function`, `get_plot_data`, `render!`) lives in
`ext/ProspectorUtilsCairoMakieExt.jl` and loads automatically once the user runs
`using CairoMakie`.

This keeps the core package light to load and install for users who only need result
extraction, while plotting remains a one-import away. The trade-off: plotting methods are
unavailable until CairoMakie is loaded, which is the intended behaviour.
