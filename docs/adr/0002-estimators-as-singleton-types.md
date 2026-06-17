# Estimators as singleton types

Point estimates from the posterior are selected with singleton types — `BestFit`,
`Median`, `WeightedMedian` (subtypes of `Estimator`) — passed to `get_mass`, `get_sfh`,
`logmass_to_masses` and resolved by dispatch, rather than a `Symbol` keyword.

This makes the choice extensible (a new estimator is a new type + `estimate` method) and
keeps `estimate`/`get_sfh` free of runtime branching. The default is `WeightedMedian` when
the chain carries dynesty importance `weights`, otherwise `Median`: dynesty produces
weighted samples, so an unweighted median would misrepresent the posterior.
