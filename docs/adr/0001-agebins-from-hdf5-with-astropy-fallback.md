# Age bins from HDF5 with astropy fallback

Prospector stores the age bins it used in `bestfit/agebins`. `get_agebins` reads them
directly when present, and only reconstructs them from redshift (via astropy Planck18,
`zred_to_agebins`) when they are absent.

We keep the astropy/PythonCall dependency for the fallback rather than reimplementing the
cosmology in pure Julia: the stored bins are always preferred, and on the rare fallback
path matching Prospector's own astropy Planck18 age exactly avoids shifting the bins
relative to the user's actual fits. Numerical fidelity outweighs dropping the dependency.
