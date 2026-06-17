# SFH branch comparison

- input: `/home/matteo/Dottorato/DATA/DAWN_Archive/prospector/results/abell2744-castellano1-v4_prism-clear_3073_18819/abell2744-castellano1-v4_prism-clear_3073_18819.h5`
- chain rows: `11277`
- chain parameters: `22`
- has weights: `true`
- default current estimator: `WeightedMedian`
- redshift: `4.6513`
- stored bestfit agebins available but ignored: `true`
- max abs agebin diff current-vs-main zred: `0.0`

| Method | logmass input | max abs SFH diff vs current default | max abs SFH diff vs main canonical | total formed mass | recent log10 SFR |
|---|---:|---:|---:|---:|---:|
| current default: WeightedMedian + zred agebins | 9.55290200 | 0.000000e+00 | 8.049401e-01 | 3.571922e+09 | 1.208103 |
| current Median + zred agebins | 9.61528665 | 8.049401e-01 | 0.000000e+00 | 4.123696e+09 | 1.077913 |
| current BestFit + zred agebins | 9.65980354 | 7.366265e-01 | 3.124261e-01 | 4.568815e+09 | 1.293474 |
| main canonical: unweighted median + zred agebins | 9.61528665 | 8.049401e-01 | 0.000000e+00 | 4.123696e+09 | 1.077913 |

## Ratio Inputs

- current default: WeightedMedian + zred agebins: `[1.490389e-01, 2.293166e-02, 4.156627e-01, -6.866298e-02, 9.856407e-01, 6.092994e-02, 4.902812e-01]`
- current Median + zred agebins: `[1.135038e-01, 2.844013e-02, 4.090131e-01, -4.650584e-02, 2.361547e-01, 3.600910e-02, 3.440773e-01]`
- current BestFit + zred agebins: `[4.693549e-02, -1.856747e-03, 9.064844e-01, 2.386033e-02, 1.849533e-01, 1.702681e-02, 2.271631e-01]`
- main canonical: unweighted median + zred agebins: `[1.135038e-01, 2.844013e-02, 4.090131e-01, -4.650584e-02, 2.361547e-01, 3.600910e-02, 3.440773e-01]`
