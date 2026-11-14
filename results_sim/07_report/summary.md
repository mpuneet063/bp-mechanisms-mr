# Results summary (SIMULATED DATA)

Exposure scaling: per 10 mmHg higher SBP. Clustering variants: 697; MR instruments: 614. Traits: sbp, pp, dbp, ldl, bmi, hr, egfr, urate, t2d, tg, hdl, crp, fi, uacr, whradjbmi, alt, na_urine, k_urine.

## Clusters (K = 5, bNMF K counts {'4': 1, '5': 91, '6': 8}; hclust vs bNMF ARI = 0.34)

| cluster | n | top traits | bootstrap Jaccard | noise keep prob | stability | to MR |
|---|---|---|---|---|---|---|
| c0 | 84 | hr+, dbp+, pp- | 0.82 | 0.81 | stable | yes |
| c1 | 111 | pp+, hr-, dbp- | 0.86 | 0.82 | stable | yes |
| c2 | 72 | bmi+, fi+, t2d+ | 0.96 | 0.93 | stable | yes |
| c3 | 61 | egfr-, urate+, uacr+ | 0.69 | 0.85 | moderate | yes |
| c4 | 369 | pp-, bmi-, dbp- | 0.84 | 0.92 | stable | yes |

## IVW odds ratios per 10 mmHg higher SBP (hclust)

| cluster | Any stroke | Ischaemic stroke | Large-artery stroke | Cardioembolic stroke | Small-vessel stroke | Intracerebral haemorrhage | Coronary artery disease | FinnGen stroke | FinnGen ischaemic stroke | FinnGen CHD |
|---|---|---|---|---|---|---|---|---|---|---|
| all | 1.34 (1.33-1.36) | 1.33 (1.31-1.35) | 1.36 (1.31-1.41) | 1.23 (1.19-1.26) | 1.42 (1.37-1.47) | 1.28 (1.15-1.42) | 1.35 (1.33-1.36) | 1.34 (1.30-1.38) | 1.33 (1.29-1.38) | 1.38 (1.35-1.41) |
| c0 | 1.35 (1.31-1.40) | 1.34 (1.29-1.38) | 1.35 (1.21-1.50) | 1.62 (1.50-1.76) | 1.24 (1.12-1.37) | 1.35 (1.01-1.82) | 1.23 (1.21-1.26) | 1.40 (1.29-1.52) | 1.42 (1.30-1.56) | 1.22 (1.15-1.29) |
| c1 | 1.49 (1.45-1.54) | 1.43 (1.38-1.48) | 1.32 (1.20-1.46) | 1.05 (0.97-1.13) | 2.10 (1.91-2.31) | 2.15 (1.59-2.89) | 1.25 (1.23-1.28) | 1.52 (1.41-1.64) | 1.45 (1.31-1.60) | 1.33 (1.26-1.41) |
| c2 | 1.24 (1.19-1.28) | 1.23 (1.18-1.29) | 1.45 (1.29-1.62) | 1.33 (1.19-1.48) | 1.29 (1.15-1.44) | 1.02 (0.73-1.43) | 1.35 (1.32-1.38) | 1.09 (0.99-1.20) | 1.18 (1.06-1.31) | 1.42 (1.32-1.52) |
| c3 | 1.29 (1.24-1.35) | 1.29 (1.24-1.35) | 1.40 (1.24-1.58) | 1.10 (1.00-1.21) | 1.53 (1.34-1.74) | 1.12 (0.78-1.62) | 1.57 (1.52-1.61) | 1.20 (1.09-1.33) | 1.28 (1.14-1.43) | 1.69 (1.57-1.82) |
| c4 | 1.33 (1.31-1.35) | 1.33 (1.31-1.35) | 1.35 (1.29-1.42) | 1.20 (1.16-1.25) | 1.33 (1.27-1.40) | 1.19 (1.03-1.36) | 1.37 (1.35-1.39) | 1.36 (1.30-1.41) | 1.32 (1.26-1.39) | 1.38 (1.35-1.42) |

## Do clusters differ? (stable clusters only)

| outcome | Q_between | df | p (chi2) | p (permutation) |
|---|---|---|---|---|
| Any stroke | 70.2 | 4 | 2.1e-14 | 0.0005 |
| Ischaemic stroke | 31.9 | 4 | 2e-06 | 0.0005 |
| Large-artery stroke | 1.8 | 4 | 0.77 | 0.728 |
| Cardioembolic stroke | 71.4 | 4 | 1.2e-14 | 0.0005 |
| Small-vessel stroke | 83.0 | 4 | 4e-17 | 0.0005 |
| Intracerebral haemorrhage | 15.1 | 4 | 0.0045 | 0.002 |
| Coronary artery disease | 236.6 | 4 | 5e-50 | 0.0005 |
| FinnGen stroke | 34.4 | 4 | 6.2e-07 | 0.0005 |
| FinnGen ischaemic stroke | 10.6 | 4 | 0.031 | 0.028 |
| FinnGen CHD | 46.4 | 4 | 2.1e-09 | 0.0005 |

## Strongest pairwise contrasts

- Coronary artery disease: c3 (egfr-, urate+, uacr+) vs c0 (hr+, dbp+, pp-): log-OR ratio 2.14 (bootstrap 95% CI 1.89 to 2.46); OR ratio 1.27 per 10 mmHg higher SBP; difference p = 4.6e-40
- Coronary artery disease: c3 (egfr-, urate+, uacr+) vs c1 (pp+, hr-, dbp-): log-OR ratio 1.99 (bootstrap 95% CI 1.79 to 2.22); OR ratio 1.25 per 10 mmHg higher SBP; difference p = 3.7e-39
- Coronary artery disease: c3 (egfr-, urate+, uacr+) vs c4 (pp-, bmi-, dbp-): log-OR ratio 1.43 (bootstrap 95% CI 1.33 to 1.53); OR ratio 1.14 per 10 mmHg higher SBP; difference p = 2.8e-19
- Small-vessel stroke: c1 (pp+, hr-, dbp-) vs c4 (pp-, bmi-, dbp-): log-OR ratio 2.57 (bootstrap 95% CI 2.12 to 3.16); OR ratio 1.57 per 10 mmHg higher SBP; difference p = 5.3e-17
- Coronary artery disease: c3 (egfr-, urate+, uacr+) vs c2 (bmi+, fi+, t2d+): log-OR ratio 1.50 (bootstrap 95% CI 1.35 to 1.66); OR ratio 1.16 per 10 mmHg higher SBP; difference p = 1.3e-15
- Coronary artery disease: c4 (pp-, bmi-, dbp-) vs c0 (hr+, dbp+, pp-): log-OR ratio 1.50 (bootstrap 95% CI 1.34 to 1.72); OR ratio 1.11 per 10 mmHg higher SBP; difference p = 1.4e-15
- Any stroke: c1 (pp+, hr-, dbp-) vs c2 (bmi+, fi+, t2d+): log-OR ratio 1.89 (bootstrap 95% CI 1.61 to 2.27); OR ratio 1.21 per 10 mmHg higher SBP; difference p = 3.4e-15
- Coronary artery disease: c4 (pp-, bmi-, dbp-) vs c1 (pp+, hr-, dbp-): log-OR ratio 1.40 (bootstrap 95% CI 1.27 to 1.55); OR ratio 1.09 per 10 mmHg higher SBP; difference p = 3.1e-14

## Pleiotropy checks

- MR-Egger intercept p < 0.05 in 4 of 60 cluster x outcome analyses.
- MR-PRESSO global test p < 0.05 in 8 analyses.
