# Results summary (real data)

Exposure scaling: per 10 mmHg higher SBP (unit blank in OpenGWAS; inferred mmHg from median |beta| = 0.253). Clustering variants: 1490; MR instruments: 363. Traits: sbp, ldl, tg, hdl, height, fg, egfr, hr, alcohol, crp, hba1c, pp, dbp, t2d, bmi, urate, alt, bodyfat, fi, smoking, na_urine, k_urine, uacr.

## Clusters (kmeans, K = 6: largest K with every cluster bootstrap Jaccard >= 0.6 (gate); weakest cluster 0.62, mean 0.68)

bNMF K counts {'9': 11, '10': 26, '11': 41, '12': 20, '13': 2}. Agreement between methods (ARI): kmeans vs hclust 0.21, kmeans vs bnmf 0.27, hclust vs bnmf 0.12. 4 of 6 clusters pass the stability gate (Jaccard >= 0.6).

| cluster | n | top traits | bootstrap Jaccard | noise keep prob | stability | to MR |
|---|---|---|---|---|---|---|
| c0 | 251 | bmi-, bodyfat-, hdl+ | 0.57 | 0.67 | unstable | no |
| c1 | 338 | bodyfat+, bmi+, hdl- | 0.66 | 0.69 | moderate | yes |
| c2 | 237 | pp+, dbp-, height+ | 0.64 | 0.67 | moderate | yes |
| c3 | 181 | height-, pp+, ldl+ | 0.63 | 0.71 | moderate | yes |
| c4 | 380 | height+, tg-, fi- | 0.59 | 0.65 | unstable | no |
| c5 | 103 | tg+, hdl-, fi+ | 0.77 | 0.78 | stable | yes |

## IVW odds ratios per 10 mmHg higher SBP (unit blank in OpenGWAS; inferred mmHg from median |beta| = 0.253) (kmeans)

| cluster | Any stroke | Ischaemic stroke | Large-artery stroke | Cardioembolic stroke | Small-vessel stroke | Intracerebral haemorrhage | Coronary artery disease | FinnGen stroke | FinnGen ischaemic stroke | FinnGen CHD |
|---|---|---|---|---|---|---|---|---|---|---|
| all | 1.36 (1.30-1.42) | 1.36 (1.30-1.43) | 1.68 (1.50-1.87) | 1.18 (1.08-1.29) | 1.50 (1.35-1.68) | 1.31 (1.16-1.47) | 1.36 (1.28-1.45) | 1.29 (1.21-1.38) | 1.30 (1.21-1.39) | 1.33 (1.23-1.44) |
| c0 | 1.03 (0.92-1.14) | 1.00 (0.87-1.15) | 1.18 (0.84-1.64) | 0.94 (0.75-1.18) | 1.17 (0.86-1.58) | 1.06 (0.78-1.43) | 1.30 (1.13-1.51) | 1.09 (0.90-1.32) | 1.08 (0.89-1.32) | 1.33 (1.10-1.60) |
| c1 | 1.33 (1.21-1.48) | 1.40 (1.24-1.58) | 1.84 (1.44-2.36) | 1.13 (0.92-1.39) | 1.54 (1.19-1.99) | 1.16 (0.90-1.50) | 1.46 (1.27-1.69) | 1.09 (0.94-1.25) | 1.12 (0.97-1.30) | 1.28 (1.08-1.53) |
| c2 | 1.48 (1.33-1.65) | 1.54 (1.35-1.75) | 2.41 (1.83-3.18) | 1.28 (1.00-1.65) | 1.72 (1.31-2.26) | 1.21 (0.87-1.68) | 1.30 (1.00-1.68) | 1.38 (1.18-1.63) | 1.34 (1.12-1.60) | 1.36 (0.99-1.86) |
| c3 | 1.51 (1.30-1.76) | 1.41 (1.20-1.67) | 1.36 (0.96-1.95) | 1.22 (0.91-1.63) | 1.50 (1.04-2.16) | 1.20 (0.84-1.71) | 1.47 (1.25-1.74) | 1.52 (1.22-1.89) | 1.53 (1.21-1.94) | 1.45 (1.15-1.84) |
| c4 | 1.44 (1.35-1.54) | 1.42 (1.31-1.53) | 1.63 (1.37-1.93) | 1.33 (1.16-1.54) | 1.58 (1.32-1.89) | 1.58 (1.29-1.92) | 1.32 (1.20-1.45) | 1.37 (1.24-1.51) | 1.37 (1.24-1.52) | 1.32 (1.19-1.47) |
| c5 | 1.31 (1.10-1.55) | 1.28 (1.01-1.61) | 2.14 (1.28-3.55) | 0.73 (0.51-1.05) | 1.20 (0.76-1.90) | 1.16 (0.68-1.97) | 1.47 (1.19-1.82) | 1.43 (1.07-1.90) | 1.55 (1.17-2.06) | 1.38 (1.06-1.79) |

## Do clusters differ? (stable clusters only)

| outcome | Q_between | df | p (chi2) | p (permutation) |
|---|---|---|---|---|
| Any stroke | 3.4 | 3 | 0.33 | 0.371 |
| Ischaemic stroke | 2.2 | 3 | 0.53 | 0.554 |
| Large-artery stroke | 6.4 | 3 | 0.093 | 0.0865 |
| Cardioembolic stroke | 6.9 | 3 | 0.076 | 0.139 |
| Small-vessel stroke | 1.8 | 3 | 0.62 | 0.618 |
| Intracerebral haemorrhage | 0.0 | 3 | 1 | 0.999 |
| Coronary artery disease | 0.8 | 3 | 0.86 | 0.902 |
| FinnGen stroke | 9.0 | 3 | 0.03 | 0.065 |
| FinnGen ischaemic stroke | 7.2 | 3 | 0.065 | 0.112 |
| FinnGen CHD | 0.7 | 3 | 0.87 | 0.922 |

## Strongest pairwise contrasts

- Large-artery stroke: c2 (pp+, dbp-, height+) vs c1 (bodyfat+, bmi+, hdl-): log-OR ratio 1.44 (bootstrap 95% CI 0.88 to 2.50); OR ratio 1.31 per 10 mmHg higher SBP (unit blank in OpenGWAS; inferred mmHg from median |beta| = 0.253); difference p = 0.16
- Ischaemic stroke: c2 (pp+, dbp-, height+) vs c5 (tg+, hdl-, fi+): log-OR ratio 1.76 (bootstrap 95% CI 0.74 to 10.50); OR ratio 1.20 per 10 mmHg higher SBP (unit blank in OpenGWAS; inferred mmHg from median |beta| = 0.253); difference p = 0.17
- Any stroke: c2 (pp+, dbp-, height+) vs c1 (bodyfat+, bmi+, hdl-): log-OR ratio 1.36 (bootstrap 95% CI 0.88 to 2.29); OR ratio 1.11 per 10 mmHg higher SBP (unit blank in OpenGWAS; inferred mmHg from median |beta| = 0.253); difference p = 0.17
- Any stroke: c3 (height-, pp+, ldl+) vs c1 (bodyfat+, bmi+, hdl-): log-OR ratio 1.43 (bootstrap 95% CI 0.81 to 2.42); OR ratio 1.13 per 10 mmHg higher SBP (unit blank in OpenGWAS; inferred mmHg from median |beta| = 0.253); difference p = 0.18
- Any stroke: c3 (height-, pp+, ldl+) vs c5 (tg+, hdl-, fi+): log-OR ratio 1.55 (bootstrap 95% CI 0.76 to 4.87); OR ratio 1.16 per 10 mmHg higher SBP (unit blank in OpenGWAS; inferred mmHg from median |beta| = 0.253); difference p = 0.21
- Any stroke: c2 (pp+, dbp-, height+) vs c5 (tg+, hdl-, fi+): log-OR ratio 1.48 (bootstrap 95% CI 0.79 to 4.67); OR ratio 1.13 per 10 mmHg higher SBP (unit blank in OpenGWAS; inferred mmHg from median |beta| = 0.253); difference p = 0.22
- Ischaemic stroke: c2 (pp+, dbp-, height+) vs c1 (bodyfat+, bmi+, hdl-): log-OR ratio 1.27 (bootstrap 95% CI 0.77 to 2.11); OR ratio 1.10 per 10 mmHg higher SBP (unit blank in OpenGWAS; inferred mmHg from median |beta| = 0.253); difference p = 0.32
- FinnGen ischaemic stroke: c3 (height-, pp+, ldl+) vs c2 (pp+, dbp-, height+): log-OR ratio 1.47 (bootstrap 95% CI 0.57 to 5.96); OR ratio 1.15 per 10 mmHg higher SBP (unit blank in OpenGWAS; inferred mmHg from median |beta| = 0.253); difference p = 0.37

## Pleiotropy checks

- MR-Egger intercept p < 0.05 in 2 of 70 cluster x outcome analyses.
- MR-PRESSO global test p < 0.05 in 43 analyses.

## Univariable vs multivariable MR (adjusted for BMI, LDL, T2D)

| cluster | outcome | OR univariable | OR MVMR (95% CI) |
|---|---|---|---|
| c1 bodyfat+, bmi+, hdl- | Any stroke | 1.33 | 1.29 (1.18-1.42) |
| c2 pp+, dbp-, height+ | Any stroke | 1.48 | 1.47 (1.31-1.64) |
| c3 height-, pp+, ldl+ | Any stroke | 1.51 | 1.50 (1.31-1.71) |
| c5 tg+, hdl-, fi+ | Any stroke | 1.31 | 1.24 (1.03-1.49) |
| c1 bodyfat+, bmi+, hdl- | Ischaemic stroke | 1.40 | 1.35 (1.21-1.51) |
| c2 pp+, dbp-, height+ | Ischaemic stroke | 1.54 | 1.52 (1.33-1.74) |
| c3 height-, pp+, ldl+ | Ischaemic stroke | 1.41 | 1.40 (1.19-1.64) |
| c5 tg+, hdl-, fi+ | Ischaemic stroke | 1.28 | 1.21 (0.97-1.50) |
| c1 bodyfat+, bmi+, hdl- | Large-artery stroke | 1.84 | 1.71 (1.32-2.21) |
| c2 pp+, dbp-, height+ | Large-artery stroke | 2.41 | 2.35 (1.73-3.20) |
| c3 height-, pp+, ldl+ | Large-artery stroke | 1.36 | 1.33 (0.92-1.93) |
| c5 tg+, hdl-, fi+ | Large-artery stroke | 2.14 | 1.74 (1.06-2.87) |
| c1 bodyfat+, bmi+, hdl- | Cardioembolic stroke | 1.13 | 1.09 (0.90-1.33) |
| c2 pp+, dbp-, height+ | Cardioembolic stroke | 1.28 | 1.27 (1.00-1.60) |
| c3 height-, pp+, ldl+ | Cardioembolic stroke | 1.22 | 1.21 (0.91-1.61) |
| c5 tg+, hdl-, fi+ | Cardioembolic stroke | 0.73 | 0.70 (0.47-1.03) |
| c1 bodyfat+, bmi+, hdl- | Small-vessel stroke | 1.54 | 1.45 (1.14-1.83) |
| c2 pp+, dbp-, height+ | Small-vessel stroke | 1.72 | 1.68 (1.27-2.23) |
| c3 height-, pp+, ldl+ | Small-vessel stroke | 1.50 | 1.49 (1.06-2.09) |
| c5 tg+, hdl-, fi+ | Small-vessel stroke | 1.20 | 1.07 (0.68-1.68) |
| c1 bodyfat+, bmi+, hdl- | Intracerebral haemorrhage | 1.16 | 1.12 (0.86-1.45) |
| c2 pp+, dbp-, height+ | Intracerebral haemorrhage | 1.21 | 1.23 (0.90-1.68) |
| c3 height-, pp+, ldl+ | Intracerebral haemorrhage | 1.20 | 1.18 (0.81-1.73) |
| c5 tg+, hdl-, fi+ | Intracerebral haemorrhage | 1.16 | 1.28 (0.77-2.11) |
| c1 bodyfat+, bmi+, hdl- | Coronary artery disease | 1.46 | 1.40 (1.22-1.60) |
| c2 pp+, dbp-, height+ | Coronary artery disease | 1.30 | 1.32 (1.12-1.55) |
| c3 height-, pp+, ldl+ | Coronary artery disease | 1.47 | 1.44 (1.19-1.75) |
| c5 tg+, hdl-, fi+ | Coronary artery disease | 1.47 | 1.34 (1.03-1.74) |
| c1 bodyfat+, bmi+, hdl- | FinnGen stroke | 1.09 | 1.06 (0.92-1.23) |
| c2 pp+, dbp-, height+ | FinnGen stroke | 1.38 | 1.37 (1.16-1.60) |
| c3 height-, pp+, ldl+ | FinnGen stroke | 1.52 | 1.50 (1.23-1.84) |
| c5 tg+, hdl-, fi+ | FinnGen stroke | 1.43 | 1.32 (1.01-1.73) |
| c1 bodyfat+, bmi+, hdl- | FinnGen ischaemic stroke | 1.12 | 1.10 (0.95-1.28) |
| c2 pp+, dbp-, height+ | FinnGen ischaemic stroke | 1.34 | 1.32 (1.11-1.56) |
| c3 height-, pp+, ldl+ | FinnGen ischaemic stroke | 1.53 | 1.52 (1.22-1.88) |
| c5 tg+, hdl-, fi+ | FinnGen ischaemic stroke | 1.55 | 1.42 (1.07-1.90) |
| c1 bodyfat+, bmi+, hdl- | FinnGen CHD | 1.28 | 1.24 (1.04-1.49) |
| c2 pp+, dbp-, height+ | FinnGen CHD | 1.36 | 1.37 (1.12-1.68) |
| c3 height-, pp+, ldl+ | FinnGen CHD | 1.45 | 1.42 (1.10-1.84) |
| c5 tg+, hdl-, fi+ | FinnGen CHD | 1.38 | 1.25 (0.89-1.76) |
