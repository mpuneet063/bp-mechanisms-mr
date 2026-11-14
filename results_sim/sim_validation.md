# Simulation validation

Planted mechanisms: {'bp_only': 250, 'vascular': 119, 'adiposity': 90, 'cardiac': 90, 'kidney': 79, 'lipid': 69}

| partition | ARI vs truth | K found |
|---|---|---|
| hclust | 0.570 | 5 |
| bnmf | 0.397 | 5 |

## hclust: MR estimates vs planted effects (per 10 mmHg higher SBP)

- correlation(estimate, planted) = 0.939; 95% CI covers the planted value in 91% of 35 cluster x outcome estimates (expect ~95%, lower for the adiposity x CAD cells with planted pleiotropy).
- c0: majority 'cardiac' (purity 96%)
- c1: majority 'vascular' (purity 97%)
- c2: majority 'adiposity' (purity 100%)
- c3: majority 'kidney' (purity 98%)
- c4: majority 'bp_only' (purity 66%)

## bnmf: MR estimates vs planted effects (per 10 mmHg higher SBP)

- correlation(estimate, planted) = 0.743; 95% CI covers the planted value in 94% of 35 cluster x outcome estimates (expect ~95%, lower for the adiposity x CAD cells with planted pleiotropy).
- c0: majority 'cardiac' (purity 58%)
- c1: majority 'lipid' (purity 100%)
- c2: majority 'adiposity' (purity 100%)
- c3: majority 'kidney' (purity 98%)
- c4: majority 'bp_only' (purity 55%)

## MVMR check (planted BMI->CAD pleiotropy in the adiposity cluster)

- c2: CAD OR univariable 1.350 -> MVMR 1.286 (planted SBP effect OR 1.284)
