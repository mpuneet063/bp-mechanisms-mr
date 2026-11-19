# Which blood-pressure mechanisms drive stroke, and which drive coronary artery disease?

Blood-pressure (BP) variants do not all raise BP the same way. Some act alongside adiposity, some through the
kidney, some through arterial stiffness. This project groups genome-wide BP variants into **mechanism
clusters** by their effects on other traits, and then uses **Mendelian randomisation (MR)** to ask whether
raising SBP through each mechanism has the same causal effect on **stroke subtypes** and on **coronary artery
disease (CAD)**.

It builds on:

* **Keaton et al. 2024, *Nat Genet*** (2,103 BP signals, n = 1,028,980) for the variants
* **Pascat et al. 2026, *Nat Commun*** for the clustering approach (Ward hierarchical clustering, with bNMF as a cross-check)

It adds three things that paper does not report:

1. **Cluster stability.** Variant bootstrap (Jaccard), a bootstrap over GWAS sampling noise, leave-one-trait-out and a K scan. Clusters below a pre-specified threshold are flagged before any MR is run.
2. **Causal estimates per cluster.** MR on ischaemic stroke subtypes, intracerebral haemorrhage and CAD, with FinnGen replication. Stroke outcomes come from GIGASTROKE (Mishra 2022) where OpenGWAS has it, otherwise MEGASTROKE (Malik 2018). CAD comes from a consortium meta-analysis (Aragam 2022, else van der Harst 2017, else Nikpay 2015).
3. **A formal between-cluster test.** Cochran's Q, a label-permutation null and pairwise contrasts, plus multivariable MR that separates "SBP through mechanism k" from each cluster's BMI/LDL/T2D pleiotropy.

Everything uses public summary statistics. No individual-level data is needed.

## Run it on Kaggle

1. Upload this folder as a Kaggle Dataset. Optionally add the Keaton et al. 2024 Supplementary Tables `.xlsx` to the same dataset.
2. In a notebook, turn on **Settings → Internet**. Add a secret named `OPENGWAS_JWT` under **Add-ons → Secrets**. You can get a free token at <https://api.opengwas.io/profile/>.
3. Run one cell:

```python
import os
from kaggle_secrets import UserSecretsClient
os.environ["OPENGWAS_JWT"] = UserSecretsClient().get_secret("OPENGWAS_JWT")
os.environ["KEATON_XLSX"] = "/kaggle/input/<dataset>/keaton_supp.xlsx"   # optional
!bash /kaggle/input/<dataset>/bp-mechanisms-mr/kaggle_run.sh
```

The script copies the repo to `/kaggle/working` and installs dependencies. It runs the unit tests and the
simulation check, fetches the data, runs every experiment and writes
`/kaggle/working/results_bundle.zip`. Use `SKIP_SIM=1` to skip the simulation check, or `FAST=1` for a quick
smoke test.

## Run it on Google Colab (progress survives runtime restarts)

Keep the repo on Google Drive. The fetch caches every API response in `data/raw/cache/`, and each finished
analysis step leaves a `.done` marker. After a runtime dies, re-run cells 1 and 3 and the run continues where
it stopped. Upload `bp-mechanisms-mr.zip` (and optionally the Keaton workbook) to your Drive first, and add
`OPENGWAS_JWT` under Colab **Secrets**.

```python
# Cell 1 - every session
from google.colab import drive, userdata
import os
drive.mount("/content/drive")
os.environ["OPENGWAS_JWT"] = userdata.get("OPENGWAS_JWT")
os.environ["KEATON_XLSX"] = "/content/drive/MyDrive/keaton_supp.xlsx"   # optional
```
```python
# Cell 2 - once: unpack the repo onto Drive
import zipfile, os
if not os.path.exists("/content/drive/MyDrive/bp-mechanisms-mr"):
    zipfile.ZipFile("/content/drive/MyDrive/bp-mechanisms-mr.zip").extractall("/content/drive/MyDrive/")
```
```bash
%%bash
# Cell 3 - every session: resumes automatically
cd /content/drive/MyDrive/bp-mechanisms-mr
SKIP_SIM=1 bash kaggle_run.sh
```
The results bundle is written to `MyDrive/results_bundle.zip`. Use a CPU runtime; the GPU is never used.

### Running the pieces yourself

| Script | What it does |
|---|---|
| `setup.sh` | install Python dependencies |
| `fetch_data.sh` | **data requirements**: resolves OpenGWAS study IDs, then streams each study's full summary-statistics file and keeps only our variants (about 1-2 h for 34 studies; one checkpoint per study). Prints the study manifest and a per-study missing-variant report. `FETCH_SOURCE=api` uses the old per-variant look-ups (~2 min per 64 variants). |
| `run_experiments.sh [--fast] [--raw DIR] [--out DIR]` | runs experiments 01–07 in order, logging to `<out>/run.log` |
| `run_simulation_check.sh [--fast]` | unit tests, then the full pipeline on simulated data with known truth |
| `kaggle_run.sh` | all of the above, in one go |

### Experiments (`scripts/`, all Python)

| Step | Output folder | Content |
|---|---|---|
| `00_fetch_opengwas.py` | `data/raw/` | study manifest, long association table, instrument list with LD-clump flags |
| `01_prepare.py` | `01_prepare/` | harmonise every association to the SBP-raising allele; variant × trait z-matrix (missingness filter, Bonferroni signal filter, \|r\| > 0.85 trait pruning, sample-size scaling) |
| `02_cluster.py` | `02_cluster/` | K chosen by bootstrap stability of k-means (K = 3–8); k-means (primary) and Ward at that K; bNMF (Udler/Smith ARD algorithm, ported and checked against R) with its own K; φ scan; silhouette; pairwise ARIs; cluster profiles and names |
| `03_stability.py` | `03_stability/` | for k-means and Ward: variant bootstrap Jaccard, noise bootstrap co-assignment, leave-one-trait-out, K ± 1, stability by K; bNMF noise bootstrap; **MR gate** |
| `04_mr.py` | `04_mr/` | IVW (random effects), MR-Egger, weighted median, MR-PRESSO and Steiger-filtered IVW, for all variants and each cluster × outcome |
| `05_compare.py` | `05_compare/` | Q_between, permutation p, pairwise log-OR differences and ratios (bootstrap CIs), stroke-vs-CAD contrast per cluster |
| `06_mvmr.py` | `06_mvmr/` | multivariable MR: one SBP exposure per cluster, plus BMI, LDL and T2D; conditional F |
| `07_report.py` | `07_report/` | figures (PNG + PDF) and `summary.md` with every headline number |

## Design decisions (and why)

* **Clustering input is each variant's per-unit-SBP profile:** trait z ÷ SBP z, with the SBP column removed. Every variant was selected for raising SBP, so raw z mostly encodes effect size, not mechanism. In simulation this roughly doubled cluster recovery for both methods (Ward ARI 0.32 → 0.66 at K = 6). `--profile raw` reproduces the Pascat-style input as a sensitivity analysis.
* **Stroke, CAD, MI, HF and AF are never clustering traits.** The February 2026 paper clustered on them. Testing clusters against traits used to build them would be circular.
* **Primary partition is k-means, with K chosen by stability.** On the real matrix (1,490 variants), Ward assignments did not survive the variant bootstrap at any K (mean Jaccard 0.3–0.5, every cluster below the gate). k-means with 20 restarts gave clusters that did (K = 5: Jaccard 0.72–0.82). The BP-variant profiles form a continuum, and Ward's greedy merges are fragile on a continuum. K is the largest value in 3–8 whose every cluster is stable (bootstrap Jaccard ≥ 0.75); failing that, the largest whose every cluster passes the 0.6 gate. In simulation this recovers the planted K = 6, where "highest mean Jaccard" merged two mechanisms (K = 4). Ward (Pascat 2026) at the same K, and bNMF with its own ARD K, are reported as sensitivity analyses. `--primary hclust` and `--k` override. **This choice was made after seeing Ward fail on the real data, and the write-up says so.**
* **The stability gate is fixed before MR** (mean bootstrap Jaccard ≥ 0.6). Clusters that fail are still estimated, but flagged in grey and excluded from the between-cluster tests. If no cluster passes, steps 05–06 run on all clusters and label every row exploratory.
* **Study coverage is checked, not assumed.** A study missing more than 20% of the BP variants (sparse imputation, e.g. HapMap-era or exome-based files) is replaced by the next-ranked candidate, up to 4 per trait; a trait with no dense study is dropped and logged in `manifest.csv`. For binary traits, case-control GWAS (log-odds) rank above biobank linear-model releases that report no case count.
* **Exposure units are inferred when OpenGWAS leaves them blank** (Evangelou 2018 does). Significant SBP effects are ~0.2–1 mmHg per allele versus ~0.01–0.05 in SD units, so the median |β| identifies the scale. Results are reported per 10 mmHg.
* **"X times more" is reported carefully.** The primary contrast is the difference in log OR (a Wald test). The ratio of log ORs comes with a bootstrap CI, and is called interpretable only when the denominator's CI excludes zero.
* **Instruments:** clustering uses BP variants independent at r² < 0.05. MR uses the subset independent at r² < 0.001 with F ≥ 10. Exposure betas come from the SBP GWAS in OpenGWAS (Evangelou 2018), and the Keaton 2024 signals extend the variant list when supplied.

## Validation on simulated data

`run_simulation_check.sh` plants 6 mechanisms (700 variants), cluster-specific causal effects and BMI→CAD
pleiotropy, writing the data in the same format as the fetch script. Fast-mode result:

* Ward clusters: 4 of 5 are 96–100% pure for one planted mechanism (overall ARI 0.57). The fifth absorbs the planted "BP-only" group, which has no profile to find.
* MR estimates correlate 0.94 with the planted effects, and the 95% CIs cover the planted value in 91% of cluster × outcome cells.
* MVMR removes the planted pleiotropy: the adiposity cluster's CAD OR goes from 1.35 (univariable) to 1.29 (MVMR), against a planted 1.28.
* The bNMF port reproduces the reference R implementation's K distribution on the same inputs.

## Limitations to state in any write-up

* **Sample overlap.** UK Biobank is in the SBP, GIGASTROKE and CAD GWAS, which biases MR towards the observational estimate. FinnGen replication is the check.
* **European ancestry only.**
* **A cluster's MR estimate is the effect of that set of variants.** It becomes "BP acting through mechanism k" only to the extent the MVMR adjustment captures the pleiotropy.
* **The Steiger filter uses a latent-scale r² for binary outcomes without prevalence adjustment.** It is a sensitivity analysis only.

## References

* Keaton JM et al. Genome-wide analysis in over 1 million individuals of European ancestry yields improved polygenic risk scores for blood pressure traits. *Nat Genet* 2024;56:778–791.
* Pascat V et al. Partitioned polygenic scores show mechanistic heterogeneity in type 2 diabetes and hypertension comorbidity. *Nat Commun* 2026.
* Udler MS et al. Type 2 diabetes genetic loci informed by multi-trait associations point to disease mechanisms and subtypes. *PLoS Med* 2018. bNMF code: github.com/gwas-partitioning/bnmf-clustering
* Sobczyk MK, Leyden GM et al. Distinct pathway-based effects of blood pressure and body mass index on cardiovascular traits. *Genome Med* 2025.
* Hennig C. Cluster-wise assessment of cluster stability. *Comput Stat Data Anal* 2007.
* Verbanck M et al. MR-PRESSO. *Nat Genet* 2018. Bowden J et al. Weighted median, *Genet Epidemiol* 2016; MR-Egger, *Int J Epidemiol* 2015.
* Sanderson E et al. Conditional F-statistics for multivariable MR. *Int J Epidemiol* 2019.
