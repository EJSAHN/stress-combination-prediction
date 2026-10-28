PLANT STRESS COMBINATIONS — CODE6 EXECUTED RECORD

OUTCOME: CONTINUE, WITH A NARROWED CLAIM
The fixed primary test (WT, mean FPKM >=1) met all four descriptive sign criteria. The equal mean of all three pair-stress profiles predicted the triple-stress profile, and its component outside the fitted single-stress span contributed positive MSE gain. Both threshold sensitivities and the prespecified mutant secondary analysis also met the same signs. This is not a significance test or publication-readiness decision.

PRIMARY RESULTS: WT, 19,180 genes, three target RNA libraries
Model | Mean MSE | Reduction from control (%)
control | 0.701011984 | 0.000
single_additive | 0.661835168 | 5.589
second_order_factorial | 0.443886461 | 36.679
within_only | 0.179714739 | 74.364
outside_only | 0.638738965 | 8.883
full_pair_centroid | 0.117441720 | 83.247

Adding the fixed outside component to the within-only prediction reduces MSE by 34.651%.
The outside component contributes 0.062273019 to total MSE gain, or 10.671% of full-centroid gain. These percentages are prediction-error comparisons, not biological variance explained.
Outside gains by target library: 0.053418980, 0.061648220, 0.071751857
Corrected outside proximity margins by control split: 0.241877385, 0.231165505, 0.208420151
The latter use different reference bases and a single-stress-library-only gene filter. They are population-geometry diagnostics, not directly comparable magnitudes of observed prediction gain. The three splits overlap.

WHAT WAS FIXED BEFORE OPENING
Primary WT endpoint; secondary RAP2.3 mutant; filters 1 (primary), 0.5 and 2; log2(1+FPKM); full numerical rank of three single-condition mean displacements; equal weights for all three pair-condition donors; six models; MSE primary and MAE secondary; disjoint-control covariance correction; four positive-sign requirements at numerical tolerance 1e-12.
All genotype/filter scenarios were frozen before any triple-condition numerical columns were parsed. Fitting requested only 21 nontarget columns per genotype. Source byte hashes and header/gene-ID inspection preceded fitting. A poison-target test replaced all target fields with nonnumeric strings and reproduced every frozen array exactly. Independent QR reconstruction checked prediction vectors and bases.
The source has 32,826 rows per matrix. The first freeze correctly stopped at duplicated IDs. INPUT_AMENDMENT_01.json, fixed before target opening, excludes all four rows for AT1G79790 and AT2G04725. No averaging or representative row selection was used.
A truncated NPZ was detected before opening. That freeze was rejected and all six scenarios regenerated with atomic persistence, full-array round trips, CRC and hash checks. The cause is unknown; the scientific calculations were unchanged. Only frozen_verified/ is the valid production freeze.

SCIENTIFIC SCOPE
This is GSE330879, a separate Arabidopsis experiment with the same research group as the discovery source: WT and rap2.3_1, each with the complete salt/excess-light/heat three-factor cube and three pooled-seedling RNA libraries per cell. Replicate labels were not used to infer paired blocks. Author FPKM (TAIR10, HISAT2/Cufflinks) is not raw count data. Absolute MSE should not be compared across expression-unit scales.
It tests a narrowed pair-to-triple hypothesis. It is not exact replication of the earlier selected-triple-to-fourway result, transfer of fitted coefficients across species, independent-laboratory validation or crop validation.
The outside component is outside an ESTIMATED observed single-stress span. Noisy means may omit real single-response directions. Interactions may also lie inside that span. The second-order predictor has exactly three times the pair-centroid outside component; a checked algebraic identity. Thus transfer does not establish pure higher-order interaction, a new superior algorithm, pathway activation, causal protection or stress memory.
Author preprocessing can couple libraries, and qualitative published context was known. There are only three biological RNA pools per condition; gene number does not provide extra biological replication. No gene-level significance test or gene bootstrap was used.

NEXT SCIENTIFIC QUESTION
Quantify whether uncertainty in the fitted single-stress span can explain the outside component. Preserve the fixed external result; do not retune against its now-open target. Exact four-factor replication and crop transfer remain separate tasks.
Tomato GSE294808 was excluded before expression acquisition because controls/processing were unresolved. Citrus GSE203331 was structurally eligible but deferred. A public repeat-layer reanalysis by Zhijun Wu exists; that is not evidence that this user has previously analysed it. Neither supplied numerical evidence to this result.

SOURCES AND FILES
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE330879
Sinha et al., Plant, Cell & Environment, online 8 October 2026. DOI 10.1111/pce.70963. https://pubmed.ncbi.nlm.nih.gov/42847751/
The workbook is a readable view of all result tables. The archive includes original FPKM matrices, executed English code, fixed plan plus amendment, frozen arrays, results, source audits and independent verification. MANIFEST_SHA256.json lists hashes. No new environment or installation was used; SciNet was unnecessary.
This is an executed analysis record. Final GitHub/Windows packaging and the requested local confirmation are deferred until the manuscript pipeline is settled.
