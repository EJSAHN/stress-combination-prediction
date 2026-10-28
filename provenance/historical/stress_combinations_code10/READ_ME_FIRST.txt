PLANT STRESS COMBINATIONS — CODE10 GAIN-PROFILE REPEATABILITY

RESULT
The pair-centroid rule retains its aggregate error advantage in both citrus
genotypes under every sampled complementary-pool refit. However, the gene-wise
gains agree much better within genotype than between genotypes. A recurring
aggregate advantage therefore does not imply a portable gene-level gain profile.
This is a post-result diagnostic of the existing study, not new validation.

There are 1,024 side-specific fits (256 partitions x 2 genotypes x 2 sides) and
3,072 comparator-specific mean-gain evaluations on 15,466 common features.
All 3,072 gains are positive; the smallest is 0.040334 squared log2(1+FPKM)
units per feature. No fit failed. The original CODE8/9 inputs, predictions,
protocols and results remain unchanged.

WHAT WAS MEASURED
Gain = comparator squared error minus pair-centroid squared error. It describes
where a prediction rule gains or loses accuracy relative to a specific baseline;
it is not a gene effect, differential expression or a causal contribution.
The source is the author-processed GSE203331 FPKM matrix, transformed by
log2(1+FPKM): Carrizo and Cleopatra, eight treatment conditions, and three RNA
pools per genotype/condition. Pools, not genes, are the repeat-measurement units.
Each pool contains multiple plants; those plants are not extra RNA replicates.

A. TARGET-POOL REPEATABILITY, CONDITIONAL ON THE ORIGINAL FITS
The three target-pool gain vectors were centered over the same 15,466-feature
universe. Their mean's observed variance was separated into an estimated
repeat-error contribution and an off-diagonal, cross-pool stable moment.
The conditional stable-moment fraction ranges from 0.979240 to 0.995662 across
the six genotype/comparator combinations. Cross-genotype Pearson adjustment
barely changes the original CODE9 associations:

Comparator                   Raw Pearson   Adjusted Pearson   Raw Spearman
second_order_factorial          0.118982           0.120981        0.063295
pair_trained_factor_OLS          0.153887           0.155260        0.053900
within_V3                       0.268145           0.269966        0.118451

Under the assumed common conditional mean and uncorrelated repeat errors,
weak between-genotype association is not accounted for by the target-pool
variability estimated here. This does not mean that 97.9–99.6% is biological
signal or that all noise is excluded. Shared predictor errors, correlated pool
errors and preprocessing effects can remain in the stable component. The
adjusted correlation is a plug-in ratio, not an unbiased truth estimate.
No Spearman correction, P values or confidence intervals were calculated.

B. REFITTING BOTH COMPLEMENTARY POOL SIDES
For each genotype, 256 unique partitions were fixed before the new calculations.
In each of the eight conditions, side 1 uses one RNA pool and side 2 the other
two. All lower-order predictors are fitted separately on each side, retaining
the original full genotype panels: 17,470 Carrizo and 16,247 Cleopatra features.
Target pools are excluded from fitting. Scores and vector associations use the
same 15,466 common features. Absolute MSE averages error against the respective
one or two target pools, rather than only against their mean.

For each draw, within_mean averages Car1–Car2 and Cleo1–Cleo2 correlations;
cross_mean averages Car1–Cleo2 and Car2–Cleo1, matching the one-versus-two pool
sizes. The table reports medians of these per-draw means, then the separately
computed median and quantiles of within_mean minus cross_mean:

                              Spearman means       Matched difference
Comparator                   Within    Between    Median    q0.025–q0.975
second_order_factorial        0.565      0.052      0.515      0.341–0.605
pair_trained_factor_OLS        0.584      0.050      0.534      0.336–0.639
within_V3                     0.714      0.119      0.594      0.471–0.651

All 256 matched contrasts are positive for each comparator and for Pearson,
Spearman and exact-top-10%-set Jaccard. These statements concern the specified
means, not every individual within-pair versus cross-pair comparison. Medians
of matched differences need not equal differences of marginal medians.
The Jaccard sets contain exactly 1,547 features, with lexical-ID tie breaking
and no restriction to positive gains. Complete pairings and metrics are saved.

The two sides within a genotype/draw are disjoint, but different partitions
reuse the same observed pools. The 1,024 fits are NOT mutually disjoint or
independent. The quantile ranges describe overlapping-partition sensitivity,
NOT confidence intervals; all-positive signs provide no calibrated P value.
Part B shows that the contrast is not solely a consequence of conditioning on
one particular lower-order fit. Its moderate within-genotype rank concordance
also shows that individual gene rankings are not perfectly stable. Part A's
near-one conditional fractions are not reliability estimates for all refitting.

INTERPRETATION LIMITS
The same prediction rule was fitted separately using each genotype's own data.
No model trained only in one genotype was transferred to predict the other.
The full fit panels differ, despite common scoring; differences in global OLS
or span estimation can therefore influence gain profiles. The fixed panels
were selected using both sides' source observations, and author preprocessing
is shared. Only three pools per cell support these diagnostics. They cannot
isolate causal genotype programs, remove reference-mapping or normalization
effects, or establish generalization across species, laboratories or crops.
Cross-genotype associations are low, not absent. CODE7's uncertainty about an
unobserved biological span remains unresolved. Established moment and partition
methods are used here; these results alone do not establish a new estimator,
a new mechanism or publication-level novelty. All parent outcomes were already
known when CODE10 was planned, even though its protocol and draws were frozen
before its new calculations.

CONTENTS AND VERIFICATION
code/attribution_reliability.py is the executed numerical implementation.
CODE10_PROTOCOL.json and PROTOCOL_FIXED_AT.json record the fixed specification;
frozen/ contains the exact partitions and common feature IDs. results/ contains
Part A moments, Gram matrices and correlations; Part B per-draw scores,
coefficients/ranks, associations and matched contrasts; and complete summaries
with finite/unavailable counts. Gene-by-draw arrays are replayable from pinned
inputs and are not redundantly stored. logs/ records execution and code hashes.
The archive root contains PLANT_STRESS_COMBINATIONS_CODE10_RESULTS.xlsx.
report/ contains the scientific figure CODE10_attribution_repeatability.pdf/.png,
plus report code and checks.
review/ contains the precomputation design review, theory/reference review,
interpretation review, synthetic identity checks and independent numerical audit.

Synthetic edge-case tests passed. All original-pool fits reproduced CODE8.
The independent audit imported no production implementation and refitted all
1,024 sides using QR-based calculations. Its 80,479 numeric/exact checks passed;
maximum absolute discrepancy was 2.13e-13. All parent pins remained unchanged.
The executed analysis took approximately 44 seconds in the available environment.

REPRODUCTION
Extract the previously delivered CODE8 archive alongside this CODE10 folder,
preserving the stress_combinations_code8/ and stress_combinations_code10/ names.
CODE8 supplies the processed matrix and six frozen assets; it is not duplicated.
Archive: PLANT_STRESS_COMBINATIONS_CODE8_ANALYSIS_RECORD.zip
SHA256: ee225c95082e3a029d3063761463110ad7e60442263ade042f2f9e9d8155a549
Download: https://chatgpt.com/library?file_id=libfile_1c6ec5d8ae0c81918a16440f444c028a&account_id=personal
This CODE10 bundle also includes the two small CODE9 dependencies at their
original relative paths: stress_combinations_code9/CODE9_PROTOCOL.json and
stress_combinations_code9/results/cross_genotype_summary.csv. The full CODE9
archive is unnecessary for numerical replay.

Use existing Python with NumPy and pandas; no environment creation, installation
or SciNet is required. Independent audit additionally uses SciPy; plotting uses
matplotlib. From the common parent directory, set OPENBLAS_NUM_THREADS,
OMP_NUM_THREADS and MKL_NUM_THREADS to 1. PowerShell example:

$env:OPENBLAS_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$py=Join-Path $env:USERPROFILE 'anaconda3\python.exe'
& $py stress_combinations_code10/code/attribution_reliability.py selftest
& $py stress_combinations_code10/code/attribution_reliability.py run --protocol stress_combinations_code10/CODE10_PROTOCOL.json --output stress_combinations_code10/replay_results

Use a fresh output directory; do not overwrite the included results. The run
rechecks the source, frozen assets, protocol and partition hashes. Windows/GitHub
packaging and local confirmation remain deferred until the pipeline is finalized.
