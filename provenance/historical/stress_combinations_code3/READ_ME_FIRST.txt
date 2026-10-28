PLANT STRESS COMBINATIONS CODE3
GSE147962. Completed 28 October 2026.

RESULT
Residual expression from the two previously identified quiet triple conditions
helps predict the numerically unopened four-way condition in this experiment.
The simple quiet-state centroid performs slightly better than the residual-only
prediction. Projection therefore provides an explanatory decomposition here;
it has not produced a superior new prediction method.

PRIMARY FOUR-WAY COMPARISON
Errors are averaged over three target RNA pools and 21,495 training-filtered
genes on the log2(1 + author-normalized count) scale. Lower is better.

Predictor                       MSE          MAE
Control mean                    0.171270     0.264924
Quiet common component only     0.171012     0.263682
Quiet residual component only   0.129876     0.234115
Full quiet centroid             0.129619     0.230367

The residual-only predictor reduces MSE by 24.17% and MAE by 11.63% relative
to the control-mean predictor. Both losses improve in each target pool. The
full quiet centroid reduces MSE by 24.32% and has slightly lower MAE. Almost
all the gain lies in the residual component under the fixed rank-two geometry.

All sixteen prespecified models remain in the workbook and CSVs. Eleven are
the unchanged CODE2 lower-order models fitted using 33 libraries. Five also
use previously explored triple conditions. These families have different
information available; this is not an equal-training-data algorithm benchmark.
The quiet conditions, PQ-HL-HS and S-PQ-HL, were selected after CODE1/CODE2
discovery. The comparison is not evidence of universal transfer to unseen
species, stress types, experiments or doses.

FUNCTIONAL RESIDUALS
All fifteen externally annotated GO programs and five fixed rank/filter/
annotation scenarios were retained. The full residual-program profile has
cosine similarity 0.8608 between the four-way target and the quiet centroid;
the five scenarios range from 0.8608 to 0.8847. Primary mean-effect signs agree
for 11 of 15 programs. Heat-response, oxidative-stress and protein-folding
gene sets retain negative residual means. All discordant programs are reported.

Directional agreement is not accurate magnitude transfer. Primary residual
program-profile RMS prediction error is 0.06335, versus 0.06262 for the raw
program profile. The target protein-folding residual mean is -0.309, versus
the frozen quiet-centroid prediction of -0.183. No improvement in absolute
program-profile prediction is claimed from the higher cosine.

Matched random gene sets are a descriptive check only. For the primary
scenario, 1,000 sets per program match retained size and joint bins of lower-
order mean expression and within-condition projected residual variance.
The draws were fixed before target opening. All fifteen comparisons were
feasible. These sets do not preserve biological coexpression or GO overlap;
their empirical percentiles are not P-values or mechanistic enrichment tests.

POST-OPEN REPLICATE-COUNT DIAGNOSTIC
The quiet centroid uses six RNA libraries, whereas the control mean uses
three. After the primary results were known, a separate diagnostic compared
all eighteen balanced three-library quiet subsets, taking one library from
one quiet condition and two from the other. Every subset improves both MSE
and MAE over the control predictor for each of the three target libraries.
Subset mean MSE ranges from 0.11887 to 0.15669, compared with control 0.17127.
This addresses the simple extra-library averaging explanation. The subsets
overlap, are not eighteen independent validation experiments, and were not
used to replace or select the frozen primary predictions. Plug-in variance
estimates based on three pools per cell remain uncertain.

WHAT WAS FIXED BEFORE OPENING
CODE3_PROTOCOL.json records the formulas, sixteen models, exact source hashes,
gene filters, ranks, GO memberships, random-set design, metrics and limitations.
Both freeze stages parsed only the 45 non-four-way expression columns. Exact
prediction vectors, bases, gene sets and matched draws were saved with hashes.
OPENING_RECORD.json records the verified hashes before the first numerical
four-way read. Separate evaluation commands then read the target columns.
No primary model, rank, filter or membership was changed after opening.

Independent numerical reconstruction confirmed the sixteen prediction vectors,
all 192 model/library metric rows, all 1,200 program/condition rows, pool-based
uncertainty, and target matched-set percentiles. Exact model/program gene order,
control mean and primary basis agree. Review evidence is under review/.

INFERENCE BOUNDARIES
There is one new factorial cell with three biological RNA pools. Published
qualitative four-way patterns were already known, and the source normalization
may use all 48 libraries. This is a numerically withheld within-experiment
check, not independent preprocessing or fully blinded external validation.
All 48 libraries from the 16 conditions have now been analyzed; future analyses cannot
treat this four-way condition as a new unopened test.

GO scores are equal-weight expression summaries of unsigned member genes, not
pathway activity, resource allocation or flux. A residual GO score can depend
on genes outside that set through the global projector. Residual/raw RMS for
a gene subset can exceed one; it is not an unexplained-variance fraction.
Pool-level SE is conditional on a fitted basis and omits its estimation
uncertainty. Leave-one-pool ranges are sensitivity ranges, not confidence
intervals. Shared-control profile similarities are descriptive.

Salt/PQ remain present during the acute heat/light challenge. There is no
washout, temporal-order or dynamical phase-transition experiment. Paraquat-
light photochemistry, bulk composition and relative RNA normalization remain
alternative explanations. The original paper already described attenuated
common responses; novelty is not established merely by observing attenuation.

NEXT QUESTION
Keep this dataset. The candidate question is whether attenuated common stress
expression coexists with a transferable, functionally nonuniform residual
response. Further work on these data is exploratory. A claim of a general
biological rule or a useful new method requires an independently chosen test
experiment and a comparison with the simple quiet-centroid baseline.

FILES AND LOCAL REPRODUCTION
Start with PLANT_STRESS_COMBINATIONS_CODE3_RESULTS.xlsx. Complete numerical
tables are in results/. Frozen vectors and assets are in frozen/. Original
data, membership tables, code, protocols and verification are bundled.
reproduction/README.txt contains the complete Downloads-to-PowerShell command.
The analysis has already run. Optional local reproduction uses an existing
Python with NumPy/pandas; no installation, new environment, network or SciNet
is required. Fresh runs and their result ZIPs stay under the D: project.
The Windows wrapper was reviewed but not executed in this Linux environment.

SOURCES
GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147962
Original study: https://doi.org/10.1111/nph.17232
TAIR annotation archive: https://zenodo.org/records/17371665
GO release: https://release.geneontology.org/2024-11-03/ontology/go-basic.obo
CODE3 uses the unchanged, hash-pinned CODE2 functional memberships.
