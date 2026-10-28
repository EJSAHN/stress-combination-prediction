PLANT STRESS COMBINATIONS CODE5
Completed 28 October 2026. GSE147962 discovery analysis and coffee-paper bridge.

DECISION
Keep the quantitative biological question. Do not advance the low-common-
response selector as a superior new method. Its primary donor choices are
identical to ordinary distance-to-control selection for every target.
The more informative result is that the previously observed four-way
prediction gain has two substantial components: directions represented by
weaker lower-order responses, and directions outside the full fitted
lower-order response span. Both remain positive in a disjoint-control check.

The analysis and independent numerical checks were performed here. The user
does not need another Windows run. GitHub and local execution packages remain
deferred until the pipeline and scientific narrative are settled. No package
installation, new environment, SciNet task or external expression download
was required. The attached coffee paper was read, not edited.

TARGET-EXCLUDED DONOR SELECTION
Five targets were considered: four triple-stress cells and the four-way cell.
Each triple could use its three immediate pair donors; the four-way cell
could use its four immediate triple donors. Gene filtering and basis fitting
used only the 33 control/single/pair libraries. The generator read only these
libraries plus permitted donors for the current target, never its target
expression columns. A separate command evaluated saved predictions.

All values had already been examined in earlier work. This is exploratory
code-level separation, not a new blind or external validation. Source-wide
normalization may couple all 48 libraries. Pair donors for triple prediction
also contributed to the fitted basis; triple donors for four-way prediction
did not. Tasks and information sets are therefore not interchangeable.

The twelve fixed models included control, additive and pairwise baselines,
all-donor and half-shrunk means, six one/two-donor selection rules, and one
residual-only prediction. Rank two was primary; rank one was a sensitivity.
Primary comparison: common_low2 versus raw_low2 and the all-donor mean.

Primary mean MSE across three target RNA pools, 21,495 genes:
Target       Control        All donors      Low-common two = Low-raw two
PQ-HL-HS      0.117474        1.204457                 0.951077
S-HL-HS       1.973082        0.631740                 1.178980
S-PQ-HL       0.180685        0.556512                 0.489048
S-PQ-HS       1.899995        0.481705                 0.847601
S-PQ-HL-HS    0.171270        0.590791                 0.129619

S = salt; PQ = paraquat; HL = high light; HS = heat stress. Smaller MSE is
better. All twelve models, both ranks, MAE and individual-library errors are
retained. Low-common two beats the all-donor mean in three targets and loses
in two. This is not a universal selection rule. Low-common one chooses
PQ-HL-HS for the four-way target and gives MSE 0.211924, worse than control.
The four-way two-donor prediction is exactly the earlier CODE3 quiet centroid,
not another independent performance result. Its primary membership is stable
under all 48 conditional single-donor-pool omissions across these tasks;
the fitted basis/filter/control were held fixed. Stability does not establish
predictive validity, as the unfavorable target results demonstrate.

WHAT THE EARLIER PREDICTIVE GAIN CONTAINS
The coffee comparison exposed a necessary distinction: a rank-two projection
remainder is not a gene-wise factorial interaction. We therefore decomposed
the unchanged CODE3 quiet-centroid effect using the complete fitted span of
the ten nonzero lower-order condition-mean displacements. Numerical rank was
ten; no target was used to choose rank or directions.

Component                                  Mean MSE gain over control
First two lower-order directions                    0.0002575034
Remaining lower-order directions 3-10               0.0200225740
Outside the full fitted lower-order span             0.0213708661
Total gain                                          0.0416509435

These gains add because the three gene-space components are orthogonal.
The latter two each give a positive gain in all three target RNA pools.
Roughly half the observed gain lies in each of those two components. This
is a decomposition of one realized prediction result, not a proportion of
biological variance, a causal partition, or a new model-selection result.

After those component results were known, a separately recorded amendment
extended the CODE4 disjoint-reference diagnostic to the complete span. It
used 21,487 genes filtered without controls, one control library for fitting
the basis and the other two as distinct donor/target references, rotating
through all three choices. Library covariance corrections were retained.

Outside-full-span inner products were 0.0275773, 0.0228990 and 0.0257001.
Corresponding population proximity margins were 0.0202795, 0.0167103 and
0.0179836. The remaining lower-order directions also had positive margins
in all three splits. Thus both parts survive this targeted reference-noise
check. The splits overlap; their margins are not independent replications,
P values, or the original CODE3 prediction gains.

The full span is estimated from noisy condition means. Basis estimation
error can place genuine lower-order signal outside that estimated span.
Conversely, higher-order interactions can lie inside it. Therefore the
accurate description is 'outside the fitted lower-order expression span',
not 'pure higher-order interaction' or proof of a new regulatory mechanism.

CONNECTION TO THE COFFEE PAPER
The attached paper is Genotype-specific non-additive heat-drought responses
and recovery hysteresis in coffee, by Baek, Park, Upadhyay, Kim, Meinhardt
and Ahn. Its gene-wise normalized-log expression contrast is
SWD_s - WW_s - SWD_T25 + WW_T25,
supported by genotype-specific DESeq2 water-by-stage interaction tests.
It provides a natural biological antecedent: departures from additive
responses vary with genotype and stress/recovery stage.

Our prediction and geometry analysis extends the question of how observed
responses inform other combinations. Its projected remainders are not the
same estimand as the coffee interaction contrast. The coffee recovery
interpretation also cannot be transferred to one-endpoint GSE147962 data.
Coffee's two-water-regime/four-stage design is suitable for a later separately
defined application, but cannot directly validate the four-binary-factor
donor-selection task. Its data are already known to the team.

The user's accepted/proofing status at npj Science of Plants is consistent
with the narrowly checked administrative record. No final publisher DOI was
located in this search; no November publication date is assumed. A primary
Crossref record verifies the Research Square preprint DOI
10.21203/rs.3.rs-9849440/v1, posted 12 June 2026. Official npj instructions
allow accepted papers and recognized preprints in references. Use appropriate
in-press/preprint status until final bibliographic information is available.
The preprint abstract lacks R1's later STAR/QC validation; do not cite that
version as evidence of those added checks.

WORKING SCIENTIFIC DIRECTION
The evidence supports studying the order and context at which transferable
expression information becomes observable. Pair-derived predictions poorly
capture the quiet triple states, whereas their observed triple centroids
contain information useful for the four-way state. Part of that information
lies outside the fitted lower-order mean span. This is a within-experiment
finding, not an established general rule or a claim of method superiority.

Keep GSE147962 as the discovery dataset. A useful next test would apply fixed
component definitions and prediction comparisons in an independently chosen
factorial experiment, including the simple centroid and control baselines.
Do not pick a new rule after looking at an external target. Tomato GSE294808
remains metadata-only in this project; its expression has not been opened.
Coffee can motivate the crop context without forcing its different design
into the same validation task.

VERIFICATION AND FILES
Independent reconstruction checked 120 prediction vectors, 32 donor-score
rows, 120 donor selections, 576 omission records, 360 library errors and
120 target means. The largest numerical difference was 6.75e-14. Replacing
each target's three columns with invalid strings, NaN and extreme values
left selection and predictions unchanged. Separate reconstruction confirmed
the component and disjoint-reference results (largest fixed-component
difference 3.64e-13; disjoint difference 3.85e-16). These checks start from
author-processed expression, not raw sequencing reads.

The Excel workbook summarizes the primary comparison, all model means,
subspace contributions, disjoint-control checks and source/method notes.
frozen/ preserves predictions, donor scores, allowed columns and selections.
results/ preserves all model/library comparisons.
subspace_resolution/ preserves the additional protocol, timed amendment,
executed scripts, component metrics and assets.
review/ preserves independent numerical and input-exclusion audit evidence.
references/ preserves the coffee comparison and primary DOI metadata.

The source expression matrix and previous-stage pinned assets are not
duplicated in this analysis record. Scripts identify them through command
arguments and recorded hashes. Final public/local reproducibility packaging
is intentionally deferred. The attached manuscript and its full extracted
text are not redistributed in this archive.
