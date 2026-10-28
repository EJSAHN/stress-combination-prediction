PLANT STRESS COMBINATIONS — CODE9 GENE ATTRIBUTION

DECISION
The frozen prediction finding remains. Its advantage over learned single-response weights is spread across many observed features and survives deletion of the largest10% contributors. The advantage over second-order factorial is more concentrated. An outcome-fitted global rescaling of that second-order vector still does not reach the original frozen centroid accuracy. Functional patterns depend on genotype and comparator; there is no established universal pathway mechanism.

STATUS
CODE9 is post-result explanation of CODE8, not new held-out validation. All CODE8 target values and aggregate outcomes were already known. The attribution protocol and external annotation categories were fixed before CODE9 gene gains were computed. No parent predictions, parameters, gene panels or decisions were changed.

PRIMARY GENE DISTRIBUTION
Genotype | comparator | fraction improved | mean gain after removing top10% | smallest ranked deletion erasing advantage
Car | pair_trained_factor_OLS | 65.255% | 0.041553170 | 17.533%
Car | second_order_factorial | 57.315% | -0.044073682 | 4.419%
Car | within_V3 | 71.958% | 0.080119024 | 34.322%
Cleo | pair_trained_factor_OLS | 70.604% | 0.087934309 | 24.835%
Cleo | second_order_factorial | 60.787% | -0.030712114 | 7.004%
Cleo | within_V3 | 78.328% | 0.118320591 | 43.023%

Gain is comparator squared error minus centroid squared error, averaged over the same three RNA pools. Its mean across features exactly recovers the corresponding CODE8 MSE difference. Positive and negative contributions are retained. Feature_GID is the original scoring unit; individual features are not independent biological replicates.
Top-k removal uses descending observed mean gains, lexical-ID ties and exactly ceiling(q*G) removals at1%,5%,10%. Only the scoring set changes; predictions are not refitted. It minimizes retained mean-gain sum for that fixed number of deletions, but not separately for each RNA pool. The minimum removal fraction is an outcome-selected score-concentration diagnostic, not estimator breakdown, experimental gene deletion or inferred biological sparsity. Net-gain retention can be negative; net-gain shares can exceed100% when other features contribute negatively.

GLOBAL RESCALING DIAGNOSTIC
Genotype | model | fitted scale | rescaled MSE | original frozen centroid MSE
Car | full_pair_centroid | 1.281760058 | 0.157064704 | 0.212010622
Car | second_order_factorial | 0.800849994 | 0.244452433 | 0.212010622
Car | pair_trained_factor_OLS | 0.867980214 | 0.391174803 | 0.212010622
Cleo | full_pair_centroid | 1.244075689 | 0.225264922 | 0.279195458
Cleo | second_order_factorial | 0.755272836 | 0.352849123 | 0.279195458
Cleo | pair_trained_factor_OLS | 0.879512220 | 0.573182843 | 0.279195458
The scale minimizes error using the already observed target mean. It is a diagnostic lower bound along each observed prediction direction, not an independently evaluated predictor. Because even the optimally scaled second-order vector has larger error than the original centroid in both primary genotypes, one global scale alone cannot explain this observed advantage. It does not identify an unknown true biological subspace or remove CODE7 measurement uncertainty.

SHARED FEATURES AND FUNCTIONS
The exact common primary universe has 15,466 features. The strict descriptive set contains 604 features with positive gain in every target RNA pool in both genotypes against all three primary comparators. 304 have admissible biological-process annotation. These are outcome-selected examples, not validated markers or differentially expressed genes.
second_order_factorial: cross-genotype Pearson=0.118982; Spearman=0.063295; all-three-pools-positive shared count=3569.
pair_trained_factor_OLS: cross-genotype Pearson=0.153887; Spearman=0.053900; all-three-pools-positive shared count=5333.
within_V3: cross-genotype Pearson=0.268145; Spearman=0.118451; all-three-pools-positive shared count=6776.
The low cross-genotype gain-rank agreement does not prove distinct biological mechanisms: baselines were fitted separately, sampling noise remains, and the comparison is conditional on the fixed observed predictions.
Translation-associated features show a relatively high fraction of positive gains versus learned single-response weights in both genotypes. This is comparator-dependent: Carrizo translation has negative net gain against second-order factorial. Cleopatra protein-folding features have negative net gain versus learned single-response weights, despite a majority having positive individual gains. Functional activation, uniform benefit or a conserved protective pathway must not be inferred from prediction-error gains.
The original citrus study already discussed RNA translation and unfolded-protein responses (doi:10.1093/hr/uhad107). Their appearance here is not a new pathway discovery. The present contribution concerns where a fixed prediction rule gains and loses accuracy.

ANNOTATION AND UNITS
NCBI AR103 GAF annotations were joined by exact Entrez ID to the older source expression-reference identifiers. No symbol guessing or orthology remapping was used. GAF and ontology dates match2023-07-27; all63,808 GAF records carry computational IEA evidence. GO.db3.18.0 supplies the matched ontology without installing R or packages. Only is_a/part_of ancestry was used, excluding regulatory edges. All source terms resolved.
Functional categories were fixed to protein folding, photosynthesis, translation, cellular oxidant detoxification, transmembrane transport and cell-wall organization/biogenesis. Categories can overlap, so their contributions cannot be summed as a partition. Any-GO and BP-only coverage are separate. BP comparisons use the unambiguous BP-annotated background; unknown, non-BP and other-BP features remain in complete tables.
Twenty-three source Entrez IDs map to46 feature rows. All those rows were excluded from functional summaries; a separate scoring-only exclusion check preserves the original benchmark. No representative was selected. The oxidant-detoxification group has only10/9 primary-panel members and is labelled small. There are no GO enrichment P values or genes-as-biological-replicates tests.

REPRODUCTION AND CONTENTS
The exact executed numerical code is code/gene_gain_attribution.py. See logs/PIPELINE_EXECUTION_RECORD.json for commands and results/RUN_STARTED.json for hashes. Outputs include all402,616 feature-comparator rows, complete deletion curves, shared-set membership, annotations, examples and scalar diagnostics. The scientific figure has PDF and300dpi PNG versions. The Excel file reports all fixed summary scenarios, including negative findings.
CODE9 requires the previously delivered CODE8 record extracted alongside this folder. That existing archive contains the original processed matrix and six frozen prediction assets; it is not duplicated in this bundle. Parent archive SHA256: ee225c95082e3a029d3063761463110ad7e60442263ade042f2f9e9d8155a549
Parent archive: https://chatgpt.com/library?file_id=libfile_1c6ec5d8ae0c81918a16440f444c028a&account_id=personal
The annotation GAF, GO.db source archive, extracted graph/terms and exact memberships are included. GO.sqlite is an excluded reproducible cache, obtainable by extracting the included GO.db archive. Full raw-read reprocessing was not performed. Existing Python with NumPy/pandas is sufficient for numerical execution; SciPy is used by the independent audit and matplotlib by report plotting. No environment creation or package installation occurred.

VERIFICATION
Independent scripts verify the raw annotation join/graph, every402,616 contribution row, ordering and deletion quantities, cross-genotype sets, direction-oracle algebra, summaries and functional/QC tables. Maximum numerical discrepancy was1.36e-12. Parent hashes remain unchanged. Exact workbook-cell and visual checks are recorded under report/. Final GitHub/Windows confirmation remains deferred.
