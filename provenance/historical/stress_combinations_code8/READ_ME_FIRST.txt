PLANT STRESS COMBINATIONS — CODE8 CITRUS VALIDATION

DECISION
The fixed pair-response average met the predeclared descriptive prediction criterion in BOTH citrus genotypes. This supports transfer of the prediction rule to this crop experiment. It does not settle CODE7 uncertainty about the unknown true single-stress subspace.

PRIMARY RESULTS (FPKM filter 1; error averaged over genes and three target RNA pools)
Genotype | genes | pair-average MSE | best ordinary comparator MSE | relative MSE reduction
Car | 17470 | 0.212010622 | 0.309363573 (second-order factorial) | 31.469%
Cleo | 16247 | 0.279195458 | 0.486561664 (second-order factorial) | 42.619%

Both primary genotypes outperform every one of the seven ordinary comparators on mean MSE; the pair average also beats all seven within EACH target library. All six genotype/filter scenarios meet the fixed mean-comparison and per-library-control criterion. Thresholds 0.5 and 2 are sensitivity analyses, not extra independent experiments or rescue choices.
Compared with the pair-trained three-weight single-response model, MSE is lower by 48.549% (Car), 52.914% (Cleo).

WHAT WAS FIXED BEFORE OPENING
Protocol, thresholds, predictors, targets and decision rule were fixed before acquisition. The file was then acquired and its header and identifiers inspected. GEO called the processed layer RPKM, but the deposited file has 48 FPKM columns plus Read_Count and TPM layers. An explicit amendment chose only FPKM, preserving all numerical choices, before numerical expression parsing. All48 sample names were independently matched to GEO and original SRA FASTQ names. No raw reads were processed.
Each genotype was fitted separately using its control, three single and three double conditions. Genes were selected using only 12 control/single libraries. Three double-condition means received equal weights. Eleven predictors, three observed spans and all six scenarios were saved and independently checked before the six triple-condition numerical columns were read. Replacing those target fields with NaN left every frozen array bitwise identical.
Additional baselines used the double conditions to fit either one global scale or three stress-specific coefficients for the single responses. No intercept, bounds, gene-wise fitting or target tuning was used. This is prediction-rule transport, not transfer of Arabidopsis gene coefficients or an orthology model.

SECONDARY GEOMETRY
Gain outside the full estimated single-mean span remained positive, including after expansion with nine single libraries and control replicate variation. Outside-V11 mean gains were 0.169562220 (Car), 0.284781417 (Cleo). These concern an observed noisy span and are not proof of a new biological dimension.
The target-fitted V3 oracle is deliberately outcome-seen. It fits the target mean inside the observed affine single-response span, providing an optimistic minimum mean error in that span, not a legitimate prospective predictor. Oracle mean errors were 0.379155282 (Car), 0.566821099 (Cleo), both larger than the fixed pair-average errors. No oracle coefficient enters any prediction.

BIOLOGICAL AND NOVELTY LIMITS
GSE203331 contains leaf RNA from three biological pools per genotype/condition, each pooled from at least10 plants; the RNA sample size is three. Carrizo and Cleopatra belong to one experiment and were aligned to the same sweet-orange reference. Experimental crop and provider differ from the earlier Arabidopsis work, but investigators overlap. This is not fully independent-laboratory replication.
Heat and high irradiance ran over five days; water stress was applied for the last24h. This is an endpoint RNA experiment, not recovery or hysteresis. Papers disagree about baseline growth irradiance (312 versus100); applied stress doses and condition identities are consistent. No inference uses that disputed baseline dose.
Weighted integration, stress dominance and attenuation have prior plant literature (especially Tan et al.2023, Nature Communications14:986, doi:10.1038/s41467-023-36517-w). This simple average is not a new mathematical algorithm. CODE8 contrasts are not exact reproductions of that study. A publication claim must center on the prospectively fixed combination-order prediction and its limits, not new arithmetic, pure high-order interaction or a universal law.
No gene-as-replicate P values, gene bootstrap, biological significance threshold or noise-identifiability claim is made. CODE7 limitations remain. Author preprocessing can couple information before our column-level separation; this was not a blinded raw-read pipeline.

REPRODUCTION
Use the included original processed matrix, mapped design, protocol and pre-numeric amendment. The English script has selftest, freeze, poison-check and evaluate commands; use fresh output directories. See EXECUTION_RECORD.json and logs. It needs only an existing Python with numpy/pandas for analysis, plus matplotlib for the report figure. No environment installation or SciNet was used. Final Windows/GitHub packaging is deferred until the manuscript claim is settled.
Independent verification uses QR rather than production SVD and independently assembled coefficient equations. The review records identify precise checked quantities and a repaired incomplete reviewer asset. Frozen production predictions were not changed.

SOURCES
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203331
https://doi.org/10.1111/ppl.13809
https://doi.org/10.1093/hr/uhad107
https://doi.org/10.1038/s41467-023-36517-w
Additional prior-art and sample-provenance details are under references/ and metadata/.
