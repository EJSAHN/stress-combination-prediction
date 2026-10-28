PLANT STRESS COMBINATIONS — CODE7 UNCERTAINTY ANALYSIS

DECISION
Retain CODE6 fixed predictive performance. The signal survives expansion through observed lower-order replicate variation. The stronger claim of a genuinely additional biological response direction remains unresolved: one fixed conditional noise scenario reaches the observed statistic within its upper simulated range. This does not prove the signal is noise.
CODE7 is post-result sensitivity analysis. CODE6 target values were already known; CODE7 settings were fixed before these new calculations. This is not new external validation or a calibrated biological hypothesis test.

EMPIRICAL EXPANSION
WT outside-V3 MSE-gain contribution: 0.062273019
Outside-V11 contribution: 0.055673527
Original outside prediction-gain contribution retained: 89.402%
Outside-V11 target-library gains: 0.047636669, 0.058754359, 0.060629553
Outside-V9 disjoint-control margins: 0.218381277, 0.216165864, 0.198171554
All six genotype/filter scenarios met the fixed empirical sign rule. V3 spans three single-condition means; V9 spans nine individual single-stress libraries; V11 adds control replicate variation. Pair/target values did not choose the bases. The exact gain partition is not biological variance explained.
Disjoint-control checks use V3/V9 only: V11 would include held-out control directions and break the covariance argument. The diagnostic uses another fixed gene panel and population-geometry estimand. Overlapping folds are not independent replications.

CONDITIONAL SIMULATION: WT, UNIT NOISE SD
Latent model | diagonal covariance fraction | conditional 95th percentile | at least observed /1000
single_mean_span | 0 | 0.048011823 | 18/1000
joint_rank3 | 0 | 0.064190945 | 56/1000
single_mean_span | 0.5 | 0.038408589 | 1/1000
joint_rank3 | 0.5 | 0.051186611 | 12/1000
single_mean_span | 1 | 0.028418540 | 0/1000
joint_rank3 | 1 | 0.040391814 | 0/1000
Observed WT statistic: 0.062273019.
Joint-rank-3/empirical-covariance q95=0.064191, with 56/1000 outcomes at least as large as observed (Monte Carlo SE of the fraction about 0.00727). This is a borderline conditional range comparison, not a biological P-value, noise attribution or calibrated acceptance/rejection. The corresponding mutant secondary scenario has 61/1000.
All 36 scenarios were run: two genotypes, two latent-mean constructions, three covariance mixtures, three SD scales, 1000 draws each. Common errors across nulls/scales mean these 36,000 scenario outcomes are not independent biological experiments. No scenario was selected to rescue the claim.
The single-span null projects noncontrol means into the observed single span. The joint-rank-3 null fits three directions to all seven noncontrol means, including the already-seen target. Its projected singles span all three directions. This is an outcome-seen nuisance construction, not a predictor.
Each condition covariance mixes its rank-at-most-two empirical residual covariance with a diagonal preserving gene variances (fractions 0, 0.5, 1). Noise SD scales are 0.5, 1, 2; scale 2 gives four times covariance. Three RNA pools cannot identify full cross-gene covariance. Neither orientations nor covariance choices exhaust the possibilities.
Mean-fit MSE divided by estimated mean-noise trace is a scale diagnostic, not a calibrated goodness-of-fit statistic. Matching scalar T does not show that the null can generate the whole expression dataset or its mean-residual directions: low-rank covariance support matters.
Three scores share each draw: refitted span/noisy control, latent span/noisy control, latent span/true control. Shared control error can create positive gain even with a known span. Their paired differences separate those effects under the imposed model.
Eight condition-mean draws suffice because target replicate variance cancels in mean error differences. These T-only simulations do not calibrate disjoint-control margins. Gaussian working log-values are untruncated and are not synthetic counts.
Zero exceedances means 0/1000 (resolution 1/1000), not zero probability. Plug-in Monte Carlo SE is zero at boundary counts without eliminating uncertainty. Strict negative-log diagnostics can count floating-point roundoff around zero; they do not establish biological invalidity.

WHAT REMAINS SUPPORTED
CODE6 fixed predictions and errors are unchanged. Useful prediction information lies outside the particular estimated single-response span and mostly survives the observed lower-order replicate cloud. These n=3 data do not establish that the unknown true single-response span must lack that information.
Do not claim pure higher-order interaction, a new biological dimension, a universal prediction law or that noise has been ruled out. Further tuning against the open target cannot create independent evidence. A stronger biological geometry claim needs more informative independent measurements or a separately justified uncertainty model.

FILES AND REPRODUCTION
The workbook preserves all empirical scenarios and conditional summaries. The archive includes raw simulation draws, original expression matrices and required CODE6 dependencies, exact English code, fixed settings, seeds, logs and independent audits. No new environment, installation or SciNet was used.
Use fresh output paths with EXECUTION_RECORD.json commands. This is an executed record; final GitHub/Windows packaging remains deferred.
