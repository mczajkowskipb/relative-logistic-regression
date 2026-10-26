# CRIT Stage H4A — frozen-blind LSH stress test; no retuning

Protocol fixed **2026-10-26**, before looking at H4A outcomes. Tracking [Stage H4 #2](https://github.com/mczajkowskipb/relative-logistic-regression/issues/2). Git reference for the *frozen* H3 estimator: commit \`cd1b3038f3b662216a19655faccbec394de5900e\`, source blob SHA \`dcffff64207ff71d7c52872a7c1f62e04617ceaf\` at \`scripts/crit_stageH3_blind_lsh_pilot.py\`.

## Evidence and question

Stage H3-fix1 identified planted exact two-rank-pair XOR interaction **19/20** times in the top six across p32–512, at nD120, interaction logit coefficient 2.4, using no known truth anchor but supervised D labels. The independent auditor recomputed all 240 saved final D scores. All prior Stage H3 data are **exposed**.

H4A asks whether **the exact same method** (18 sample hash bits; 24 tables; 4096 scoring batch; exact top-six, no oracle) retains promising planted-XOR recovery when *p, nD and signal strength* change. Do not alter these parameters based on H4A results.

## Predeclared experimental factors

- p∈{32,128,512}; nD∈{60,120,240}; XOR logit amplitudes a∈{1.0,1.6,2.4}, and a=0 Bernoulli null. Eight new independently seeded D datasets per (p,nD,a) combination. **3×3×4×8 = 288 studies in 36 cells**.
- Seed \`20280617+nD*101+p*313+amplitude_index*100003+replicate*47\` with \`amplitude_index=0,1,2,3\` for 0,1.0,1.6,2.4. Data X are independent Gaussian features, generated using NumPy \`default_rng\`; label distribution \`Pr(Y=1)=logistic(a*sign(X1−X2)*sign(X3−X4))\`.
- Frozen H3 estimator file is imported **unchanged**; only the module's global \`N\` is assigned the current nD, so that signature shape and score normalization follow the input D sample count. This is a compatibility adjustment to experiment sample size, **not** a scoring/LSH algorithm update.
- All D samples available to hash signatures and rank candidate interactions, no V or test E; true pair identities used only in scientific evaluation and the previous H3 auditing instrumentation.
- Primary descriptive endpoint: **both exact planted rank pairs as an unordered XOR product present in D top-six**, with seedwise hit and hash-collision flags. Under null, the same two pairs are **arbitrary reference labels**, not causal truth or false acceptance. Also report incomplete lists as legitimate retrieval failures, scores, number of collisions/scored candidates, elapsed time and process peak RSS.
- Wilson 95% intervals are descriptive uncertainty around a cell's 8 independently drawn samples, not confirmatory inference on a universally applicable method or clinical effect.
- No p1024 because the current H3 hash-join can create huge repeated candidate fanout. H4A explicitly limits p at 512 and should not be described as transcriptome-scale.

## Strict version/control and PASS requirements

- Recheck the frozen H3 estimator git blob hash *byte-exact* before running.
- GitHub-hosted Actions, public synthetic only, no SSH or patient data.
- Distinct output \`results/crit_stageH4A_blind_robustness_v1/\`; never overwrite prior CRIT files or artifacts.
- Independent Python evaluator regenerates X and y for all 288 new seeds **without importing the estimator**, independently recomputes every published top-six pair-product score, verifies disjoint genes, no duplicated candidates, truth flags and source metadata.
- Required \`STAGEH4A_288_FROZEN_LSH_STUDIES=PASS\`, \`STAGEH4A_36_PREDECLARED_CELLS=PASS\`, \`STAGEH4A_288_SOURCE_SEED_AND_FEATURE_AUDITS=PASS\`, \`STAGEH4A_36_CELLS_AUDITED=PASS\`, and \`STAGEH4A_BATCH_EXIT=0\`.
- The independent audit logs the dynamically computed number of candidate scores rather than presuming all shortlists have six members.

## Interpretation and decisions

**GO to an isolated full original-RLR D/V/E integration prototype only if** blind retrieval works acceptably over a meaningful range of *weaker* signals and realistic p/n combinations and resource consumption is bounded, with the original-RLR source and patient-level baseline parity preserved. There is no predeclared numerical success threshold for all those cells; decisions remain *exploratory feasibility* and should not be misreported as confirmatory inference. If retrieval fails as conditions get harder, report failure transparently, and design a new method in a **separate new protocol**, never tweak hash bits/tables based on these same simulations and then claim independent replication.

For eventual scientific efficacy or null-error control, this H4A test is NOT the answer: no independent E generalization, V selection, RLR offsets, biomedical cohorts, or alpha calibration is evaluated. Full pipeline null simulations (hundreds of independent runs) and held-out E results are future work.
