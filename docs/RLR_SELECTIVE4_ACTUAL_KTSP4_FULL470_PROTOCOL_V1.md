# Four-pair selective one-sided rank voting — original RLR 11-cohort matched protocol
Date 2026-10-27. **The author's original RLR is UNPUBLISHED and all frozen manuscript/core files are unchanged.**

## Scientific motivation and guardrails

SAVS V1–V4 completed on original paper's 11 CompCancer firstfold inner D/V cohorts, with independent patient-level audits. Actual 11-cohort V1 original rank-margin `m10` tree F1 .7660 versus matched global four-rule F1 .8415; under fixed `m10pct=max(10,ceil(.1*(p-1)))` root tree F1 .7669 versus global four-rule F1 **.873149**, Macro-F1 **.891151** and mean **1.747361** active votes. Original RLR50 manuscript's CompCancer structural active means are 28.09 GLOBAL nonzero fitted coefficients and 25.55 LOCAL nonzero per sample under the ternary m10 rule, nominal 50. The V1 global 4-rule method is **not** yet proven competitive on the outer E repeated 470-fold protocol.

V3 fixed supervised-root/pair shortlist leakage in V2's D CV via nested S-only root nomination, H-only selection, final D refit, V audit. V4 allowed root gene/anchor reuse in different child comparisons and independently proved patientwise invariance. [V4 all 11 CI #40567165512](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40567165512) **PASS**, 132 root-method V records independently replayed. Mean V F1 m10pct: hub-sharing conditional .774187, no-hub honest tree .784761, root-first tree .766859, global FOUR .873149. Conditional model path (root + child) tests mean 7.27; global FOUR tests exactly four per patient. V4 no robust benefit. [standalone hub mechanics CI #40568133585](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40568133585) PASS (previous CI setup missing source dependencies was repaired). **Therefore do not proclaim a victorious model tree**.

## NEW fixed benchmark question

Can a tiny, directly interpretable **GLOBAL sparse selective vote with four original within-patient rank-only directed inequality rules, of which ~1.75 typically fire, compete in positive F1 and Macro-F1 with REAL switchBox kTSP4 or even kTSP50, with much smaller model and smaller per-patient explanation?** Not inherently a new decision-tree publication result; the original unpublished RLR article is not delayed/edited. Candidate global selective voting method is an ablation/next hypothesis requiring real unseen experimental superiority before any paper claim.

Methods (all 11 matched same heldout subjects **per fold**):
1. TSP1 (original paper).
2. REAL switchBox `SWAP.KTSP` K4,10,20,50 **native predictions**, and TRUE **realized** pair count read from saved selected-pairs, possibly < nominal K.
3. Original RLR K10,20,50: original `d=2,m=10,alpha=.5,cv.glmnet lambda.min`, original supervised switchBox pair shortlist, all train only.
4. Standard value LR alpha=.5 with MAD-1000 supervised-unsupervised training prefilter matching original paper baseline.
5. **SELECTIVE4_M10**: four one-sided, class-directed rules, one-sided rank gap threshold fixed m10.
6. **SELECTIVE4_M10PCT**: threshold max(10,ceil(.10*(p-1))), p = probe count in original training assay gene universe; this variant was identified as more promising on previous D/V and is therefore *development-led*, NOT newly prospectively prespecified before any cohort was seen.

For the two new methods, same global class-weighted selection: choose up to 2 positive-class directed one-sided inequalities, up to 2 negative-class inequalities from top350 original switchBox candidate pool, gene-disjoint, training contrast magnitude x sqrt(binary rule prevalence variance), prevalence >=.08 and <=.92. For a new patient each oriented rule `h_j=1[sign_j*(rank(gene_right)-rank(gene_left))>=m]` contributes ONLY when h=1, else abstains. Patient vote sum maps to one of {-1,0,1}, converted to empirical-Bayes smoothed training positive-class posteriors with fixed strength=4. No absolute values, per-patient normalization to prescribed range, ratios, or learned scalar cutoffs.

**Complexity fields separately measured:** `nominal_panel` (four rules or true realized kTSP K); `mean_tested` (number of rules needing evaluation to determine votes; four for global selective method, K for native kTSP, number globally retained RLR coefficients for original RLR); `mean_active` (number of genuinely firing one-sided votes for global-four, K for kTSP, globally retained RLR coefficients with nonzero ternary rank feature for original RLR). `frac_active_le2` monitors how often ≤2 one-sided rules fire. **All four tests must still be evaluated to know whether they fire**. Do NOT claim lab measurement of only those genes; full expression array gene universe is still measured. For original RLR, both +1 and −1 rank states are locally active, only 0 is inactive.

## Repeated 470 outer E split audit

Original CompCancer protocol `9 cohorts x 10 repeats x5 folds + (West,Yeoh) x2 repeats x5 = 470 patient test folds per classifier`. [Firstfold end-to-end CI #40569478318](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40569478318) launched for all eleven original cohort names, real matched KTSP and RLR baselines. Each fold has original T-only feature/pair selection, exact same E subject IDs/labels, saved native kTSP labels, patient posterior, true active voting counts, train/sample ID lists, original RLR core commit hash and source SHA manifest. Separate patient-level independent Python audit recalculates F1/Macro-F1/BA/AUC/Brier/logloss and abstentions, with follow-up 470-level repeated test-patient coverage. If firstfold technical audit passes, the same immutable method code can be evaluated across 470; no outer E-driven hyperparam changes permitted.

**Source provenance limitation:** GitHub hosts publicly reconstructed Schliep expression matrices, not original paper-processed `CompCancer_truth_9x10`/ `CompCancer_truth_2x2` exact ARFF; checksum and probe-level parity unverified. Every result is exploratory/development/replication, NOT independent confirmation of exact historical paper inputs.

## Publication decision

- A global four-rule vote would have notable simplicity only if its F1/AUC on original 11×470 external splits and stability genuinely support a strong **complexity–discrimination tradeoff versus kTSP4 and kTSP50**, not just an anecdotal tiny patient example.
- If global4 wins but a true tree does not, publish it as sparse selective rank-rule panel or ablation; do not misleadingly call it a decision tree or imply a tree-specific causal benefit.
- If performance falls well below true kTSP50, leave original unpublished RLR ready for publication; no amount of rhetorical 'active rules' rescues unsupported classification superiority.
