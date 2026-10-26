# Blind Discovery of Interactions Between Within-Sample Gene-Ordering Relations: A Reproducible Feasibility Study

**Manuscript draft v0 — 26 October 2026. DEVELOPMENT-ONLY EVIDENCE. NOT SUBMISSION-READY.**  
Provisional title; target journal and author order not yet decided. This is a distinct experimental methods manuscript, not a revision of the locked paper *Relative Logistic Regression: A Probabilistic Refinement of Hard-Voting Pair-Based Omics Classification*.

> **Publication-status rule:** Every numerical claim below must remain linked to a frozen source commit, an auditable run, and a precise evaluation role. Synthetic D-only screen retrieval must NOT be described as held-out predictive performance, calibrated false-positive control, or clinical evidence. Sections flagged **REQUIRED BEFORE EMPIRICAL CLAIM** are intentionally unwritten.

## Abstract — v0 (provisional)

Classifiers based on within-sample relative expression orderings can yield concise, interpretable molecular signatures. Conventional candidate-ranking procedures, however, can miss interactions that carry no marginal class information at the level of a single pairwise comparison. We investigate whether interactions between two binary ordering relations can be retrieved without supplying either participating relation as an oracle.

We represent each within-sample ordering comparison as a binary signature across training samples. For XOR-like interactions, the observed outcome sequence is approximately the binary exclusive OR of two such signatures, potentially up to a complement. This observation enables outcome-aware Hamming hash joins that search promising disjoint relation pairs without enumerating every relation-pair combination. An implementation with 18-bit signatures and 24 hash tables was frozen before an independent synthetic robustness matrix. In an initial developmental pilot, the algorithm recovered both planted relations in the six highest-ranked interaction candidates in 19 of 20 synthetic datasets comprising 32–512 features. In a subsequent 288-study developmental stress test across three feature counts, three training sizes, and four interaction amplitudes, exact planted-interaction recovery was 70/72 at the strongest simulated amplitude (2.4), 36/72 at amplitude 1.6, and 4/72 at amplitude 1.0. All final shortlisted scores in the stress test were independently replayed from regenerated data.

These findings establish feasibility of blind rank-relation interaction *retrieval* under sufficiently strong synthetic effects, but also reveal pronounced sensitivity to effect size. The current system has not yet been evaluated as a complete independent-test biomedical classifier; additional evaluation with frozen parent-model offsets, separate discovery/confirmation/test data, null calibration, and independent cohorts is required.

**Keywords:** relative expression ordering; rank-based learning; interpretable classification; XOR interaction; relational features; locality-sensitive hashing; methodological feasibility.

## 1. Introduction

The order of two gene-expression measurements within one sample defines an interpretable relation with potential resistance to some monotone sample-specific scale transformations. Top-scoring pairs and related rank-based approaches exploit this property for molecular classification [R1,R2]. More expressive rules include pair-comparison decision trees and broader hierarchical rank structures [R3,R4]. As a complementary line of work, the separately locked Relative Logistic Regression (RLR) manuscript explores regularized probabilistic integration of selected ordering relations [R5].

An important limitation of marginal selection arises when the class signal is carried primarily by the *interaction* of two relations. Let a signed ordering indicator be

\[
z_{ab}(x)=\begin{cases}
+1,&x_a>x_b\\
-1,&x_a\le x_b.
\end{cases}
\]

A pure XOR-like model may have \(\Pr(Y\mid z_{ab})\approx\Pr(Y)\) and \(\Pr(Y\mid z_{cd})\approx\Pr(Y)\), while the product \(z_{ab}z_{cd}\) is predictive. Screening individual ordering relations by their marginal class association cannot be assumed to retrieve either member of such a structure.

This paper studies a specific computational question: **can an unknown interacting pair of rank-ordering relations be found from labeled discovery samples without giving the algorithm either true relation?** The proposed hash-join retrieval is deliberately separated from downstream conditional-tree fitting and from published RLR model evaluation. Our evidence is primarily synthetic and mechanism-oriented; the paper must not conflate identifying planted rules on discovery data with improving held-out biomedical predictions.

### Intended contributions — pending final validation

1. A clear formulation of **joint** rank-relation discovery when the individual members of the interaction have weak or zero marginal predictive value.
2. A binary signature hash-join candidate generator for disjoint order relations, with orientation-invariant handling of the two true relations; *no oracle anchor is used in blind retrieval*.
3. Reproducible experimental separation of candidate-pool coverage, interaction shortlist ranking, and oracle partner-completion capabilities, including revealing failure cases.
4. A pre-frozen robustness grid quantifying the sharp dependence on synthetic effect size, sample count and feature dimension.
5. **Required for broader empirical claims:** a strictly held-out original-RLR comparison, statistically calibrated null behavior, and biological cohort evaluation.

## 2. Related work and distinction from existing RLR

Original pair-ranking methods include TSP and k-TSP [R1,R2]. Pair-ordering decision trees, including the earlier work by Czajkowski and Krętowski [R3], predate the present study; therefore **the use of gene-pair comparisons in a tree is not itself novel**. Surveys of gene-pair biomedical methods [R4] and a broader taxonomy of biomedical rank-based methods [R6] establish existing alternatives and emphasize variation in model, search and scalability.

The archived RLR manuscript [R5] addresses a different question: improving the *probabilistic combination* of already selected relative-ordering predictors. Our present question is **interaction-aware discovery of relational predictors**, including pure XOR cases where ordinary marginal pair ranking fails.

We do not yet claim that our specific LSH implementation is unique among all XOR, Boolean-function discovery, feature-interaction or approximate-nearest-neighbor algorithms; a dedicated comparison to general-purpose interaction search and binary LSH is needed before a final novelty statement. Nor do we assert practical transcriptome-wide scalability: creating signatures for all unordered gene pairs already requires \(O(n p^2)\) memory or streamable computation.

## 3. Methods

### 3.1 Signed within-sample order relations

For each sample \(i\) and unordered feature pair \((a,b)\), define the binary ordering signature \(S_{ab,i}=\mathbb{1}[X_{ia}>X_{ib}]\), with a fixed canonical orientation \(a<b\). Each pair induces a binary vector over the discovery set \(D\) of size \(n_D\). Two disjoint feature pairs \((a,b)\) and \((c,d)\) generate an interaction by XOR of their vectors, equivalent to a signed product under the \(\{-1,+1\}\) encoding.

A planted scenario used in simulation is
\[
\operatorname{logit}\Pr(Y_i=1\mid X_i)
=\beta\, z_{12}(X_i) z_{34}(X_i),
\]
where the four gene measurements are sampled independently from standard-normal distributions. The null uses \(\beta=0\). This generator supplies a mechanism test, not a biologically realistic population model.

### 3.2 Discovery-only blind Hamming hash joins

Let \(M=p(p-1)/2\) be the number of canonical rank comparisons. H3 constructs \(M\) binary signatures of length \(n_D\). In each of 24 independently seeded hash tables, 18 training individuals are selected as bit positions. Hash codes are stored for all \(M\) signatures. For each candidate relation \(S_a\), bins are queried using \(h(S_a)\oplus h(Y)\) and its 18-bit complement, seeking \((S_a,S_b)\) whose XOR agrees with \(Y\) or its complement at those coordinates. Candidate pairs sharing a gene are excluded.

The entire discovery procedure uses only the \(D\) features and labels. Candidate interaction pairs returned by the hash joins are evaluated using the absolute signed agreement
\[
  T(a,b)=
  \left|\,\frac{2}{n_D}\sum_i
      \mathbf{1}\!\left[(S_{a,i}\oplus S_{b,i})\ne Y_i\right]-1
  \right|.
\]
The six highest-scoring **unique unordered** disjoint relation-pair combinations are retained. The known planted true relation IDs are inspected **only in retrospective recovery measurement**, not in hash construction, proposal generation, scoring or selection.

All hash counts, source seeds, actual candidate IDs, scores and source code are archived. Selection and scoring on the same labeled \(D\) may overfit; \(T\) is not a held-out accuracy estimate or a p-value.

### 3.3 Experimental roles

**Stage F:** frozen original CRIT-v9 versus *oracle-given* true relations in 12 reused synthetic null/XOR studies. Tests whether the model form and V acceptance could exploit the planted signal if discovery succeeded. **No blind retrieval claim.**

**Stage G/H0:** exhaustive all 496 unordered rank pairs at 32 genes, followed by otherwise unchanged CRIT-v9 D selection and V gates; 24 reused and 24 new synthetic development datasets, each with an independent E set. Corrected orientation-invariant recovery in XOR was 12/12 (5 direct, 7 gate/partner swapped). **Results are developmental and designed after observing earlier failures.**

**Stage H1:** isolated, simplified D-only screening with a 320-pair budget and a full 496-pair control at p32; fixed random, per-study random and graph-based proposal strategies at p32–512. These are **not matched complete-RLR classifier performance comparisons**.

**Stage H2:** fixed synthetic data with one planted true oracle relation supplied, exhaustive D-only partner search over gene-disjoint pairs, p32–1024. Establishes the conditional partner-search ceiling, **not feasible blind prediction**.

**Stage H3-fix1:** blind 18-bit/24-table Hamming hash-join prototype with four D-only replicates for each p32,64,128,256,512 under strong \(\beta=2.4\) and null.

**Stage H4A (primary method robustness evidence):** frozen unchanged H3 estimator hash across an independent 3×3×4 grid: \(n_D\in\{60,120,240\}\), \(p\in\{32,128,512\}\), \(\beta\in\{0,1.0,1.6,2.4\}\), eight independent seeds per cell (288 discovery sets). Outcome is exact planted pair-product present in the top-six, not predictive test logloss. Regenerated data and candidate scores are audited by a separate implementation without calling the estimator.

### 3.4 Preplanned missing tests — REQUIRED BEFORE CLASSIFIER CLAIMS

The blind LSH stage must be integrated in an isolated full pipeline with exact original published RLR parent predictions and strict within-D out-of-fold residuals. Model nomination must use **D labels only**, select a single candidate via its own frozen D-only rule, pass through the existing one-shot V confirmation, and be evaluated on a completely untouched E test set. The primary endpoint will be \(\Delta\text{logloss}_{E}\) relative to frozen ORIGINAL_RLR, alongside Brier score, acceptance, fallback parity and full provenance.

Formal correction acceptance / familywise null calibration is not measured by an arbitrary reference pair's absence from a null shortlist. Before claiming controlled false acceptance, rerun the **entire selection and confirmation pipeline** over a sufficiently large set of independently generated null datasets, publish the binomial uncertainty intervals, and specify whether each method arm or the five-arm family is the target.

Any independent biomedical comparison must use cohort(s) not repeatedly consulted while developing CRIT-v7 through H4A. Predeclare exact datasets, feature normalization, patient-independent partitions, baselines and primary endpoint.

## 4. Results — verified frozen development evidence

### 4.1 Why CRIT-v9 missed planted XOR and the oracle could fit it

The Stage F frozen candidate-pool audit found the planted gene01>gene02 gate absent in **6/6** XOR replicates, while the planted gene03>gene04 partner relation was present in **6/6**. Discovered CRIT-v9 corrections were accepted 0/6. When given the known true structure, an oracle conditional tree passed V 6/6 and yielded mean heldout E logloss change **−0.2167817** relative to original RLR. The oracle global-interaction control likewise passed 6/6 and yielded **−0.2877593**. These are oracle interventions on *exposed synthetic development*, not the accuracy of blind H3.

A complete 32-gene pool in Stage G made all 496 rank pairs available; the existing CRIT-v9 nomination/ranking/confirmation then found a structurally correct true pair-product, considering gate–partner swap symmetry, in **12/12** XOR development experiments. H0 independently checked 5 direct plus 7 swapped structures; all 12 conditional trees were accepted on V and improved E logloss. The Stage G mean conditional-tree E logloss changes were −0.21678631 (six reused Datasets) and −0.21590897 (six new development Datasets). These results support the computational explanation that missing candidate coverage, not inability to represent the interaction, caused the v9 XOR failure.

### 4.2 The first unknown relation is the principal discovery challenge in this generator

In H1, for p32 the graph-based 320-pair proposal had exact planted pair-product in the simple D top-six screen in **7/8** XOR datasets, but reached only **2/8** at p64, **1/8** at p128, **0/8** at p256 and **1/8** at p512. Fixed-random and per-study-random methods generally had even lower recall. The different synthetic seeds and designs across stages must not be treated as a matched statistical head-to-head comparison.

By contrast, Stage H2 supplied one true rank pair as an oracle and exhaustively ranked all disjoint partner pairs. The true partner ranked **first in 16/16 oracle orientations** at every investigated dimension from p32 through p1024 (8 replicates for each orientation). A sham starting pair did not retrieve the reference target in top-six in any of 8 replicates at those dimensions. This supports a distinction between blind initial-pair discovery and partner completion **under the synthetic generator with an oracle anchor**.

### 4.3 Blind XOR retrieval with no planted anchor

H3-fix1 recovered the precise planted unordered pair-product among the top six in **19/20** discovery-only XOR datasets across p32–512 at \(n_D=120,\beta=2.4\). By dimension, the counts were **4/4** at p32, p64, p128 and p256, and **3/4** at p512. The missed p512 case had no true-pair hash collision in the 18-bit/24-table indexing. Independent audit replayed all **240** finalist scores on regenerated D data. No original RLR prediction, confirmation V or test E was involved.

### 4.4 Robustness of the frozen method to weaker interaction strength

All H4A simulations used a **byte-identical** H3 estimation source and new seeds. The primary observational counts below pool over three nD sizes, **24 simulations per (p,amplitude)**; pooled counts are descriptive, not homogeneous clinical response-rate estimates.

| Synthetic interaction coefficient \(\beta\) | p=32 | p=128 | p=512 | All p and nD |
|---|---:|---:|---:|---:|
| 0.0 (null; arbitrary reference relations) | 0/24 | 0/24 | 0/24 | 0/72 |
| 1.0 (weak) | 2/24 | 2/24 | 0/24 | **4/72** |
| 1.6 (moderate) | 11/24 | 16/24 | 9/24 | **36/72** |
| 2.4 (strong) | 24/24 | 23/24 | 23/24 | **70/72** |

In the p512, nD60, beta1.6 cell, the true interaction reached a hash collision in 5/8 replicates but entered top-six in only 1/8. Thus the limitation concerns **both retrieval and final competition among scored candidates**. In a p512, nD120, beta1.0 cell, collision occurred in 1/8 but top-six was 0/8.

The independent auditor regenerated all 288 H4A datasets and recomputed all **1728** finalist scores from saved candidate gene IDs without importing the H3 search function. The CI runner passed unit tests, source parity and all archive gates. The 0/72 arbitrary-reference null entry is **not** an estimate of the full classifier's false-positive or false-correction rate. It is also not evidence that the highest-scoring null interactions are statistically innocuous.

### 4.5 Computational costs and limits

H3 represents each of \(M=p(p-1)/2\) gene pairs as an \(n_D\)-bit signature. The pair-level representation is already \(O(n_D p^2)\); while hash joins avoid straightforward exhaustive \(O(M^2)=O(p^4)\) scoring, fan-out of collisions can be substantial. For p512, \(M=130\,816\) rank-pair signatures; in the H3 pilot, roughly 1.55 million table-level candidate interactions could be scored per discovery dataset (including duplicates across independent hash tables). The H4A measured mean elapsed times were approximately 0.007 s at p32, 0.09–0.10 s at p128 and 2.3–2.6 s at p512 on GitHub-hosted CPUs under the tested nD settings. These are implementation- and platform-specific timings and do **not** include complete original-RLR fitting and testing.

## 5. Discussion

The experimentally important distinction is between **representation**, **search**, and **statistical validation**. The synthetic XOR interaction is easy for an oracle-supplied relation product and can also be fitted by a conditional rank tree when all underlying rank pairs are exposed to the discovery mechanism. Blind discovery is harder, because neither constituent pair is marginally informative. Binary signature hash joins offer an outcome-aware way to propose joint interactions from labeled D without needing an oracle pair.

A frozen pilot shows encouraging recovery of high-amplitude pure XOR even with several hundred genes. Equally important, a larger predeclared stress test demonstrates substantial performance collapse for weak effects and identifies separate hash-collision and post-collision ranking limitations. Describing the algorithm simply as "solving high-dimensional XOR" would be inaccurate.

Further limitations include the synthetic Gaussian feature generator, perfect canonical rank measurements, no biological correlation/background structure, feature counts far smaller than routine whole-transcriptome assays, label-adaptive screening on the same D data, lack of heldout predictive comparison for blind H3/H4A, limited per-cell replication, and lack of multiplicity-calibrated full-pipeline error guarantees. These caveats constrain any claim of biomedical benefit.

### Decision after current evidence

The appropriate next step is **not** to continue retuning H3 against these already inspected datasets. It is to freeze a model integration and conduct a limited, well-audited independent heldout D/V/E evaluation against unchanged RLR. A broader journal performance claim additionally requires fresh biological cohort evidence and a null-calibration study. If such evidence is unfavorable, the current material supports at most a narrower, method-focused feasibility article explicitly highlighting weak-signal limitations.

## 6. Conclusion — provisional

We documented a reproducible blind hash-based approach to retrieving interactions of interpretable within-sample ordering relations. In controlled synthetic experiments it reliably recovered strong XOR structures without a true pair supplied as an oracle, while detection of weaker effects remained unreliable. The proposed search is a feasibility result with clear limitations, not yet a validated or statistically calibrated biomedical classifier.

## Evidence ledger: frozen sources, runs and claim scope

| Source / experiment | Frozen evidence | Legitimate claim |
|---|---|---|
| Locked parent RLR publication archive | \`manuscript_lock/RLR3_submission_final.tex\`; main README; canonical original code | Distinct previous probabilistic rank-pair method, NOT a new H3 outcome |
| Stage F oracle / missing gate | Stage F SHA \`94d2bc049bf987d08d87637847427c00cf1acee9\`; audited local report | Oracle mechanism and v9 missing-pool diagnosis only |
| Stage G / H0 exhaustive, 32 genes | \`6ca44c30a439dbdb97b3a45f7e8f45432da3bd59\`; H0 SHA \`21986f8ece37ae5056ea1123360e12c74bcf60d3\`; archived PASS | Conditional rule XOR structure recovery 12/12 in developmental simulation |
| Stage H1-fix1 bounded proposals | [GitHub Actions PASS 40285624939](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40285624939) | D-only proposal coverage / ranking |
| Stage H2 known-anchor oracle | [GitHub Actions PASS 40286700787](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40286700787) | Partner completion ceiling with supplied true anchor |
| Stage H3-fix1 blind LSH | SHA \`cd1b3038f3b662216a19655faccbec394de5900e\`; [PASS 40288356046](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40288356046) | Blind D-only planted XOR retrieval 19/20 at strong effect |
| Stage H4A frozen robustness | SHA \`14f6fe34a3bab4dd90a3705f66dbb886eb101451\`; [PASS 40289616641](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40289616641) | 288 new D-only simulations; robust-score replay and effect-size dependence |
| **REQUIRED BEFORE JOURNAL EMPIRICAL CLAIM** | **NOT YET AVAILABLE** | Untouched heldout E prediction, true original RLR parity, null alpha calibration, new biological cohort inference |

## References to verify and finalize

- [R1] Geman D. et al. (2004). *Classifying gene expression profiles from pairwise mRNA comparisons.* (Find exact bibliographic details/DOI at reference-lock.)
- [R2] Tan A.C. et al. (2005). k-Top Scoring Pairs classification (verify exact title/journal/DOI).
- [R3] Czajkowski M., Krętowski M. (2011). *Top scoring pair decision tree for gene expression data analysis.* Advances in Experimental Medicine and Biology 696:27–35. DOI 10.1007/978-1-4419-7046-6_3. https://pubmed.ncbi.nlm.nih.gov/21431543/
- [R4] *Applications of gene pair methods in clinical research: advancing precision medicine.* Molecular Biomedicine (2025). https://doi.org/10.1186/s43556-025-00263-w (verify full author details).
- [R5] Czajkowski M. et al. *Relative Logistic Regression: A Probabilistic Refinement of Hard-Voting Pair-Based Omics Classification.* Archived manuscript in this repository; publication status not assumed.
- [R6] Czajkowski M., Jurczuk K., Krętowski M. (2026). *Rank-based relational methods for interpretable biomedical modeling: a taxonomic review.* Artificial Intelligence Review. DOI 10.1007/s10462-026-11679-3.

## Final submission placeholders (do not fabricate)

- [ ] Full frozen RLR+LSH heldout protocol executed; inspect exact E logloss/Brier gains vs baseline.
- [ ] Primary inferential endpoint, repeated-test uncertainty and null calibration specified.
- [ ] Independent unexposed biological cohort(s), provenance and auditable patient-level predictions.
- [ ] Empirical baselines and novelty comparison beyond classical pair trees and existing feature interaction algorithms.
- [ ] Citation details/permission/author contribution approval and appropriate journal formatting.
- [ ] Explicit conclusion revisited after all studies, including negative results.
