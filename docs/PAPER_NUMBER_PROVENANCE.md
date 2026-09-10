# Paper-number provenance

`paper_numbers_manifest_v3.csv` and `.xlsx` are the numerical source map for the current manuscript lock.

Rules:
1. Primary values come from the canonical frozen/reference outputs.
2. CompCancer is one 11-dataset benchmark assembled from a 9-dataset x 10-repeat shard and a West/Yeoh x 2-repeat shard; inference is performed after dataset-level aggregation.
3. Candidate-pool and rank-margin screens are separate exploratory states and are not numerically interchangeable with the frozen primary benchmark.
4. PTB fixed metrics are from the predefined discovery-data split (103 training / 44 test), not from GSE83456.
5. PTB LOCO is descriptive.
6. Corrected PTB ORA uses the 103-sample training split, a training-only MAD-2000 universe, and official MSigDB 2026.1.Hs GMT files bundled under `reference_outputs/ptb_ora/gmt/`.
7. If prose/table values disagree with verified canonical outputs, correct the manuscript rather than altering the results.
