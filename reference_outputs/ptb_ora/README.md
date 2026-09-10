# Corrected PTB fixed ORA

This directory contains the final training-only PTB fixed over-representation analysis retained in the paper.

- Predictive models were **not** refit.
- Gene lists were read from the frozen canonical PTB fixed model state.
- Background: `01_train_data.arff`, 103 samples, training-only MAD top 2000.
- Gene sets: official Broad MSigDB 2026.1.Hs Hallmark, GO Biological Process and Reactome GMT files bundled under `gmt/`.
- The published figures and tables are under `figures/` and `tables/`.

In `tables/gene_list_summary.csv`, `source_file` entries beginning with `release_reference/` are logical paths in the authors' frozen private working archive. They identify the exact audited model-artifact locations used to extract selected/active genes; the full frozen split-level tree is not duplicated in the public repository.
