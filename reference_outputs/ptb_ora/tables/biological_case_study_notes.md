# PTB fixed biological case study

 ## What was done
 - Main case study: PTB fixed-split.
 - Gene lists extracted for:
   - kTSP_10
   - kTSP_20
   - RLR_K10_m10_d2_pool1000
   - RLR_K20_m10_d2_pool1000
 - For kTSP, genes were taken from `selected_pairs.csv`.
 - For RLR, genes were taken from `active_pairs.csv` (or from non-zero coefficients if needed).
 - The background universe was **not** the whole genome. It was reconstructed from the PTB training set using the same training-only MAD prefilter (`top 2000 features`) that had been used before pair scoring.
 - Over-representation analysis (ORA) was performed separately for each model using Hallmark, GO Biological Process, and Reactome gene sets.
 - Only terms with FDR < 0.05 were retained.

 ## Intended interpretation
 This analysis should not be presented as validation of a biomarker panel.  The appropriate interpretation is that RLR does not necessarily introduce a completely different biological theme,  but may sharpen the signal toward a more coherent inflammatory / interferon / host-response core relative to classical kTSP shortlist genes.

 ## Key output files
 - `background_universe_mad2000.csv`
 - `gene_list_summary.csv`
 - `genes_kTSP_10.csv`
 - `genes_kTSP_20.csv`
 - `genes_RLR_K10_m10_d2_pool1000.csv`
 - `genes_RLR_K20_m10_d2_pool1000.csv`
 - `enrichment_full_all.csv`
 - `supp_table_top10_enrichment.csv`
 - `fig_ptb_fixed_ora_k20.png/pdf/svg`
 - `supp_fig_ptb_fixed_ora_k10.png/pdf/svg`
