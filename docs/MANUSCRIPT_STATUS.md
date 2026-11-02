# Current CMPB manuscript — single-author submission status

**Updated:** 2026-11-02  
**Manuscript:** *Relative Logistic Regression: A Probabilistic Refinement of Hard-Voting Pair-Based Omics Classification*  
**Manuscript author:** **Marcin Czajkowski** (sole author)  
**Affiliation:** Faculty of Computer Science, Bialystok University of Technology, Bialystok, Poland  
**Correspondence:** m.czajkowski@pb.edu.pl  
**Target journal:** *Computer Methods and Programs in Biomedicine* (CMPB)  
**Publication status:** prepared for submission; **no publication DOI, acceptance, or journal publication claimed**.

## Scope of the 2026-11-02 repository update

The revised MAIN PDF and Supplementary Material supplied to the research conversation both list **Marcin Czajkowski as the sole author**. This update aligns the public repository’s software citation metadata and landing page with that single-author manuscript.

The following are **not changed** by this authorship update:

- RLR implementation, pair nomination/selection and ternary encoding logic;
- the frozen CompCancer, PTB and legacy8 evaluation designs;
- score, statistical test, data-source provenance, reference outputs and SHA256 hashes;
- authorship of earlier **distinct** referenced multi-author publications, or the historical Git commit authors;
- MIT license, existing third-party package acknowledgements and historical repository records.

**The current journal-submission PDFs and LaTeX are not stored in this public reproducibility repository.** The repository intentionally contains code, data-acquisition documentation, frozen configurations, source/input checksums, reproduction scripts and reference outputs. The final PDF/LaTeX submission files are managed separately through the journal editorial process. In particular, the obsolete `manuscript_lock/` RLR3 snapshot was removed from the current main tree; older Git commits may still contain it.

For citation of this software, use the metadata in `CITATION.cff`. Because the paper has not been reported as published here, **do not invent a journal volume, article number, publication date or DOI**.

## Reproducibility

Run `python3 scripts/check_release.py --quick` for the repository source/config release checks and consult the main `README.md` for exact ARFF checksum preflight and Supplement S7.1 matching-pairs ablation. This metadata update does not assert that any computational experiment was rerun.
