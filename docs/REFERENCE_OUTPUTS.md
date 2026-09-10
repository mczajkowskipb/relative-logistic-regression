# Reference outputs

The stage-6 collection verified the frozen result-tree sizes used as the numerical source of truth:

- CompCancer: 4,230 split-level metric files and 4,230 test-prediction files;
- PTB fixed: 9 + 9;
- PTB LOCO: 36 + 36;
- legacy8: 3,600 + 3,600.

The public publication repository intentionally does **not** duplicate the complete split-level frozen result trees. Instead it contains the full publication code, portable canonical configurations, environment capture, input/reference checksums, data-acquisition documentation, the paper-number provenance manifest, locked manuscript assets, and the corrected PTB ORA outputs. The complete frozen working archive is retained privately by the authors. With checksum-matched processed inputs, the publication states can be regenerated from the public scripts and configurations.

The corrected PTB ORA is bundled completely under `reference_outputs/ptb_ora/`, including the exact MSigDB 2026.1.Hs GMT inputs, tables and PNG/PDF/SVG figures used by the corrected analysis.

`paper_numbers_manifest_v3.csv` / `.xlsx` maps each locked paper number to its canonical source path or corrected ORA output. Paths beginning with `release_reference/` are logical provenance paths within the authors' frozen internal archive; they are retained to identify the audited source state and are not expected to exist as directories in a fresh public clone.
