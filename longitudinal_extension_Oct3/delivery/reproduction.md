# Reproduction and independent review

This stage is self-contained under C:\Users\dylan.bagley\Desktop\outputs\longitudinal_extension_Oct3. Run the commands from PowerShell. They use the bundled Python runtime and stage-local HDF5 reader; no package installation is required to reproduce the analysis. The plotting script reads the pre-existing Matplotlib installation under strengthening_execution_Oct2\dependencies and writes only into this stage.

~~~powershell
$root = 'C:\Users\dylan.bagley\Desktop\outputs'
$stage = Join-Path $root 'longitudinal_extension_Oct3'
$python = 'C:\Users\dylan.bagley\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe'
Push-Location $root
& $python (Join-Path $stage 'analysis\validate_h5_structure.py')
& $python -m unittest discover -s (Join-Path $stage 'analysis\tests') -v
& $python (Join-Path $stage 'analysis\run_pilot.py')
& $python (Join-Path $stage 'analysis\plot_pilot_figure.py')
& $python (Join-Path $stage 'preservation\verify_preservation.py')
Pop-Location
~~~

The protocol lock is analysis\pilot_protocol.lock.json; it records SHA-256 f32b335d2bc88a1571f4045452033d422049304b501738c032c966bb8c45a9ac and predates all matrix downloads and expression parsing. The six public H5 inputs are listed with retrieval status, byte counts, and SHA-256 values in data\download_manifest.tsv and .json. The download total is 100,468,485 bytes. The pilot script checks the saved protocol hash before analysis and regenerates the numeric result tables and run manifest.

The fixed pilot is descriptive pre/post-venetoclax RNA detection for PM01, PM07, and PM10 and targets TNFRSF17, GPRC5D, and SLAMF7. The independent unit is participant (three total); individual barcodes are nested subsamples. The depth sensitivity jointly thins target and remaining RNA counts with a multivariate hypergeometric draw at D=1,318, using 100 predeclared seeds. These replicates describe sampling sensitivity and are not participant-level confidence intervals.

Numeric source tables are in analysis\results. Figures are paired_joint_category_change.svg/.png and paired_double_undetected_raw_vs_depth_sensitivity.svg/.png in that directory. The first figure shows deposited raw-count detection (>0), not depth-adjusted estimates. The second shows raw versus mean joint-thinning sensitivity values.

Independent checks recorded in logs\assistance_log.tsv include: preservation verification of 301 file references with two pre-existing status-suffix reconciliations and zero failures; parent standalone reruns of the six-file structure validator, eight focused tests, final pilot, and figure generator; direct raw-H5 recomputation of all 18 primary joint profiles and all 44,981 barcode rows; and independent checks of the hypergeometric summaries against exact probabilities and 9,000 replicate-delta records. PNG dimensions and approximately 300 dpi metadata were inspected.

The contribution decision is in delivery\path_decision.md, not the feasibility folder: the final decision was grouped with delivery materials although the plan had placed the decision under feasibility. The matching STATUS.md and decision log point to the actual location. No full extension analysis is authorized by the failed contribution gate. The only remaining options are metadata-only tests in delivery\adjacent_gap_data_cards.md; each requires a new novelty-overlap check before any further data analysis.