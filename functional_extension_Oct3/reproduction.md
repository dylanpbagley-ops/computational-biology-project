# Reproduction and verification

This was a source-selection exercise only. There is no analysis pilot, random seed, new empirical dataset, or new software dependency. Python 3.12.14 and NumPy 2.3.5 were present; existing read-only h5py 3.16.0 was used only to inspect the prior HDF5 schema. No source archive was downloaded. The only attempted remote data-folder access was blocked by the browser access layer. Wagh's paper figure and supporting files provide the identified measured combination panel; Williamson's public CATNAP input concerns measured single-antibody susceptibility, while its combination values are calculated.

From `C:/Users/dylan.bagley/Desktop/outputs`, run:

```powershell
$env:PYTHONDONTWRITEBYTECODE = '1'
$taskPython = 'C:/Users/dylan.bagley/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe'
& $taskPython functional_extension_Oct3/validate_artifacts.py
& $taskPython longitudinal_extension_Oct3/preservation/verify_preservation.py
```

The stage validator checks finite JSON, consistent TSV widths, complete stage-hash coverage, unchanged directly inspected source hashes and every entry in both prior output manifests. The longitudinal preservation verifier reads its captured baseline and checks 63 protected inputs, the prior snapshot, 155 previously delivered outputs, directly used prior files and two documented status-suffix reconciliations. The separate stage preservation record verifies all 201 files in the longitudinal output manifest and all 155 files in the earlier delivery manifest. Expected exits are 0. No pilot command or test suite exists because no candidate passed its gates.
