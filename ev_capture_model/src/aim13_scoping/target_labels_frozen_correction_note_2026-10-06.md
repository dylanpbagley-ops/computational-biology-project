# Correction note, October 6, 2026

The file `target_labels_frozen.tsv` previously contained 40 of the 44 operational
targets; AXL, ERBB3, F3, and GPNMB were missing (F3 is one of the 11 toxic targets in
the scored panel). The TSV has been regenerated from `target_labels.py`, which is the
source of truth; all 44 rows are now present. **No outcome label changed**: for the 40
previously present targets, every outcome, agent, and expected-normal-site field is
identical to the prior version. The single exception is the SLC34A2 (NaPi2b) evidence
citation, which was expanded to match `target_labels.py` (it adds the September 2026
documentation note on upifitamab rilsodotin's independent discontinuation; the label
itself is unchanged). Evidence fields are documentation-only and never enter any
analysis. New file SHA-256:
`148bbd3b7314bec9f79a0925219b110302dc9eb0b7e2f38b05d3833ae80ba9f6`. The four
precommitment hash records that pinned the 40-row TSV are superseded as of this date;
the original file and hashes are preserved under `strengthening_Oct2/before/`. No
analysis result is affected: all analyses read `target_labels.py`, not the TSV.
