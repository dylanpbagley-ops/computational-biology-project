# Post-pilot reporting clarification for selected protocol v2

**Date:** 2026-10-03. **Timing:** written after the small pilot outputs were inspected and after reviewer feedback; before the broader 370-file retrieval and full-expression analysis. This addendum preserves `selected_protocol_v2.md` unchanged and does not alter the training observations, estimands, candidate thresholds, cutoff policies, code assignments, or primary 10% budget.

The full pilot's table-level exports retain the exact finite numeric threshold in `threshold_per_10k` and report `threshold_state` separately as `FINITE` or `ABOVE_MAX`. `ABOVE_MAX` means that the exact cutoff exceeds the largest training-reference value; it does not mean that every tumor cell is below threshold. The finite number remains available for independent reproduction.

The held-out code summary separates `heldout_healthy_codes_above_normal_budget` from `heldout_mm_codes_with_coverage_above_illustrative_fraction`. The same numerical fraction is a prespecified illustrative MM coverage summary, not a normal-reference violation. `NR` is retained when a metric does not apply or no data exist; it is not converted to zero.

This clarification follows the existing strict positive rule (`value >= threshold`) and numeric units of normalized marker mRNA molecules per 10,000 captured mRNA molecules. It adds no post-outcome exclusion, tuning, threshold choice, or clinical interpretation. The pilot output regeneration solely applied this serialization clarification to the already analyzed pilot data.
