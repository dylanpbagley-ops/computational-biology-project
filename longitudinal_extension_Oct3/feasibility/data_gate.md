# G2 — paired-data feasibility decision

## Identity, order, and exposure

The primary Table 1 in Leblay et al. (Blood 2024, DOI 10.1182/blood.2023020276; https://pmc.ncbi.nlm.nih.gov/articles/PMC10797556/, indexed full text lines 225–236) separates scRNA/scATAC specimen columns from scDNA specimen columns and is headed Prevenetoclax/Postvenetoclax. It links these RNA samples:

| Participant | Pre-venetoclax RNA | Post-venetoclax RNA | GEO records |
|---|---|---|---|
| PM01 | P1393 | P1555 | GSM5971019, GSM5971023 |
| PM07 | P1389 | P1479 | GSM5971021, GSM5971025 |
| PM10 | P1217 | P1438 | GSM5971020, GSM5971024 |

PM01's scDNA pre specimen is P1426, not its scRNA pre specimen P1393. It is not substituted. PM09's table lists P1498/P1698, but the selected GEO family record contains no P1698 RNA sample; PM09 is excluded rather than guessed.

These are independently collected specimens with a documented participant key and order relative to venetoclax. They are not lineage-tracked cells. Exact collection dates and intervals are NR. The GEO titles say pre/post treatment generically; the paper's table supplies the venetoclax-specific mapping.

## Public measurement records

A 5,344-byte public GSE199359 family SOFT file was saved under data/source_metadata with compressed and uncompressed SHA-256 values in its manifest. Its six sample records provide exact supplementary paths and say Primary CD138+ MM cells; the associated file description is a CellRanger filtered feature-barcode HDF5 matrix. All six public paths returned HTTP 200 to HEAD and GET. The downloaded files total 100,468,485 bytes (95.8 MiB); the largest is 22,301,164 bytes (21.3 MiB), below both caps. Exact bytes, hashes and timestamps are in data/download_manifest.json/.tsv. Exact paths, sizes, and response metadata are in GSE199359_pair_samples.tsv and GSE199359_h5_file_probes.tsv.

All six individual GEO sample protocols identify 10x reagent-kit CG00052 v2. The series-level overall-design text says 3-prime v3. This metadata discrepancy is preserved; the six individual sample records agree with each other, but exact chemistry remains uncertain pending file/header and source supplement review. The paper states scRNA-seq reads were processed with CellRanger v3.0.2 against GRCh38. The supplementary methods were not accessible in this bounded pass.

## Cell-label boundary

The source identifies a CD138-positive selected fraction and reports cell-state analyses, but the GEO sample records do not supply per-barcode malignant/clone annotations. The downloaded filtered matrices are described as CD138-selected RNA barcodes/cells only. Do not call every sequenced cell malignant, infer labels from TNFRSF17/GPRC5D/SLAMF7, or claim target loss is protein loss or functional escape. No paired surface-protein assay for these six records was found.

## Decision

**G2: pass for the scope-reduced, descriptive RNA pilot.** All six HDF5 matrices validated as integer Gene Expression-only count matrices; each contains one exact feature for each predeclared target. The 100-replicate joint depth sensitivity completed at D=1,318. The malignant-state gate does not pass; no full analysis, independent validation, or generalization claim is supported by three participants. The later strong-contribution gate failed because participant-level changes were heterogeneous, PM10 post depth was markedly lower, and the scope remains three participants/CD138-selected RNA with no per-barcode labels or matched protein.

