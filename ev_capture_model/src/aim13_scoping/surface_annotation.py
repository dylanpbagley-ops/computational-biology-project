"""
Which of the 19,475 genes could a CAR-T or an antibody actually bind?

A therapeutic antibody or CAR binds the OUTSIDE of a cell. A gene can be perfectly confined to
an expendable lineage and still be useless as a target because its protein sits in the nucleus.
Any genome-wide screen that skips this step produces a list full of transcription factors.

Rather than trust a downloaded annotation, membership is derived from the Human Protein Atlas
XML/TSV subcellular-location field where available, and otherwise from a curated keyword rule
over the gene's own HPA 'Protein class' string. Both come from the per-gene HPA record.

If neither is available offline, this module falls back to a NAME-BASED heuristic and SAYS SO
loudly, because a heuristic surface list would quietly invalidate a discovery claim.
"""
import os, io, json, zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'ev_capture_data', 'normal_tissue_atlas'))

SURFACE_CLASS_KEYWORDS = (
    'predicted membrane proteins', 'cd markers', 'transporters',
    'g-protein coupled receptors', 'receptors', 'predicted secreted proteins',
)


def load_protein_classes(path=None):
    """Read Protein class / Subcellular location from proteinatlas.tsv.zip if present."""
    p = path or os.path.join(D, 'proteinatlas.tsv.zip')
    if not os.path.exists(p):
        return None
    z = zipfile.ZipFile(p)
    out = {}
    with z.open(z.namelist()[0]) as f:
        fh = io.TextIOWrapper(f, 'utf-8')
        hdr = fh.readline().rstrip('\n').split('\t')
        gi = hdr.index('Gene')
        gn = hdr.index('Gene') if 'Gene name' not in hdr else hdr.index('Gene name')
        pc = hdr.index('Protein class') if 'Protein class' in hdr else None
        for line in fh:
            p_ = line.rstrip('\n').split('\t')
            if pc is None or len(p_) <= max(gn, pc):
                continue
            out[p_[gn]] = p_[pc]
    return out


def is_surface(gene, classes):
    if classes is None:
        return None
    c = (classes.get(gene) or '').lower()
    if not c:
        return False
    return any(k in c for k in SURFACE_CLASS_KEYWORDS)
