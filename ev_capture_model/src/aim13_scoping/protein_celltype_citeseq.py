"""
DOES THE CONFINEMENT MEASURE TRANSFER FROM RNA TO PROTEIN, AT CELL-TYPE RESOLUTION?

This is the question results section 26 left open and section 31 could not close.

Section 26 established: the safe-vs-toxic result holds when restriction is measured across
normal CELL TYPES (RNA 0.868, antibody staining 0.785) and fails at whole-ORGAN resolution in
BOTH technologies (RNA 0.626, mass-spectrometry protein 0.650). Section 31 tried to settle it
with HPA's Deep Visual Proteomics and could not: 7 of the 8 targets that reference fails to
detect are the safe blood-lineage antigens the framework rests on.

THE PROJECT ALREADY OWNED A PROTEIN, CELL-TYPE-RESOLUTION REFERENCE and had not used it for
this. The healthy-donor CITE-seq panel from section 16: 161,764 cells, 8 donors, **228 surface
proteins measured by ANTIBODY** -- the same physical modality a CAR or a therapeutic antibody
uses to bind -- annotated into 31 blood cell types (B naive/memory/intermediate, plasmablasts,
every T subset, NK, monocytes, dendritic cells, erythrocytes, platelets, HSPC).

WHAT IT CAN AND CANNOT TEST.
  CANNOT: the safe-vs-toxic outcome test itself. The panel is blood only and the toxic class is
  almost entirely solid-tissue targets, leaving too few toxic cases to test.
  CAN: the question underneath it. If protein-measured confinement agrees with RNA-measured
  confinement over the SAME cell types, then the RNA-based result is a statement about the
  protein too, and section 26's scope condition is about RESOLUTION rather than about RNA. If
  they disagree, the RNA result does not transfer, and that is a genuine limitation.

The comparison is like-for-like: same cell types, same measure, one antibody-based and one
transcript-based.
"""
import os, io, json, zipfile, collections
import numpy as np
from scipy.stats import spearmanr
import target_labels as TL
from adt_gene_map import MAP

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'ev_capture_data'))
ADT = os.path.join(DATA, 'healthy_pbmc_citeseq', 'healthy_adt.npz')
ATLAS = os.path.join(DATA, 'normal_tissue_atlas', 'rna_single_cell_type.tsv.zip')

BLOOD_RNA = ('b-cells', 'plasma cells', 't-cells', 'nk-cells', 'monocytes', 'macrophages',
             'granulocytes', 'neutrophils', 'eosinophils', 'basophils', 'cdc', 'pdcs',
             'erythroid cells', 'erythrocytes', 'megakaryocytes', 'langerhans cells',
             'kupffer cells', 'microglia', 'mast cells', 'hofbauer cells',
             'erythrocyte progenitors', 'neutrophil progenitors', 'thymocytes',
             'megakaryocyte-erythroid progenitor cells')


def confinement(v):
    v = np.clip(np.asarray(v, float), 0, None)
    tot = v.sum()
    if tot <= 0:
        return None, None
    top = float(v.max())
    return float(top / tot), int((v > 0.10 * top).sum())


def main():
    d = np.load(ADT, allow_pickle=True)
    X, names = d['adt'], [str(x) for x in d['names']]
    ct = np.array([str(x) for x in d['celltype_l2']])
    keep = ct != 'Doublet'
    X, ct = X[keep], ct[keep]
    types = sorted(set(ct))
    print('protein reference: %d cells, %d proteins, %d blood cell types, %d donors'
          % (X.shape[0], X.shape[1], len(types), len(set(str(x) for x in d['donor']))))

    means = np.vstack([X[ct == t].mean(0) for t in types])
    prot = {}
    for j, nm in enumerate(names):
        g = MAP.get(nm)
        if not g or g in prot:
            continue
        s, n = confinement(means[:, j])
        if s is not None:
            prot[g] = {'top_share': s, 'n_types': n, 'home': types[int(np.argmax(means[:, j]))]}
    print('ADT markers mapped to unique genes: %d' % len(prot))

    z = zipfile.ZipFile(ATLAS)
    vals = collections.defaultdict(dict)
    want = set(prot)
    with z.open(z.namelist()[0]) as f:
        fh = io.TextIOWrapper(f, 'utf-8')
        hdr = fh.readline().rstrip('\n').split('\t')
        gi, ci, vi = hdr.index('Gene name'), hdr.index('Cell type'), hdr.index('nCPM')
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if p[gi] in want and p[ci].lower() in BLOOD_RNA:
                try:
                    vals[p[gi]][p[ci]] = float(p[vi])
                except ValueError:
                    pass
    rna = {}
    for g, dd in vals.items():
        ks = list(dd)
        s, n = confinement([dd[k] for k in ks])
        if s is not None:
            rna[g] = {'top_share': s, 'n_types': n,
                      'home': ks[int(np.argmax([dd[k] for k in ks]))]}
    nblood = len(next(iter(vals.values()))) if vals else 0
    print('RNA confinement computed over the same %d blood lineages for %d genes'
          % (nblood, len(rna)))

    shared = sorted(set(prot) & set(rna))
    xa = np.array([prot[g]['top_share'] for g in shared])
    xb = np.array([rna[g]['top_share'] for g in shared])
    rho, p = spearmanr(xa, xb)
    rho_n, p_n = spearmanr([prot[g]['n_types'] for g in shared],
                           [rna[g]['n_types'] for g in shared])
    same_home = sum(1 for g in shared
                    if prot[g]['home'].split()[0].lower()[:2] == rna[g]['home'].lower()[:2])

    print('\n' + '=' * 86)
    print('DOES CONFINEMENT TRANSFER FROM RNA TO PROTEIN AT CELL-TYPE RESOLUTION?')
    print('=' * 86)
    print('  genes measured both ways: %d' % len(shared))
    print('  confinement (top_share):  Spearman rho = %.3f   p = %.2e' % (rho, p))
    print('  number of cell types:     Spearman rho = %.3f   p = %.2e' % (rho_n, p_n))

    print('\n  FOR SCALE, the agreements already reported elsewhere in this project:')
    print('    RNA cell-type vs RNA bulk-tissue          rho = 0.512  (section 23.7)')
    print('    RNA cell-type vs PROTEIN mass spec        rho = 0.346  (section 26.3)')
    print('    RNA cell-type vs PROTEIN antibody stain   rho = 0.429  (section 26.3)')
    print('    RNA blood vs PROTEIN antibody, SAME cell types -> the number above')

    LAB = {t[0]: t[3] for t in TL.TARGETS}
    tg = [g for g in shared if LAB.get(g) in ('safe', 'toxic')]
    if tg:
        print('\n  therapy targets measured both ways:')
        print('  %-10s %-7s %8s %8s  %-20s %s'
              % ('GENE', 'outcome', 'protein', 'RNA', 'protein home', 'RNA home'))
        for g in tg:
            print('  %-10s %-7s %8.3f %8.3f  %-20s %s'
                  % (g, LAB[g], prot[g]['top_share'], rna[g]['top_share'],
                     prot[g]['home'][:20], rna[g]['home'][:20]))

    order = sorted(prot.items(), key=lambda kv: -kv[1]['top_share'])
    print('\n  most lineage-confined proteins in the antibody data:')
    for g, v in order[:10]:
        print('    %-10s %.3f  %2d types  home %s' % (g, v['top_share'], v['n_types'], v['home']))
    for i, (g, _) in enumerate(order, 1):
        if g == 'CD19':
            print('    CD19 ranks %d of %d mapped proteins.' % (i, len(order)))

    json.dump({'n_shared': len(shared), 'rho_top_share': float(rho), 'p_top_share': float(p),
               'rho_n_types': float(rho_n), 'p_n_types': float(p_n),
               'n_blood_cell_types_protein': len(types), 'n_blood_lineages_rna': nblood,
               'protein': prot, 'rna_blood': rna},
              open(os.path.join(HERE, 'protein_celltype_citeseq.json'), 'w'), indent=1)
    print('\nwrote protein_celltype_citeseq.json')


if __name__ == '__main__':
    main()


LINEAGE = {'b': 'B', 'plasmablast': 'plasma', 'plasma': 'plasma', 'cd4': 'T', 'cd8': 'T',
           'treg': 'T', 'mait': 'T', 'gdt': 'T', 'dnt': 'T', 't-cells': 'T', 'thymocytes': 'T',
           'nk': 'NK', 'nk_cd56bright': 'NK', 'ilc': 'NK', 'nk-cells': 'NK',
           'cd14': 'myeloid', 'cd16': 'myeloid', 'monocytes': 'myeloid',
           'macrophages': 'myeloid', 'kupffer': 'myeloid', 'microglia': 'myeloid',
           'neutrophil': 'myeloid', 'neutrophils': 'myeloid', 'granulocytes': 'myeloid',
           'eosinophils': 'myeloid', 'basophils': 'myeloid', 'mast': 'myeloid',
           'cdc1': 'DC', 'cdc2': 'DC', 'pdc': 'DC', 'asdc': 'DC', 'cdc': 'DC', 'pdcs': 'DC',
           'langerhans': 'DC', 'eryth': 'erythroid', 'erythrocytes': 'erythroid',
           'erythroid': 'erythroid', 'platelet': 'megakaryocytic',
           'megakaryocytes': 'megakaryocytic', 'hspc': 'stem'}


def lineage_of(name):
    n = name.lower().replace('-', ' ').replace('_', ' ')
    for k, v in LINEAGE.items():
        if n.startswith(k) or k in n.split():
            return v
    return None


def home_agreement():
    """The share and the home are two different claims. The two-factor rule of section 23.2
    depends mostly on the HOME -- is the cell type carrying this antigen expendable? -- and only
    secondarily on the share. So the two are checked separately."""
    r = json.load(open(os.path.join(HERE, 'protein_celltype_citeseq.json')))
    prot, rna = r['protein'], r['rna_blood']
    shared = sorted(set(prot) & set(rna))
    agree = tot = 0
    disagreements = []
    for g in shared:
        a, b = lineage_of(prot[g]['home']), lineage_of(rna[g]['home'])
        if a is None or b is None:
            continue
        tot += 1
        if a == b:
            agree += 1
        else:
            disagreements.append((g, prot[g]['home'], rna[g]['home'], a, b))
    print('\n' + '=' * 86)
    print('THE HOME VERSUS THE SHARE: two different claims, checked separately')
    print('=' * 86)
    print('  the two references put the antigen in the SAME broad lineage for %d of %d genes '
          '(%.0f%%)' % (agree, tot, 100 * agree / max(tot, 1)))
    print('  compare: the SHARE agrees only at Spearman rho = %.3f' % r['rho_top_share'])
    LAB = {t[0]: t[3] for t in TL.TARGETS}
    tg = [g for g in shared if LAB.get(g) in ('safe', 'toxic')]
    ta = sum(1 for g in tg if lineage_of(prot[g]['home']) and
             lineage_of(prot[g]['home']) == lineage_of(rna[g]['home']))
    print('  among the %d therapy targets measured both ways, the home agrees for %d.' % (len(tg), ta))
    print('\n  the disagreements (first 12):')
    for g, ph, rh, a, b in disagreements[:12]:
        print('    %-10s protein %-22s -> %-14s   RNA %-18s -> %s' % (g, ph[:22], a, rh[:18], b))
    print("""
  READING, and an earlier draft of this paragraph OVERSTATED it. Home agreement is 37% overall,
  not "good". It is better among the therapy targets (9 of 13), and much of the overall mismatch
  is a ROSTER problem rather than a disagreement: the two references do not contain the same cell
  types. The PBMC panel has no tissue macrophages, no neutrophils and no Kupffer cells, so a
  neutrophil protein has nowhere correct to land and gets assigned to whatever blood cell type
  carries the most background. Several listed disagreements (CD177 -> platelet, ANPEP -> cDC1,
  CD207 -> erythrocyte) are exactly that.

  What survives as a fair statement: antibody and transcript agree only weakly about HOW
  CONCENTRATED an antigen is (rho = 0.167), and agree better about WHERE it is when both
  references actually contain the right cell type. Neither is a strong endorsement, and the
  honest conclusion is in results section 32.""")
    json.dump({'home_agree': agree, 'home_total': tot,
               'targets_home_agree': ta, 'targets_total': len(tg)},
              open(os.path.join(HERE, 'protein_celltype_home_agreement.json'), 'w'), indent=1)


if __name__ == '__main__' and os.environ.get('HOME_AGREE'):
    home_agreement()
