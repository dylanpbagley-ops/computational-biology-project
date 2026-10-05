"""Documented DAMAGED normal tissue for each toxic target, mapped to HPA cell types (frozen before scoring; biology and the label evidence text only).
The damaged tissue comes from the 'evidence' column of target_labels.py / target_labels_validation.py, not from any expression value."""
DAMAGED = {
 'ERBB2':   ['alveolar cells type 1', 'alveolar cells type 2', 'transitional alveolar cells', 'cardiomyocytes'],   # fatal pulmonary event; cardiotoxicity
 'EGFR':    ['basal keratinocytes', 'suprabasal keratinocytes', 'alveolar cells type 1', 'alveolar cells type 2', 'transitional alveolar cells'],  # skin and lung
 'CEACAM5': ['colonocytes'],                                                                                       # severe colitis
 'CA9':     ['cholangiocytes'],                                                                                    # bile duct epithelium (liver enzyme toxicity)
 'CD47':    ['erythrocytes', 'erythrocyte progenitors'],                                                          # anaemia
 'CD44':    ['basal keratinocytes', 'suprabasal keratinocytes'],                                                 # fatal skin toxicity
 'NCAM1':   ['schwann cells', 'other brain neurons', 'brain excitatory neurons', 'brain inhibitory neurons'],    # peripheral neuropathy
 'EPCAM':   ['hepatocytes', 'cholangiocytes'],                                                                    # hepatotoxicity
 'FOLR1':   ['ocular epithelial cells', 'conjunctival goblet cells'],                                            # corneal keratopathy
 'FAP':     ['fibroblasts', 'fibro-adipogenic progenitors'],                                                     # marrow stroma
 'F3':      ['ocular epithelial cells', 'conjunctival goblet cells'],                                            # ocular toxicity (boxed warning)
 'GUCY2C':  ['enterocytes', 'colonocytes', 'goblet cells', 'paneth cells', 'enteric stem cells', 'enteric transient amplifying cells', 'tuft cells'],  # intestinal mucosa (validation set)
 'PSCA':    ['urothelial cells'],                                                                                 # cystitis (validation set)
}
