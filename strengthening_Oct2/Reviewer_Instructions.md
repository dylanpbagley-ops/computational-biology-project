# Independent endpoint review instructions

Advisory working aid, October 2, 2026. The accompanying `reviewer_form_score_free.tsv` contains 30 examples in a randomized order, candidate primary-source links, and blank judgments. It omits expression values, model scores, original outcome labels, and expected injury sites. Agent identities may still suggest outcomes to an expert; the form is score-free, not a guarantee of complete blinding. Source links were assembled from earlier AI extraction and are leads to verify, not endorsements of every source statement.

Do not send `review_key_do_not_send.tsv`, the diagnostic results, or the historical label-review sheet to an independent reviewer before their initial judgments. The reviewing person should have relevant clinical/pharmacology expertise, or identify where specialist input is needed. No form has been sent by the assistant, and no reviewer agreement is claimed.

## First resolve the endpoint before deciding labels

The historical binary definitions mix severe events, dose-limiting events, approval, and managed depletion. These do not form a consistently exclusive partition. A trial can be tolerated at one regimen and have dose-limiting toxicity at another. An approved product can have serious on-target adverse effects. The same antigen can support different outcomes with different agents. A grade 3 event can be excluded from a trial's protocol-defined DLT category.

The recommended repair is one record per **agent–study/regimen–event**, with separate fields for human evidence, event severity, dose-limiting status, affected normal tissue, author-attributed mechanism, and attribution confidence. Record absence of a reported event separately from affirmative evidence of tolerability. Unknown attribution stays unknown. Intentional therapeutic depletion and its clinical management must be recorded explicitly, not inferred from a blood-lineage label.

Choose a precisely defined primary endpoint with the reviewer before running an amended outcome analysis. Do not change that definition after seeing scores. Two possible endpoints answer different questions: documented dose-limiting on-target normal-tissue injury in a specified setting; versus observed serious on-target injury regardless of approval/management. Do not merge them into whichever rule preserves the best AUC. If binary classification is not defensible, retain the separate evidence fields and report that limitation.

## Row review

For each row, identify the actual agent and study; split a mixed-agent entry when needed. Check human versus preclinical evidence, event, severity, protocol definition, attribution, competing explanations, and adequacy of follow-up/exposure. Record the relevant source location and uncertainty. A boxed warning or a target's expression in an organ is not by itself proof that a particular injury was on-target. An abstract's silence about mechanism is not proof of off-target injury.

Record identity corrections separately from endpoint decisions. CEACAM5's Parkhurst study is TCR-T, not CAR-T; FAP's cited marrow/cachexia evidence is preclinical; the FOLR1 corneal attribution is contradicted by the cited mechanistic observations. Those are source corrections, but neither of the latter two automatically creates an adequate safe clinical record.

The initial score-free review should be returned before reconciliation with historical labels. Ideally obtain two independent reviews; record agreements, disagreements, and adjudication without claiming that two reviewers alone establish truth. One expert review is still useful if a second is unavailable. The assistant's previous knowledge of scores means its triage cannot count as blinded independent review.

## Expansion and evaluation

The existing source pool contains 166 rows across 79 target strings, including two dual-target strings; these are not 79 eligible independent clinical targets. New-source eligibility should use the same written endpoint and a defined sampling frame, with a cutoff date and explicit exclusions. Keep a record of all screened candidates and all rejected rows. Preserve whether targets appeared in development, earlier validation, or the locked predictions. Do not promote reused labels to independent validation by giving them a new file name.

When multiple studies share a target, keep them together in any training/test partition and account for clustering. Do not pretend repeated agents with the same expression score are independent new genes. If one target has inconsistent outcomes across settings, this is evidence about the limits of a gene-only predictor, not a reason to retain just one convenient study.

A wet lab is not required for this review. Access to a suitable reviewer and primary clinical literature is the immediate dependency.
