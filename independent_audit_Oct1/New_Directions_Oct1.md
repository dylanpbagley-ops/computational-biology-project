# Strategic assessment: possible next directions

Prepared October 1, 2026. This is advisory analysis, not student-authored text for a fair submission. No historical results, labels, or locked predictions were changed.

## Recommendation

There is a plausible route from the present work to a substantially stronger project. The available evidence does not establish that this route will work, or support numerical award probabilities. The most valuable change is to replace retrospective gene-level safety ranking with a narrow prediction that can be tested against independent functional measurements.

Investigate a construct-specific CAR-T performance project first. The proposed question is whether measured surface-antigen distributions, combined with measured killing and expansion responses, predict sustained tumor control better than mean antigen expression or a static kill threshold. Select one antigen and one experimental system. A model that predicts short-term killing must not automatically be treated as predicting repeated-challenge control, clinical relapse, or toxicity.

If a suitable lab already has a validated EV capture assay, experimentally tested distribution-aware EV capture may be the preferable route. The decisive resource is an established assay with reproducible biological measurements, not access to a building or a general-purpose mentor. If neither route has usable data or an executable assay, choose a different question in a mentor's established system early.

## Why the present headline should not govern the redesign

The present normal-tissue score lacks drug construct, dose, route, exposure, and tumor burden. Several existing outcome labels combine biologically different interventions. These are structural limits of the question, beyond the statistical uncertainty of the original 30-target panel. A larger retrospective panel might improve a carefully bounded association study; it cannot supply missing drug-specific information by itself.

The old work remains valuable as preliminary research, code, provenance, and evidence about failure modes. It should not constrain the new study to preserve the strongest-looking old number. The corrected mathematical frontier concerns a fixed measurement and specified response class; it is not a universal impossibility result for therapies that alter antigen expression, trafficking, or recognition.

## Candidate routes

| Route | Potential contribution | Main risk | Early decision |
|---|---|---|---|
| Construct-specific CAR-T sustained performance | Predict a measured functional outcome, identify when a static antigen threshold fails, and test whether distribution information changes a useful decision | Raw protein distributions and matched functional outcomes may be unavailable; growth and killing may not be separately identifiable | Inspect actual files and metadata; reproduce one published experiment; determine whether a credible independent test is possible |
| Distribution-aware EV capture with new experiments | Predict yield and enrichment in mixed populations; test controlled heterogeneity and validate a selected capture condition on new biological preparations | Assay reproducibility, source identification, nonspecific adsorption, and vesicle measurement may dominate the signal | Demonstrate reproducible source-specific capture and an independently measured background before optimization |
| Agent-specific clinical safety benchmark | Produce a carefully reviewed dataset and measure the limits of atlas-derived safety prediction within a coherent modality | Too few independent outcomes, inconsistent reporting, selection bias, and disease/lineage confounding | Count usable studies after blinded outcome review and exclusions; establish whether a genuinely independent evaluation exists |
| New question within an established mentor assay | Maximize opportunity for a clear intervention, causal test, and new result | May have little connection to existing code; novelty needs separate review | Find a concrete unanswered question with available materials, instruments, expertise, and a short pilot |

These are strategic comparisons, not calibrated probabilities of either scientific success or winning.

## Evidence behind the CAR-T route

- Majzner et al. measured antigen-density-dependent responses and demonstrated that CAR architecture changes the recognition threshold. This supports construct-specific modeling, not a universal response curve. [Primary study](https://pmc.ncbi.nlm.nih.gov/articles/PMC7939454/).
- Barisa et al. compared three anti-B7H3 constructs, measured antigen densities, and performed repeated tumor challenges. Immediate cytotoxicity did not fully capture sustained function. The article states that source data are provided and some raw data are available on request. The required raw single-cell files have not been verified here. [Primary study](https://www.nature.com/articles/s41467-025-61427-4), [accessible full text](https://d-nb.info/1378246683/34).
- An AML modeling study supplies an inspectable repository containing model code and a data archive. Its data are useful leads for functional-model benchmarking, but do not by themselves establish matched antigen-density distributions. The repository uses Julia and includes a 2.4 GB data archive. [Repository](https://zenodo.org/records/12823159).
- Existing work already models CAR-T design and kinetics. A new differential-equation model alone is insufficient novelty. [CARCADE](https://www.frontiersin.org/journals/molecular-biosciences/articles/10.3389/fmolb.2022.849363/full).
- Single-cell antigen selection and Boolean logic-gate searches are established. A new atlas ranking alone would be an incremental direction unless it tests a distinct limitation and provides meaningful validation. [Single-cell CAR switches](https://www.nature.com/articles/s41587-023-01686-y), [SCAN-ACT](https://link.springer.com/article/10.1186/s13073-025-01514-9).

## What would make the preferred route strong

1. Establish a distinct contribution relative to the closest studies. Distribution sensitivity and antigen-low escape are already known; the project must quantify an unresolved decision or failure condition.
2. Obtain measurements in comparable units. RNA, antibody-barcode counts, CyTOF ion counts, and molecules per cell cannot be interchanged without justified calibration.
3. Compare a small sequence of models: mean-only, static density-response, distribution-aware static, and a minimal dynamic extension only where the measurements support it.
4. Evaluate measured killing, expansion, or surviving antigen-low populations. A ranking change computed entirely inside the new model is not a functional validation.
5. Separate development from evaluation at the biological level: donor, experimental run, construct, or study as appropriate. Different cells or wells from the same donor are not automatically independent biological validation.
6. Measure prediction errors, uncertainty, and decision consequences. Choose a primary metric and a scientifically meaningful improvement before evaluation. Use pilot data to determine sample size; do not choose an arbitrary number of cells as a substitute for biological replication.
7. Test a new experimental condition if a suitable lab can execute it. Published results already read during problem selection are discovery evidence, not blinded prospective validation.

The strongest possible finding would show both where the new method improves decisions and where its assumptions fail. A negative result is scientifically useful but does not automatically have the same competitive impact as a useful, independently validated prediction.

## EV experimental alternative

An executable experiment could compare populations with similar mean accessible-marker density but different heterogeneity, then predict capture across a small range of ligand densities. A controlled bead or liposome system could isolate the physical mechanism; actual EV mixtures would be needed for an EV application claim. Whether comparable size, concentration, accessible-marker density, and source-specific readouts can be achieved must be established experimentally.

For mixed-source EVs, measure recovery of target-source particles and capture of background-source particles separately. Total fluorescence or total captured particles alone cannot demonstrate selective tumor-EV enrichment. Source labels require controls for free dye, dye transfer, and label-dependent capture. New biological preparations are needed to evaluate reproducibility.

Bayesian optimization is optional. It is valuable only if it demonstrably improves a real measured objective relative to a fair baseline under the same experiment budget. Optimization algorithms cannot validate an unmeasured objective. Previous EV work has already established capture bias and affinity-capture methods: [capture bias study](https://link.springer.com/article/10.1186/s12951-021-00987-1), [EV-Elute](https://pmc.ncbi.nlm.nih.gov/articles/PMC11653091/).

## Proposed decision schedule

- **First two weeks:** confirm actual lab resources and mentor expertise; inspect the CAR-T candidate files; reproduce one relevant result; make a table of the closest competing methods and the proposed difference. Decide whether either main route has both inputs and an independent test. Do not spend these weeks constructing a large model first.
- **By early November:** complete one end-to-end pilot against an appropriate baseline. Decide whether there is a distinct, testable contribution worth pursuing. If measurements cannot identify the mechanism or support the claim, narrow or change the question.
- **November–December:** run the locked evaluation or experimental confirmation, assess uncertainty and failure cases, and reproduce the result on a separate biological preparation or dataset when available.
- **January:** prioritize the completed study, Dylan's own writing, notebook, figures, and interview understanding. Maintain the existing project as a fallback rather than needing to rescue an unfinished expansion.

The school deadline is recorded as late January in the project files. AzSEF currently lists senior judging on March 26, 2027. Improvements after the school fair should follow the applicable fair rules and approvals. [Official timeline](https://www.azscience.org/visit/events/arizona-science-engineering-fair-azsef/educators/).

## Competition judgment

A focused project with a useful original result, adequate controls, independent validation, and demonstrable student understanding can plausibly be competitive. Nothing in the current evidence establishes an inherent ceiling that makes an ISEF award impossible. Conversely, adding complexity, a prestigious lab, or a cancer application cannot establish high odds. Qualification and an award after qualification are separate outcomes.

The official judging criteria reward a focused question, appropriate methodology, sufficient data, reproducibility, creativity, and the student's understanding and independence. Presentation matters substantially, but good presentation cannot repair inadequate evidence. [Official criteria](https://www.societyforscience.org/isef/grand-award/criteria/).

Current unknowns: available lab and mentor capabilities, weekly time, budget, access to complete candidate data, and the final distinct contribution. These determine the recommended branch. The earlier ASU/mentorship discussion is historical context and has not been verified as a current arrangement.
