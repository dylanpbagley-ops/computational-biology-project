# The Pivot, Explained Simply

## What we're changing and why

Your project has spent months building a real, working piece of math: a model that predicts whether something gets "captured" or "killed" based on how many target molecules (receptors) sit on its surface — and the key discovery is that this relationship is **not a straight line**. Doubling the receptor count doesn't just double your chance of capturing it; it can multiply that chance much more steeply, because of how multiple weak attachments combine. That means the *average* receptor count across a population is a misleading number — two populations with the same average can behave very differently depending on how spread out (heterogeneous) their receptor counts are.

You built and tested this on extracellular vesicles (EVs) — tiny particles cells shed into blood — as a cancer early-detection idea: pick which surface marker (protein) to use to fish tumor EVs out of blood. The problem your parents correctly flagged: even if this works perfectly, early cancer detection from blood is bottlenecked by biology, not math. Very few tumor EVs actually make it into the bloodstream at an early stage. So a smarter marker-picking method can only move the needle a little, because the real limitation is upstream of anything your model touches.

## The pivot: same math, different target, much bigger lever

We already did deep research (that's the PIVOT-FOR-NOVELTY document) into where else this exact kind of math — "the shape of a population's receptor distribution, run through a nonlinear response curve" — would matter more. The clear winner: **CAR-T cell therapy target selection.**

Quick background: CAR-T is a real cancer treatment where a patient's own immune cells are re-engineered to hunt down and kill cancer cells that display a specific surface protein (an "antigen") — CD19 is the most common one, used against leukemia and lymphoma. It works great at first, but a big fraction of patients relapse because some cancer cells had *low* amounts of that antigen — too low for the immune cells to reliably detect and kill. This is called "antigen escape."

Here's the key fact we found: killing efficiency in CAR-T is also **nonlinear/thresholded** in antigen density, almost exactly like your EV capture model — below roughly 2,000 CD19 molecules per cell, killing basically stops working; above it, it works well. But when researchers currently decide *which antigen to target* for a new CAR-T design, they use crude measures: "what fraction of cells are positive at all" or "what's the average antigen level" — not "given the full spread of antigen levels across this tumor's cells, what fraction will survive because they're below the kill threshold?" Nobody has built a tool that does that second thing for target *selection*. We checked directly (not just assumed, and re-checked again this session with a fresh search) and reconfirmed this gap is real and still open.

This matters because: unlike early cancer detection (bottlenecked by biology no matter how good your math is), CAR-T target selection is a decision that gets made once, up front, when a therapy is designed — and picking a better target could mean **fewer patients relapse from escape**, which is a real, meaningful clinical lever, not just a marginally better blood test.

## What actually changes technically (and what doesn't)

The physics you already built is a near-exact match once you swap the words:

| Your EV model | Becomes |
|---|---|
| Vesicle | Tumor cell |
| Receptor count on a vesicle | Antigen count on a tumor cell |
| Capture antibody surface | CAR T-cell |
| "Probability of being captured" | "Probability of being killed" |
| Tumor vs. normal EV populations | Tumor cells above vs. below the kill threshold |
| Picking the best EV marker | Picking the best CAR-T target antigen |

Your core equation, your "superselectivity" (nonlinear sensitivity) math, and your dual-marker AND-gate logic (which already models "capture only if BOTH markers are present") all transfer directly — the AND-gate piece maps onto "logic-gated" CAR-T designs, a real hot topic where a T-cell is engineered to only attack if it sees two antigens at once (used to spare healthy cells that only have one antigen).

What's new: we swap in real published CAR-T data (measured antigen-density distributions per cell, and real dose-response killing curves from actual lab experiments) in place of the EV data, and rebuild the "kill probability" version of the model. We also validate it against a real public dataset (a genuine single-cell dataset from real CAR-T patients, confirmed publicly downloadable) to show the model reorders which target looks best compared to the standard "just look at the average" approach — the same kind of honest, checkable claim your project already made for EVs (Aim 9), but now for a decision that changes what an actual therapy gets built around.

## What happens to all your existing work

Nothing gets thrown away. The EV work stays in the project as the original proof that "distribution shape beats the average" is a real, validated phenomenon — it's the foundation the CAR-T work builds on and cites. The write-up will present it as: "we proved this effect is real using EV data, then applied it to a place where it matters more."

## The honest limits (so you go in with eyes open)

This will NOT cure CAR-T's other real problems — T-cells still get exhausted, still struggle to reach solid tumors, and antigen density isn't the only reason CAR-T fails. The honest claim is "a better way to prioritize which target to build a therapy around," not "a cure for antigen escape." We will state that plainly in the write-up, the same way you've insisted on stating every EV limitation plainly.

---

**Next response, once you approve:** I'll do the actual work — pull the real datasets, adapt the code, build and test the new "kill probability" model, run the validation against real data, and produce the updated write-up materials — following the full implementation plan I'm handing you alongside this document.
