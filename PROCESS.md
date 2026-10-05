# How this project's handoff process works (read once, reference as needed)

## The problem this solves

Every session ended by writing a new, complete handoff file (`EV-Capture-Project-Handoff-vN.md`) that included the FULL text of every prior version's "what changed" narrative. That was good practice for completeness — nothing was ever lost — but it meant the file grew roughly 15-20% bigger almost every session: 49KB at v3, 224KB by v15. Every new Cowork session had to read the whole thing before doing any work, even though most of it (the "What NOT to revisit" list, old bug-catch stories, closed decisions) would never be acted on again.

Concretely, this was costing real context budget on every single session, for content that was 95% settled history rather than 5% "here's what to do next."

## What changed (decided by the user, August 26, 2026)

Three options were presented: (1) split into a small living document plus an archive, (2) trim the single file in place by deleting/compressing old sections, (3) keep the full single file exactly as-is but tell new sessions not to read it by default. **The user chose option 1.**

## The new structure

- **`EV-Capture-CURRENT.md`** — the file you attach to a NEW Cowork session going forward, instead of the old `EV-Capture-Project-Handoff-v15.md`. Small (~15-20KB), contains only: who you are and how to work with the assistant, the project in one sentence, what changed in the MOST RECENT session, current state (done/not done), immediate next steps, known risks, what not to revisit (as a short list with pointers, not full explanations), a compact literature table, and a file index. No session-by-session narrative history.

- **`EV-Capture-Project-Handoff-v3.md` through `v15.md`** — kept exactly as they are, untouched, on disk. This is the archive. Nothing was deleted. If `CURRENT.md` says "see `v9.md` for the full sigma-bug story," that means: this fact is settled and doesn't need re-verifying, but if you ever want the blow-by-blow of how it was found and fixed, it's in that file.

- **Going forward, there is no `v16.md`, `v17.md`, etc. as full re-writes.** Instead: `CURRENT.md` gets updated in place each session (its "What changed this session" section gets REPLACED, not appended to — see below), and if a session's work is substantial enough to be worth preserving in full detail (a real bug catch, a real derivation, a framing decision with real back-and-forth), that detail goes into a NEW small archive file named for what it's about (e.g. `Gillespie-Cooperative-Model-Buildlog.md`), not appended to a single ever-growing handoff.

## How to maintain this going forward (for the next session that touches this file)

1. **At the start of a session**, only `EV-Capture-CURRENT.md` needs to be attached. If the task requires knowing the full history of a specific settled decision (rare — usually only if the user asks "why did we decide X" or a review/critique questions something), open the specific archived file `CURRENT.md` points to.

2. **At the end of a session**, update `CURRENT.md` in place:
   - Replace the "What changed this session" section with a new one describing THIS session's work (the old one's key facts should already have been folded into "Current state," "Known risks," etc. — don't just accumulate more "what changed" sections back-to-back, or this file will slowly turn back into the old growing document).
   - Update "Current state: done vs. not done" to reflect what's now true.
   - Update "Immediate next steps," renumbering as needed.
   - Update "Known risks" and "What NOT to revisit" only if something genuinely new needs adding, or something old needs marking resolved/removed.
   - If a real, detailed story is worth preserving (a bug and its fix, a hard framing decision, an external critique's diagnosis), write it as its own short new archive file, and add ONE line to `CURRENT.md`'s relevant section pointing to it — don't inline the full story into `CURRENT.md` itself.

3. **Never delete the old `v3.md`–`v15.md` archive files** — they're the historical record and cost nothing sitting on disk unread. The savings come entirely from NOT reading them into context by default, not from deleting them.

4. **If `CURRENT.md` itself starts creeping back up in size** (say, past ~30KB), that's the signal to do a "consolidation pass": move anything that's purely historical narrative (not currently-actionable) out to a new short archive file, and leave only a pointer behind — the same discipline as step 2, applied retroactively.

## Why this is expected to actually save context, not just look tidier

A new session that attaches `CURRENT.md` (about 15-20KB) instead of `v15.md` (224KB) starts with roughly 90% less handoff-document text to process before it can begin working — and what's left is entirely made of information that's actually relevant to deciding what to do next, not settled history it would have to read past. The archive is still there in full if a specific old story is ever needed, but it's opened on demand (via `Read`) rather than loaded into every session automatically.
