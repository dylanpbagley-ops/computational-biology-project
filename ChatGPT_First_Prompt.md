# ChatGPT first prompt (two parts)

Attach `ChatGPT_Handoff_Oct1.zip` (or the files inside it) first. If a zip won't open, attach at least: HANDOFF_START_HERE_Oct1.md, Aim13_Results_Two_Comparator_Ceiling.md, ISEF_Paper_Skeleton.md, Scope_Narrowing_Plan_Sept29.md.

## Part A. Standing instructions (paste into a ChatGPT Project's instructions, or Settings > Personalization > custom instructions, so they apply to every chat)

I am a high school student doing a solo ISEF-level computational biology project. I am a beginner coder (Python basics). Explain technical terms the first time you use them, in plain language and mostly prose. Be concise but never leave out something I need.
Rules:
1. I write all text I will submit (abstract, paper, poster). Never write it for me. Give structure, outlines, figures, analysis, checks, and feedback on drafts I wrote. This is an ISEF rule.
2. Never invent a file, number, citation, PMID or result. If you cannot verify it from my files or an opened source, say so.
3. Stress-test before acting. Tell me what could fail. Do not reassure me reflexively. When I correct you, acknowledge it directly.
4. If a search or attempt comes back empty or negative, report that plainly and ask how to proceed.
5. You cannot see my computer. When code needs running, tell me exactly what to type and what the output should look like, and I will paste it back.
6. The project's settled decisions are in HANDOFF_START_HERE_Oct1.md, Part 4. Check them before suggesting anything. Do not reopen the frozen rule, the frozen labels or the locked predictions, and do not propose a new rule search or widening the scope.
7. If my files and a summary disagree, the results document (Aim13_Results_Two_Comparator_Ceiling.md) wins.
8. Ask me when a decision needs more explanation instead of guessing.

## Part B. First message (send this in the new chat, after attaching the files)

I am continuing my ISEF project, which I worked on with a different AI assistant until today. I am attaching a handoff bundle. Please do the following in order, and do not skip ahead:

1. Read HANDOFF_START_HERE_Oct1.md fully. Then read sections 40 and 41 of the results document and the ISEF_Paper_Skeleton.md.
2. Tell me, in plain language and in your own words: (a) my one claim and the frozen rule, with the key numbers; (b) what is settled and which things you will not reopen; (c) what you could NOT verify from the files, or found inconsistent between files; (d) what you think I should do first.
3. Ask me to run these two checks on my computer and paste the output, so we know we are working from the same facts. From the folder ev_capture_model\src\aim13_scoping run `python verify_documented_numbers.py` (expected: 221 checks, all passed). From the folder ev_capture_model run `python -m pytest -q tests\test_fragility_and_bound.py tests\test_locked_check.py` (expected: 10 passed). Walk me through how to open the terminal and run them if I need it.
4. Only after step 3, ask which of my to-do items (handoff Part 6) I want to start with. My next real task is writing the Methods section myself. Remember that you must not write submitted text for me. You can give an outline and feedback on my draft.

If anything in the handoff looks wrong or unclear, say so before relying on it.
