---
name: status_report
description: Give a structured end-of-task status report — tables, real timestamps, pass/fail status, and a clear bottom line. Use whenever Robert asks for a status report, a summary of a multi-part task, or "where are we" on RiOMar work.
---

# Status report

Produces the same style of report used throughout RiOMar sessions: dense,
scannable, structured with tables rather than prose, every timestamp and
count pulled live rather than guessed, and a one-line bottom-line verdict at
the end. Default output is plain markdown text in the chat reply — reach for
an Artifact only if Robert asks for something shareable/browsable.

## When to use this

Invoke whenever asked for:
- "give me a status report"
- "what's the status of X"
- a summary after a multi-part task (renaming sweep, pipeline sanity check,
  a batch of fixes) — anything with several independent sub-tasks that each
  have their own pass/fail/found-a-bug outcome
- "where are we" / "are we ready to run this"

Not a fit for a single, simple change — that just gets a one-sentence
confirmation, not a report.

## Before writing the report

1. **Get the real timestamp.** Run `date '+%Y-%m-%d %H:%M %Z'` (or the
   platform equivalent) and use that verbatim. Never guess or extrapolate a
   date from how much has happened in the conversation — long sessions can
   span far less (or more) wall-clock time than they feel like. If any
   earlier text in the session asserted a date, cross-check it against this
   command before reusing it.
2. **Re-verify, don't recall.** For every claim that will appear in the
   report (a file exists, a function is defined, a test passed, a count of
   something), re-run the actual check now — `grep`, `ls`, `Rscript -e
   'parse(...)'`, `python3 -c 'import ast; ast.parse(...)'`, whatever proves
   it — rather than trusting what you remember concluding five tool calls
   ago. State only what you just confirmed.
3. **Separate "fixed" from "found but not fixed."** If the task included an
   explicit, named thing to fix, and you found *other* unrelated bugs along
   the way, do not silently fix those too unless asked — report them
   clearly as a separate, flagged finding (what's broken, why, what it would
   take to fix), so nothing gets mistaken for already handled.

## Structure

Use whichever of these sections are relevant to the task — skip any that
don't apply, don't pad. Keep prose inside a section to 1-3 sentences;
prefer a table over a paragraph whenever there are 3+ comparable items.

```
## Status Report — <date from step 1> <time> <timezone>

### 1. <first sub-task area>
<one-line framing, then a table if there are 3+ comparable rows>

### 2. <second sub-task area>
...

### N. <anything found but not fixed, flagged as its own item>
**Found, not fixed** (<why not in scope>): <what's broken> — <concrete
consequence, e.g. "this will crash Stage 3 once it reaches step X">.

### Bottom line
<one or two sentences: is the thing actually ready/working, what's the
single next action, and any explicit recommendation for how to proceed>
```

## Table conventions

- Use a real markdown table (`| col | col |`) any time you're reporting
  status across a numbered/named set of things (scripts, files, functions,
  zones) — not prose like "Figure 1 is fine, Figure 2 is fine, but Figure 3
  has an issue."
- A status column should use a short, consistent vocabulary across the
  whole table: ✅ / ❌, or `OK` / `MISSING` / `STALE`, or `wired` / `gap` —
  pick one scheme per table and hold it, don't mix.
- For "old name → new name" or "before → after" tables, put the old value
  first, new value second, and a third column for what it actually produces
  or does — the reader should be able to tell *why* the change happened
  from the row alone, not just *that* it happened.

## What NOT to do

- Don't narrate the process ("first I checked X, then I looked at Y, then I
  realized Z") — a status report states current, verified facts, not a
  chronological account of how you got there. Save process narration for
  in-the-moment updates during the work itself, not the final report.
- Don't restate things the user already watched you do earlier in the same
  turn (e.g. don't re-explain a fix in full prose right after showing the
  diff) — the report is the compressed summary, not a second telling.
- Don't inflate a one-line confirmation into a full report structure. If
  the task was small, say so in a sentence or two and stop.
- Don't fabricate or round timestamps, counts, or file sizes "for the
  reader's convenience" — every number in the report must trace to a
  command you actually ran this turn.
