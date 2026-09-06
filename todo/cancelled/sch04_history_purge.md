# Purge the Fujikura `sch04` Tables from Git History

**Date:** 2026-08-31
**Purpose:** The two brochure-derived critical-current tables `share/material/fesc-sch04.hdf5` and
`share/material/fysc-sch04.hdf5` were removed from the working tree on 2026-08-31 because they are
digitisations of copyrighted Fujikura plots and must not ship with the open-source release. A `git rm`
leaves both blobs reachable in history, so anyone who clones the public repository can still check
them out. This plan rewrites the three affected branches to remove the blobs, and lists the
coordination that rewrite requires.
**Module:** meta (repository hygiene / release)
**AIs involved:** Claude (scoping + plan)
**Status:** OPEN — awaiting Christian. The working-tree removal has landed; **nothing in history has
been rewritten and no command below has been run.** Decided 2026-08-31, Christian: delete outright
and purge history.

> **Scope guards:**
> - IN scope: the two blobs, on `main`, `devel` and `claude`, on `origin` only.
> - OUT of scope: the five figshare-derived REBCO tables and `bscco-2223` — those are CC BY 4.0 or
>   figure digitisations with recorded sources and stay exactly as they are.
> - OUT of scope: `tmp/bscco/**`, which holds the generator copies. `tmp/` is `.gitignore`d
>   (`.gitignore:17`), was never tracked, and therefore needs no history work — but it does hold
>   readable copies on this machine, so decide separately whether to keep them.
> - The rewrite changes every commit SHA from `568a0fff` forward. That is the cost, and it is why
>   this is a plan rather than something done in passing.

---

## 1. Current State

| Fact | Value | How checked |
|---|---|---|
| Blobs entered history | `568a0fff` ("backup"), `b99ba890` ("more bugfixes"), both 2026-08-29 | `git log --follow` |
| Removed from working tree | 2026-08-31, this session (staged `git rm`, not yet committed) | `git status` |
| Branches carrying them | `main`, `devel`, `claude` — all three also on `origin` | `git branch -a --contains b99ba890` |
| `origin` | `https://belfem.lbl.gov/gitlab/codes/belfem.git` — **carries the blobs on all three branches** | `git cat-file -e origin/<br>:share/material/fesc-sch04.hdf5` |
| `backup` (GitHub, `cmesse/belfem_backup`) | **does not** carry them on `main`/`devel`/`claude` as of the last fetch | same test against `backup/*` |
| Size | 9 228 172 B + 9 220 125 B ≈ 18 MB | `ls -la` |
| Referenced by code, decks, tests, schema | **No.** Only `share/material/README.md` named them | tree-wide `grep` |

**Bottom line:** two blobs, two commits, three local branches, one remote. The GitHub backup mirror
appears clean, which means the exposure is confined to the LBL GitLab — re-verify with a fresh
`git fetch backup` before relying on that.

## 2. Why `git rm` Alone Is Not Enough

`git rm` removes the file from the tip tree only. Both blobs stay reachable from `568a0fff` and
`b99ba890`, so `git checkout 568a0fff -- share/material/fesc-sch04.hdf5` still produces the file in
any clone, and any release cut as a git clone or a GitHub mirror carries them. Only a history
rewrite removes the object from every reachable tree.

If the public release is instead cut as a **tarball from a clean tip**, or as a fresh repository with
squashed history, the rewrite is unnecessary — see O1.

## 3. Ordered Steps

- [ ] **R1 — Commit the working-tree removal.** The `git rm` and the `share/material/README.md`
      edits are staged but uncommitted. Land them first so the rewrite has a clean tip to work from.
- [ ] **R2 — Back up the repository before rewriting** (after: R1). A full copy, not a clone —
      `cp -a` the whole checkout including `.git` to a location off this tree. A rewrite is not
      undoable once the reflog expires, and this checkout is shared.
- [ ] **R3 — Confirm no collaborator has unpushed work on `main`, `devel` or `claude`** (after: R2).
      The rewrite invalidates every SHA from `568a0fff` forward; anyone holding local commits on top
      of the old history has to rebase by hand. Gregory (`src/homology`) and Sirous both work in this
      tree — ask before rewriting, not after.
- [ ] **R4 — Verify the GitHub backup mirror is genuinely clean** (after: R2):
      `git fetch backup && git cat-file -e backup/main:share/material/fesc-sch04.hdf5` — a non-zero
      exit on all three branches means no GitHub exposure and no second rewrite. If it *is* present
      there, R6 has to be repeated against `backup`.
- [ ] **R5 — Rewrite with `git-filter-repo`** (after: R3, R4). Not `filter-branch` (deprecated,
      slow, and its `--index-filter` leaves replace-refs behind). `git-filter-repo` refuses to run
      on a repo with remotes configured unless `--force`, and it removes `origin` afterwards by
      design — re-add it in R6.
      ```bash
      # from a FRESH clone, not this shared checkout:
      git clone --no-local /home/christian/codes/belfem /tmp/belfem-purge
      cd /tmp/belfem-purge
      git filter-repo \
          --path share/material/fesc-sch04.hdf5 \
          --path share/material/fysc-sch04.hdf5 \
          --invert-paths
      ```
      A fresh clone is used deliberately so a failed rewrite costs nothing and the shared checkout is
      never in a half-rewritten state.
- [ ] **R6 — Gate the rewrite before pushing** (after: R5). All four must pass in `/tmp/belfem-purge`:
      - [ ] `git log --all --oneline -- share/material/fesc-sch04.hdf5` prints nothing (same for `fysc`).
      - [ ] `git rev-list --objects --all | grep sch04` prints nothing.
      - [ ] `git diff <old-tip> <new-tip>` over the three branch tips is **empty** — the rewrite must
            change history, not content.
      - [ ] The tip tree still holds the six remaining tables and `share/material/README.md` reads
            correctly.
- [ ] **R7 — Force-push the three branches to `origin`** (after: R6). `--force-with-lease` per
      branch. GitLab may protect `main`; unprotect, push, re-protect. Announce before and after —
      every collaborator must then re-clone or hard-reset, since a plain `git pull` will merge the
      old history straight back in and restore the blobs.
- [ ] **R8 — Expire the objects on the GitLab side** (after: R7). A force-push leaves the old
      commits reachable from GitLab's internal refs (and from any open merge request) until
      housekeeping runs. Trigger repository housekeeping in the project settings, and close or
      rebase any MR that still points at the old history. Until this is done the blobs remain
      downloadable through the GitLab web UI by SHA.
- [ ] **R9 — Re-point this shared checkout at the rewritten history** (after: R7). Simplest safe
      path is a fresh clone into a new directory and moving the old one aside once every worktree is
      accounted for. Note the build trees (`build/`, `cmake-build-debug/`) are shared across
      worktrees and will need reconfiguring if paths move.
- [ ] **R10 — Decide the fate of `tmp/bscco/**`** (after: R1). Untracked and `.gitignore`d, so not a
      release risk, but it holds the digitised source points, the generator, and the brochure-derived
      intermediates on disk. Delete or archive off-tree.

## 4. Open Design Questions

- **O1 — Is the rewrite actually required, or does the release form make it moot?**
  If the public release is a tarball or a fresh squashed repository, the blobs never reach the
  public and R5–R9 can be dropped, leaving only R1 and R10. If the release is this GitLab repo made
  public, or a GitHub mirror of it, the rewrite is required. **This question decides whether the
  rest of the plan runs at all — answer it before R2.**
- **O2 — Should the two tables be preserved anywhere?** Decided 2026-08-31, Christian: delete
  outright, no copy kept. Recorded here so a future session does not "helpfully" restore them from
  history. The tables can be rebuilt from the brochure if they are ever needed under a licence that
  permits it.
- **O3 — Does anything else in history carry the same problem?** This plan covers only the two
  tables. A release-hygiene sweep for other copyright-encumbered blobs (digitised figures, vendor
  data, `literature/` material that may have been committed by accident at some point) is a separate
  task and is *not* covered here.

## 5. Definition of Done

- [ ] R1 committed; `share/material/` ships six tables and a README that matches.
- [ ] O1 answered.
- [ ] Either R5–R9 complete with R6's four gates green, or O1 answered "release form makes it moot"
      and that answer recorded here.
- [ ] R10 resolved.
- [ ] `git rev-list --objects --all | grep sch04` empty in a **fresh clone of `origin`**.
