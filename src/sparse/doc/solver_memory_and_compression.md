# Solver Memory and Compression {#sparse_solver_memory_and_compression}

**Date:** 2026-08-16
**Purpose:** Choose `library` and `compression scheme` for a magnetic or
thermal solver block, estimate whether the run fits the machine, and know
when block-low-rank compression is safe
**Module:** sparse

This page is for a deck author writing `solver { linear magnetic { … } }`
or `linear thermal { … }`. It describes the **unified compression
contract of the 2026-08 campaign**: both STRUMPACK and MUMPS are exact
direct solvers unless the deck says `compression scheme : blr`. The
wrapper change landed 2026-08-16; binaries built BEFORE it silently
compress on MUMPS — see §4 and the pitfalls table.

## Short rule

Start production h-φ runs with STRUMPACK in exact mode. Switch to MUMPS
when STRUMPACK fails or robustness matters more than wall-clock. Never
enable BLR compression to "make it fit" before reading §4 — it is an
accuracy-for-memory trade, not a free speedup.

Messe et al. 2023 report exactly this split for BELFEM h-φ systems:
STRUMPACK typically up to ~2× faster (their Table 2 measures 2.1–2.8× on
one workstation), MUMPS more robust on ill-conditioned systems. Mixed h-φ
systems are not a job for an iterative method on the magnetic block
(Messe et al. 2023, §2.7 — the direct-LU rationale); the choice is which
direct library, not whether to factor.

## 1. The choice: STRUMPACK or MUMPS

| Need | Use | Why |
|---|---|---|
| Wall-clock, machine fits the factor | `library : strumpack` | ~2× faster (Messe et al. 2023, Table 2) |
| The run must finish; pivot trouble; unexplained direct-solver failure | `library : mumps` | The robust fallback (Messe et al. 2023) |
| Tight nonlinear target (`tolerance : 1e-11`) | either, **exact** | Exact factor + outer GMRES at rel_tol 1e-10 reaches ~1e-15 (measured 2026-08-18; on an ill-conditioned matrix the raw factor alone gives ~4 digits); BLR may not |
| Magnetic and thermal direct on one node | don't | Measured swap exhaustion — see §6 |

Preferred first pass:

```text
solver
{
    linear magnetic
    {
        library : strumpack ;
        // no compression line: exact factorization
    }
    nonlinear magnetic
    {
        tolerance : 1e-11 ;   // Messe et al. 2023 §4 — checkerboarding guard
    }
}
```

Robust fallback: change the one word `library : mumps ;`. Do that when a
STRUMPACK run died with `ZERO_PIVOT`/`NO_CONVERGENCE`, or when the
residual floors without an obvious physical reason.

## 2. Memory model

Peak RAM is **per rank**, not the node average. The kernel OOM-kills the
job from the fattest rank; averages are not the risk number.

**What sits on rank 0** (the expensive rank before any imbalance):

- the **full mesh** (workers do *not* carry a copy — they receive a
  partition plus aura, `cl_FEM_Kernel.cpp`, `partial_mesh()`);
- the **assembled global Jacobian** and a second global-sized system
  matrix (`cl_FEM_DofMgr_SolverData.cpp`, `allocate_matrices` /
  `collect_matrices`);
- its share of the factor — often the largest share;
- STRUMPACK only, per initialize: the MC64 matching **gathers the full
  matrix onto rank 0**. Leave `matching : on` (default) — the h-φ
  Jacobian hits exact zero pivots without it.

**What sits on every worker:** submesh (owned + aura), local assembly
structures, its factor share. MUMPS receives the matrix centralized on
the host (`ICNTL(18)=0`); STRUMPACK redistributes from rank 0.

**Fill imbalance is real** (Messe et al. 2023 §2.7 reported none
observed; production since has). Committed measurements:

| Observation | Measured |
|---|---|
| tapestack3d, ~832k dofs, 8 ranks, METIS_NodeND (pre-NodeNDP) | one rank ~23 GiB, siblings ~2.8 GiB |
| same case, `metis nodendp : true` (now the default) | rank 0 ~4.8 GiB, workers 1.8–2.1 GiB — factor *size* unchanged, *shape* flattened |
| remeshed tapestack3d ~2.5M dofs, 4 ranks, + thermal MUMPS | 17.8/13.4/10.5/9.8 GiB = 51.6 GiB RSS, swap exhausted |

`metis nodendp` changes the shape of the factor, not its size. Leave it
on. It does not replace a machine that is too small.

**OpenMP threads allocate workspace on the rank that is already
fattest.** Pin `OMP_NUM_THREADS` (the production script exports 2); a
hand-launched run that inherited more threads was part of an OOM kill
(campaign observation, 2026-08-16 — STRUMPACK traced `T=4`).
`hatch_turtle()` warns on rank×thread oversubscription; it does not
estimate the extra RAM.

## 3. Reading the solver's estimate before a long run

BELFEM prints no parsed memory line of its own; both libraries print on
their own streams when the log is loud enough. Raise verbosity:

```text
prterun -np 4 belfem -v 5      # STRUMPACK's own logger needs level 5
                               # MUMPS prints its analysis at level 4+
```

| Library | What you get | How to read it |
|---|---|---|
| STRUMPACK | `estimated memory usage (exact solver) = … MB`, printed after reordering, immediately before factorization | **A SUM over all ranks**, not per-rank. Divide by ranks only as a lower bound; the fattest rank carries more (see the imbalance table). A separate per-rank peak line exists upstream but is compiled out unless the library was built with flop counting. |
| MUMPS | its analysis table (needs `-v 4`+) | MUMPS estimates memory before allocating, but BELFEM calls `JOB=6` (analysis+factor+solve in one call; a `-9` workspace retry repeats factorization+solve as `JOB=5`), so allocation follows immediately. Watch the first solve and kill the run if the number is impossible. Since the 2026-08-30 retry work, `INFOG(1:80)` is copied back to C++ and drives the wrapper's error decode; the wider memory-report entries are still easiest to read from MUMPS's own stdout. |

A workable budget for the **fattest** rank:

```text
peak ≈ its share of the factor (estimate/ranks is the FLOOR, imbalance raises it)
     + rank-0 mesh + assembled global matrices   (workers: submesh only)
     + OpenMP workspace                          (scales with OMP_NUM_THREADS)
     + OS / desktop / other jobs
```

In practice the most reliable probe is empirical: launch, watch per-rank
RSS through the first factorization (`ps -o rss -C …`), and abort before
committing the night if the margin is gone.

## 4. Compression (BLR)

Block low-rank factorization stores off-diagonal blocks of the fronts as
low-rank products; the **drop tolerance is accuracy you are selling for
memory**. With BLR active, STRUMPACK's solve becomes GMRES over the lossy
factor — the delivered accuracy is the drop tolerance, not ~2e-16.

**Unified contract (this campaign):**

| Deck value | STRUMPACK | MUMPS |
|---|---|---|
| `automatic` (default) | off — exact | off — exact |
| `off` | off — exact | off — exact |
| `blr` | BLR, tolerance = `compression cutoff` (relative) | BLR, same cutoff into CNTL(7) (**absolute** — see below) |

> **Binaries built before 2026-08-16:** the old MUMPS wrapper
> mapped every value except `off` — including the default — to
> BLR-automatic with CNTL(7) = 1e-8. On such a binary, a MUMPS deck
> without `compression scheme : off` is silently lossy at 1e-8
> *absolute*. If an old MUMPS run's residual floors unexplained, suspect
> this first, and write `off` explicitly there.

**Why off by default:** on 2026-07-06 STRUMPACK's own auto-BLR (relative
1e-4) silently stalled Newton — the nonlinear residual sat on a
rank-dependent floor while every linear solve reported success. HTS decks
drive ε_n to 1e-11 to suppress checkerboarding (Messe et al. 2023 §4);
the exact-factor solve reaches that target through the outer GMRES
(~1e-15 at rel_tol 1e-10, measured 2026-08-18 — on an ill-conditioned
matrix the raw factor alone delivers only ~4 digits), and a
lossy factor eats the remaining headroom invisibly. The nonlinear loop cannot tell a physics
plateau from a compressed-solver floor.

**The two libraries do not read the tolerance the same way.** STRUMPACK's
BLR tolerance is *relative*; MUMPS's CNTL(7) is *absolute* ("not relative
to the input matrix … or the block norms" — MUMPS 5.7.3 §5.19; MUMPS's
default preprocessing scales the matrix, which softens but does not
remove the difference). One deck number therefore means two different
strictnesses. The number is a starting point; the exact-vs-BLR residual
comparison is the measurement.

**When BLR is sensible** — all of these, together:

1. the exact factor does not fit, *after* trying rank count, pinned
   threads, and `metis nodendp` (on by default);
2. the nonlinear `tolerance` is loose compared to 1e-11 — or you have
   knowingly given up on 1e-11;
3. the cutoff sits **2–3 decades below** the nonlinear tolerance
   (at-or-above the tolerance is a hard error; between the tolerance and
   two decades below it, a warning — exactly two decades passes
   silently):

   | Nonlinear `tolerance` | Cutoff to try | Verdict |
   |---|---|---|
   | 1e-6 | 1e-8 … 1e-9 | reasonable first BLR test |
   | 1e-8 | 1e-10 … 1e-11 | measure against exact mode |
   | 1e-11 | 1e-13 … 1e-14 | usually leave BLR off — roundoff is close |

4. a short `-v 5` probe confirmed the compressed estimate actually
   dropped enough to matter.

**The headroom check (this campaign, two tiers):** pairing `blr` with a
cutoff at or above the nonlinear tolerance — zero or negative headroom,
whether you stated the cutoff or inherited the 1e-8 default — is a
**hard error** at setup: the same impossibility the PETSc gate refuses
(`input_file_reference.md` §4.1). A cutoff below the tolerance but with
less than two decades of headroom prints a **warning** and proceeds —
that is the legitimate memory-bound trade, eyes open. The old silent
default against a 1e-11 target now errors: that pairing is the
2026-07-06 failure with the numbers shifted.

## 5. Leaf size

STRUMPACK's BLR leaf size — the dense block at the bottom of the
compression tree — is auto-scaled from the matrix dimension when `blr` is
requested (`strumpacktools.cpp`): 128 below 100k rows, 256 below 250k,
512 above. Larger leaves: fatter BLAS kernels, less compression; smaller
leaves: more compression opportunity, more overhead. The heuristic exists
so a 20k-dof block is not tiled like a 2M-dof Jacobian. There is no deck
key, deliberately.

**MUMPS has no leaf-size equivalent, anywhere.** Its BLR clustering is
internal, derived from K-way partitioning of the ordering graph (which is
why the MUMPS manual insists on METIS/Scotch orderings for BLR). The
remaining MUMPS BLR controls — variant choice ICNTL(36), contribution-
block compression ICNTL(37) — are not exposed by BELFEM and are not leaf
sizes.

## 6. Practical sizing workflow

1. **Write the deck exact** — `strumpack`, no compression line.
2. **Pin threads** — `export OMP_NUM_THREADS=2` (or launch via the
   production script that exports it). Never inherit from the shell.
3. **Probe at `-v 5`** — one factorization; read the estimate (§3),
   remember it is a sum; watch per-rank RSS.
4. **Size from the fattest rank** plus the rank-0 overlay plus desktop
   margin. If one rank is 2× its siblings, that rank is the budget.
5. **If it does not fit**, in order: free the desktop; re-check threads;
   change `np` **and re-probe** (fewer ranks = larger shares but a
   sometimes-kinder front — empirical, not arithmetic); try
   `library : mumps` and read *its* estimate (robust, not automatically
   smaller); only then `blr` with a §4-compliant cutoff.
6. **Never two direct factors on a tight node.** The measured case:
   magnetic STRUMPACK + thermal MUMPS on a ~62 GB-class workstation →
   51.6 GiB RSS, swap exhausted, and the thermal answers were wrong
   before the memory died. Keep thermal on PETSc unless a spare
   node is measured.

**Case study (tapestack3d production, 2026-08, campaign observations):**
~2.5M magnetic free dofs. At 4 ranks the fattest rank's factor share was
observed near 21 GiB with peers at 10–13 GiB; with the rank-0 overlay and
a live desktop the node had no margin, and the window produced repeated
OOM kills — among them one traced to inherited `T=4` OpenMP threads and
one to the dual-direct experiment above. The mitigation was `np=3` plus
pinned threads, chosen by re-probing, not by assuming fewer ranks are
smaller.

## 7. Memory controls: one deck key, the rest wrapper policy

| Feature | Status |
|---|---|
| MUMPS out-of-core (ICNTL(22)) | never set; in-core only (the error table knows `-90` because MUMPS can raise it) |
| MUMPS memory cap (ICNTL(23)) | **deck key `memory budget`** (MB per process, MUMPS only) — the one row in this table that *is* a deck key. Stated, it applies from the first factorization. Unstated, the wrapper measures the machine once at instance creation — before MUMPS has allocated anything, so the number is the total the instance may take, which is what `ICNTL(23)` bounds: available memory (`MemAvailable` capped by the tightest cgroup limit at the standard mount; free + inactive pages on Darwin) / ranks on the node / live MUMPS instances × 0.5, minimum over all ranks — and applies it only after a first `INFOG(1) = -9` / `-8`: at or above MUMPS's own estimate (`INFOG(16)`, `INFOG(36)` under BLR) it becomes `ICNTL(23)` and the factorization repeats; below it the step goes to the controller at once, no ladder. `-19` (the cap cannot be met) is the give-up code. Unmeasurable on any rank → the ladder below, as before |
| MUMPS workspace relaxation (ICNTL(14)) | starts at 30 % in the wrapper. On `INFOG(1) = -9` or `-8` *after* the cap above is in place — or instead of it, when the machine could not be measured — and on `-17` / `-20` (MPI send / reception buffer too small, sized from `ICNTL(14)` alone, so the ladder is their only remedy) the wrapper retries the factorization itself, doubling the relaxation up to a 480 % ceiling (30 → 60 → 120 → 240 → 480) before handing the failure to the timestep controller. MUMPS documents that `-9` can recur with `ICNTL(23)` set and still asks for a larger `ICNTL(14)` then. The escalated value and the measured cap persist for later solves on the same solver instance; `free()` restores 30 % and clears the measured cap. The relaxation is not a deck key |
| MUMPS BLR variant / CB compression (ICNTL(36)/(37)) | never written |
| MUMPS `INFOG` memory report as a C++ accessor | `INFOG(1:80)` is copied back (the full library width; the BLR estimate is `INFOG(36)`). The error code, the supplementary value and the memory estimate drive the retry; the soft-fail box prints estimate, cap and shortfall. No named accessor exists yet |
| STRUMPACK HSS/HODLR | no deck value selects them — but `--sp_compression` on the command line is applied **after** the deck and can re-enable them. Do not pass it on a Newton loop. |
| STRUMPACK BLR leaf size | auto-scaled (§5), no key |

`relative tolerance` in a linear section is **not** the compression
cutoff: for PETSc it is the Krylov stop, for STRUMPACK the exit
threshold of the outer GMRES over the exact factor (always applied
since 2026-08-18, class default 1e-10; before that stated-only, with
an unstated deck inheriting the library's 1e-6). Different knob,
different section of §4.1.

## 8. Pitfalls

| Pitfall | Symptom | Action |
|---|---|---|
| MUMPS deck without `compression scheme : off` on a pre-unification binary | silent 1e-8-lossy factors | write `off` explicitly |
| BLR under `tolerance : 1e-11` | Newton stalls/floors with "successful" linear solves | exact factors, or a cutoff ≤ 1e-13 (two-plus decades — passes the gate silently; 1e-12 warns, ≥ 1e-11 errors) |
| Sizing from estimate × anything | OOM on one rank while the average looks fine | fattest rank + overlays (§3) |
| Assuming mesh replication, so fewer ranks free RAM | `np=3` can be *larger* per rank than `np=4` | workers hold partitions; re-probe |
| Unpinned OpenMP threads | extra workspace on the fat rank | export the script's value |
| Second direct factor on the thermal block | swap exhaustion + wrong thermal answers | PETSc thermal |
| `--sp_compression` on the command line | silently overrides the deck | don't |
| `matching : off` to save the rank-0 gather | `ZERO_PIVOT` on the h-φ Jacobian | leave it on |

## References

- Messe et al. 2023, *Supercond. Sci. Technol.* **36** 114001 — solver
  comparison and Table 2; §2.7 direct-LU rationale; §4 checkerboarding
  and ε_n = 1e-11. (Full citation: `doc/literature_references.md`.)
- Amestoy, Duff & L'Excellent 2000 (MUMPS); Ghysels et al. 2016
  (STRUMPACK) — as cited by Messe et al. 2023.
- MUMPS 5.7.3 Users' Guide §5.19 (BLR: ICNTL(35), CNTL(7) semantics,
  internal clustering) — source of the manual facts above.
- Input contract: `doc/input_file_reference.md` §4.1 and
  `doc/input_schema.yaml`.
