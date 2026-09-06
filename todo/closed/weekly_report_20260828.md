# Weekly Report, 22 to 28 August 2026

**Date:** 2026-08-28
**Purpose:** what the project did in the seven days before the release date, distilled from the
68 devlog entries of that window
**Scope:** open-source tree only; `nonfree/` work is recorded in `nonfree/devlog`

---

## 1. The week in one sentence

This was the week the release stopped being a target and became an object. BELFEM learned to
build as one shared library and install itself, gained a unified solver executable, retired 40
debt-register rows, and, in the same seven days, took delivery of a working ngspice parser, a
gauge-operator battery one element from complete, four new fluid property models, and the root
cause of three defects that had been silently corrupting results in production.

## 2. By the numbers

| | 21 Aug | 28 Aug |
|---|---|---|
| devlog entries in window | — | **68** |
| commits (all branches) | — | **45** |
| debt register: live rows | 40 | **41** |
| debt register: retired rows | 48 | **88** |
| highest debt ID | DR-91 | **DR-132** |

Forty-one new debt IDs were filed (DR-92 through DR-132) and forty rows retired, so the live set
grew by one. "Retired" counts struck rows wherever they sit: 84 in `debt_register_closed.md` plus
4 struck but not yet archived, because archiving lags striking by about a day. That is the honest
shape of the week. The register is holding level, not because the work is finished, but because
the work is finding things as fast as it closes them.

## 3. Release engineering

The build system was rebuilt around the release. `libbelfem` is now **one** library, assembled
from per-module OBJECT libraries. It is static by default and shared under
`USE_SHARED_LIBS=ON`; `make install` deploys the headers, `share/`, `examples/` and the named
executables. Both platforms were gated by execution, first macOS x86_64 (RPATH via
`@loader_path`, install-time re-sign, compiled-in data directory), then Linux
(`libbelfem.so.0.9.0`, clean `ldd`, an installed deck running with `BELFEM_DATA` unset, and a
`DESTDIR`-staged tree resolving through `$ORIGIN`). The `::` in the build-tree RUNPATH was
chased down and cleared as CMake's own toolchain fencing.

Two things fell out of that work. First, the user-material plugin templates could not compile
against an installed tree at all: their include path pointed at `include/`, while the headers
install under `include/belfem/<module>/`. Second, their documented workflow had **never** worked.
`BELFEM_DIR` was a `CACHE` variable without `FORCE`, so the placeholder written during the first
configure shadowed every later edit. DR-37's backend-free compile contract became a checked-in
build gate in the same arc, and its long-open R13 step was then closed by deletion after the
round proved that the dead block had been in the wrong file all along.

Separately, the two h-φ solver drivers began collapsing into one. The new `belfem` executable
picks magnetic-only or coupled h-φ/T **from the deck** (an unlabeled thermal solver section means
coupled), defaults a missing temperature to 77 K in both modes, and refuses thermal BCs without a
thermal solver. `hphirun` and `hphiTrun` are frozen through the release rather than deleted.
It landed after a full three-vendor plan and code round, which followed the gantry `np=8`
first-timestep SEGV. That crash was root-caused to a jump through an uninitialized `mFunMKF`,
and the "works in the debugger" theory collapsed into an executable mismatch.

## 4. The solver at the quench front

The week's longest thread was the quench front, and it ended with a sharper picture rather than a
solved problem. The Picard branch of `SolverData::solve` now solves in **increment** form, so an
iterative solver's relative tolerance measures the increment rather than the absolute 77 K field
offset. The A/B is decisive: 4.4e-2 K drift becomes 1.4e-9 K, and the thermal floor moves from
about −52 dB to −108.7 dB through an iterative solver. A reset seed hole was then found, fixed
and validated in production within the hour. Soft reset restored fields but never re-seeded the
dofs that the Picard residual actually reads, so the fix retired the whole restart ladder. The
Coulomb-gauge penalty `chi = 1e-4` was ruled default ON, with the input contract updated
in-session, and its Newton cross-term was planned, audited, implemented and audit-corrected in a
single evening.

**The result that matters most is a negative one.** With the gauge on and thermal moved to MUMPS,
the configuration ladder walled at 4250.9 ms with thermal converged and magnetic Picard
non-contractive at 2 µs: Picard is structurally insufficient at the quench front. A
BDF1-on-reject proposal was then jury-rejected on in-house counter-evidence, but the deck-level
BDF1 A/B **outran BDF5 plus Newton** through the front once the reset defect was gone, which
vindicates the scheme the paper started with. Nobody expected the low-order scheme to win.

The gauge campaign supporting all of this also finished its element work. **Nine of ten**
G-operators are complete (TET4, TRI3, PENTA6TS, QUAD4TS, HEX8, TRI6, TET10, HEX8TB, HEX8TS), each
through a four-stage blind-audited ladder with executed break gates. Only LINE3 remains, and it is
a contract decision rather than an implementation. The campaign paid for itself immediately:
`EF_HEX8::C()` had been returning the **negated** curl for years, invisible because HEX8 had zero
battery coverage and `K = CᵀC` is sign-invariant.

Still open: the `tapestack3d` j/jc noise at t = 2.125 s was root-caused as a ratcheting current
reorganization on a magnetic system that runs exactly singular, with roughly nine unpinned φ gauge
constants, advanced by always-accept Picard semantics that certifies the pre-update iterate and
time-steps the unchecked one. The mitigation remains open (DR-126/127/128), the formulation
decisions are Christian's, and the material table and piecewise resistivity law were both refuted
as the trigger, unanimously across the three voices.

## 5. Three silent corruptions

Three defects were found this week that had been producing wrong or missing results in production
without ever announcing themselves.

`face_key_3d` had been keying faces by the three lowest-indexed nodes of the full facet set. That
key is non-injective at order 2, which means **no order-2 gmsh mesh had ever loaded**.

`Basis::mNumberOfSources` was a `uint8_t`, so the gantry's 464 cuts needed 465 sources. The value
465 wrapped to 209, the **original node was dropped**, and the cut constraint degenerated to a
constant. That was the dead yoke arm Christian had read off a plot. The proof was a law rather
than an argument: 32 distinct N values, every group pinned to exactly `((N+1) mod 256)·I`. His own
question, "what if we go up to 16 bits, I only used `uint8_t` to save memory", got a measured
answer. The byte was never saved, because every one of those counters sits in alignment padding
ahead of a pointer.

Third, an orphaned `[DIAG]` dummy `Matrix` send in `Postprocessor::recover_fields()`, live since
July, had been poisoning the MPI comm fabric. A rank owning zero postprocessor nodes sent twice
what root collects once; the unmatched size header then truncated the *next* collect. It crashed
the shipped `examples/RLC_Circuit` under its own `./Allrun` at np=4 and the gantry at np=10, and
it is fully deterministic on partition layout, which is why np=2 and np=8 looked immune. Both
decks are red-to-green gated.

Around these: thin-shell facet nodes were kept out of the master re-derivation path, and the
parallel regression that fix caused was found and fixed the same evening. Pure half-cut periodic
ties landed the cap-corner fix. Inert coil blocks were kept out of the magnetic equation, which
let `costheta` run from its mesh for the first time. One mesh global per boundary condition
replaced one per tape, collapsing 464 identical entries to one while 2D_Undulator's genuinely
different ±0.0467 A pair survives. A double jury removed roughly 1200 lines of dead code from
`src/homology` (−1201/+135). Plain QUAD4 was unblocked in the Calculator, and Christian's
rectangular-only rule for QUAD and HEX edge elements was documented against Monk and Boffi. The
progress bar's one-burst behaviour was traced to *two* independent buffers, only one of which a
plausible story would have caught.

## 6. Materials, gas physics and the circuit path

The materials and gas side ran almost as a parallel project. The completed `EoS_Nitrogen` had six
derivative defects fixed, and its saturated liquid density problem was traced to a single wrong
digit in `mJ(19)`. Lemmon & Jacobsen transport for N2 and O2 arrived as one fluid-parameterized
class. H2 viscosity came from Muzny 2013 **with** its 2022 erratum, and H2 conductivity from
Assael 2011 with spin-isomer awareness. All reproduce the published verification tables to the
printed digits.

Elasticity data for pure metals arrived and was immediately jury-reviewed into a rebuild. The
Callaway optical channel was found dimensionally dead, and the non-metal thermal expansion branch
was refitted (α at 20 K fell by 300×, 7× and 15× for magnesia, Hastelloy and YBCO). Both gastables
and gasmodels went const-correct, with the compiler as oracle. All seven materials documentation
files were rewritten from a fact sheet, passed through blind audits, and received a prose sweep,
with every `file:line` citation now generated from anchors.

The circuit path went from question to shipped in a day. DR-39's parser question became
**ngspice parser v1**, comprising number parsing, a netlist parser, a factory, and a hybrid
`circuit{file:}` deck path. `examples/circuit` now runs from `tapestack.cir` in production, with
node indexing proven identical to the classic deck. Three P0s from the plan audit and one
pre-existing P0 in `compute_MNA_matrix` were caught on the way.

Two physics adjudications worth recording: the gantry's J/Jc of order 1000 was ruled a modeling
and deck problem rather than a postprocessor defect (the solver had already quenched a 115 A-class
tape driven at 340 A), and the deck was realigned to its design paper. And `nonlinear thermal
{ max iterations }` was found silently inert in fully-coupled mode, now enforced under a
conservative predicate.

## 7. What the method caught this week

Three findings are worth keeping as method rather than content.

- **A closure record is a claim about the tree, and only the tree settles it.** A landing-claim
  sweep over 38 register rows asked one mechanical question: does the tree contain the edit that
  the status cell says landed? It found 11 checkable claims and 1 false claim. The failure had a
  recognizable shape. It was the only landing split across two edits decided in the same breath,
  one of them a hazard guard on which the auditors had disagreed.

- **A gate re-run against a drifted deck is not the same gate.** Preparing DR-87's overnight run
  revealed that the deck's magnetic solver had been switched from STRUMPACK to MUMPS during an
  unrelated comparison, while the row's evidence rested on the STRUMPACK configuration. The run
  would have produced hours of numbers unable to support the strike for which they were run.

- **An unreachable gate can be discharged by substituting a deck that exercises the same code
  path.** DR-30's nominated corc gate was structurally unreachable. A bulk-superconductor deck put
  a number on the same claim (+10.8 % redundant aura reads, 0 stale elements).

The register itself was rebuilt to match. It was split into a live file and an archive, its
233-line pass-by-pass preamble was replaced by a rules document, and every open row received a
freeze tag: `[F]` before the design freeze, `[P]` after it, and `[W]` for latent hazards with no
trigger. Two rows carry `[F]` today. The `todo/` directory was swept from 47 active files to 35
on 26 August, and has grown again since as new plans landed. The two FVM plans moved to
`nonfree/todo/`, because the module they describe no longer exists in this tree.

## 8. What is owed

The release-blocking list is short: DR-89 (the magnetic solve exits at STRUMPACK's library-default
absolute tolerance) and DR-65 (a shipped example dead since April). Everything else is `[P]` or
`[W]`. The larger debt is not in the register but in the gate column: a substantial fraction of
this week's work is **reviewed, not verified**. Owed gates include the `make check-fast` runs
behind several landed fixes, the gantry `np=10` reruns, the coupled-mode reference match for the
new `belfem` executable, the one-key A/B ladder for the jjc noise, the macOS gate for the plugin
templates (DR-131), and the production-deck forced parity for increment-form Picard. None of these
is expected to fail. All of them are still owed.

---

*Sources: `devlog/dl20260822_*.md` through `devlog/dl20260828_*.md` (68 entries),
`todo/debt_register.md`, `todo/debt_register_closed.md`, git log. Codex ran the prose sweep;
Grok ran a fact and completeness audit, which refuted four claims in the first draft (the live-row
count, the DR-87 drift direction, the driver-collapse count, and the section 1 summary) and
identified the DR-100 omission. Every refutation was re-verified against the tree before it was
applied.*
