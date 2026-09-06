# tapestack3d j/jc Noise Onset at t = 2.125 s — Investigation Report

> **CLOSED 2026-08-30 — SOLVED, moved out of the active set by the register/todo currentness sweep.**
> The file's own Status has said CLOSED since 2026-08-29: the fixes it produced (DR-126 gauge autopins,
> DR-127 certified exit) were built and the deck rerun with the noise gone, run-verified. Its only residue,
> the low-field axis question, is owned by `todo/dr128_spap_low_field_axis.md`.

**Date:** 2026-08-28
**Purpose:** Root-cause investigation of the sudden j/jc noise in the fine-cadence
tapestack3d run (piecewise law), with jury round and recommended mitigation
**Module:** fem/kernel (controller, dof manager), fem/maxwell, physics/materials
**Method:** read-only diagnosis; blind three-AI jury (Codex + Grok) on a
pre-registered brief; every auditor citation re-verified against the tree.
Full record: `tmp/ai_exchange/review_tapestack3d_jjc_noise.md` (ephemeral).
**Evidence level:** numeric probes on the run's own output + static source
traces. Nothing in sections 1-8 is run-verified beyond those probes; the
executable gates are listed at the end.
**Status (2026-08-29): CLOSED — see Addendum 3.** The fixes this investigation
produced (DR-126 gauge autopins, DR-127 certified exit) were built and the deck
rerun: the noise is gone, run-verified.

---

## 1. What actually happened (measured, not inferred)

Run: `cmake-build-debug/tapestack3d`, `belfem` built Aug 27 23:07
(efe46051 + uncommitted), deck: h-φ/T coupled, 8-tape stack, sp-ap table,
`resistivity type : piecewise`, tolerance 1e-7, `algorithm : Newton`,
`anderson stabilization : on`, BDF1, save every 25 ms.

The "noise" is not cosmetic. Between keyframes 84 and 86 (t = 2.100 → 2.150 s):

| t [s] | max j/jc | max B [mT] | tape-1 ⟨\|Jx\|⟩ [A/m²] | I_1 [A] |
|---|---|---|---|---|
| 2.075 | 0.047 | 3.14 | 3.6e6 | −6.2e-5 |
| 2.100 | 0.049 | 3.28 | 3.8e6 | −6.6e-5 |
| 2.125 | 0.061 | 3.54 | 2.2e7 (×5.8) | −6.5e-6 |
| 2.150 | 0.143 | 9.08 | 1.5e8 (×7.0) | −2.6e-4 |

- max|B| **triples in one 25 ms step** while the transport current rises
  smoothly 9.4 → 10.4 A (sigmoid).
- The event is localized: the **bottom two tapes** reorganize — tape 1's mean
  Jz drops 18 % while its transverse current explodes ×40; tape 2's mean Jz
  doubles. Tapes 3–8 move ~10 %.
- The tape-1 transverse current had been growing at 1.04×/frame for the whole
  preceding window, then switched to ~6–7×/frame: a per-step ratchet, not a
  single bad solve.
- The free axial cohomology generator current I_1 dips at 2.125 s and jumps
  ×40 at 2.150 s. I_0 stays clean.
- The solver reported **2-iteration Picard convergence throughout the onset**;
  the first timestep rejection comes ~300 ms later (t = 2.45 s, and it is the
  *thermal* watchdog). The Picard-2 exit-residual floor had degraded from
  −156 dB early in the run to −80 dB at onset (tolerance −70 dB; dB = 10·log10).
- Everything is deep sub-critical (max j/jc ≤ 0.15, T = 77 K everywhere), so
  ρ_HTS ~ 1e-30 Ω·m: the HTS stiffness block is numerically zero.

The clean comparison (coarse dir, same byte-identical mesh, riva law) reaches
max B = 9.1 mT only at t ≈ 2.70 s — the noisy run's post-jump state matches the
comparison's state ~0.55 s later on the magnetization curve.

## 2. Verdict on the three candidate causes

**Bad material model (sp-ap) — refuted as the trigger.** The sp-ap B-axis
starts at log10 B = −2, i.e. **10 mT**; the entire pre-jump run sits below it,
where `JcFunctionDatabase::eval` clamps (jc, n are B-independent, dB-tangent
exactly 0). At the clamped edge the angular profile at 77 K is smooth: jc
anisotropy 4 %, n ∈ [24.5, 26.7], log10-jc second differences ≤ 3e-5. The clamp
is *internally consistent but physically wrong for a mT self-field deck* — a
real hygiene item (§5.3), not this event's cause.

**Error in the piecewise law — refuted as the trigger.** In the entire visited
state space (J/jc ≤ 0.15, n ≈ 25, 77 K) the piecewise and riva laws agree to
~1e-28 relative; the piecewise-specific machinery (Bézier blend, n < 4
degeneracy, complex-sqrt branch) starts at J > 1.26·jc and was never reached.
Through the onset only Picard ran, and `h_picard` consumes ρ only — no
tangents (`mt_maxwell_h.cpp:117-153`; promotion to Newton needs ε < 1e-4 *and*
iteration > 1, which a 2-iterate exit never satisfies). Both auditors
independently searched for an evaluation site where the laws split at these
states (thin-shell, side connectors, ghost) and found none — ghost consumes
saved `element_rho`, side connectors are metal-only by assertion. The known
value/tangent mismatch in `drho_piecewise_dJ` below j1 is real but Newton-only
(law debt list). Caveat kept open: nodal Exodus extrema do not *bound*
Gauss-point states; an instrumented run should confirm (§6).

The riva-vs-piecewise correlation itself is **confounded**: the clean long run
used a different build (c695af6c, `hphiTrun`, Blaze) than the noisy run
(efe46051, `belfem`, Armadillo); its deck is not preserved; the two present-day
decks also differ in `algorithm` (edited to Picard at 02:38 during this
session) and in `save every` — and save cadence is a *numerical* discriminator,
because the controller clamps Δt to land on save points
(`cl_FEM_Controller.cpp:2536-2559`). Same-build A/Bs from earlier that evening
(19:48 build) failed for **both** laws with different decks.

**Undiscovered defect / numerical — confirmed as the mechanism family**, in a
sharper form than the brief proposed. Three interacting elements:

### 2a. The magnetic system carries ~9 unpinned gauge constants (measured)

Prompted by Christian's mid-round question ("could the buffer layer be
unpinned?"): the magnesia buffer layer of each tape is pure φ
(`DomainType::Buffer`; ghost facets deliberately skip Buffer interfaces,
`cl_ThinShellFactory.cpp:2118`). Probing the run's own output:

- Each tape's buffer φ is a **single per-patch constant** (tape 1, frame 84:
  80768.5 ± 0.1 over the whole patch — the ±0.1 is the physics). Eight patches,
  eight **independent** constants.
- The constants are **redrawn every solve** (sign flips included); their
  envelope grew from 4.9e-3 at t = 0.05 s through 5.2e2 at 1.25 s to 6.1e4 at
  2.10 s and −3.9e7 at the onset frame.
- **The air constant floats too:** φ at the bearing corner node equals the
  whole far-air offset and wanders between −1.9e6 and +3.4e4 across frames;
  *no* air node holds φ = 0 across frames. The dof table prints exactly one
  prescribed magnetic dof — whatever it pins, the written air φ level is not
  pinned. Why the deck's `bearing { nodes : 9 }` is ineffective is an **open
  code question** (§6).
- The h-field interface coupling pins the buffer φ *gradient* (the ±0.1 rides
  along correctly) but can never pin the constant — Christian's prior
  assumption that connection to pinned air suffices is refuted by measurement.

Nine floating constants make the matrix exactly singular every step; MUMPS
null-pivot handling explains the 1e17–1e18 condest baseline from step 1
(spikes to 1.9e21 are intermittent and correlate with, but do not date, the
onset — a 2.0e19 spike at t = 1.85 s sits on a healthy step).

### 2b. The physical near-null modes

With ρ_HTS ~ 1e-30, the inter-tape current split and the axial generator are
resistively unpinned — determined only by inductance against a
1e-7-tolerance iteration. This is the mode family that actually moved (tape
1/2 reorganization, I_1 excursion).

### 2c. The ratchet: always-accept Picard semantics (both auditors, verified)

The reported Picard residual is the **pre-update entry-state residual**: the
solver computes r = A(x_k)x_k − b, solves, applies the (relaxed or
Anderson-mixed) update, then restores the saved r for reporting
(`cl_FEM_DofMgr_SolverData.cpp`, comment: "Picard reports the PRE-update
residual ( always-accept semantics, Messe 2023 §4 )"). The controller commits
the update unconditionally (`cl_FEM_Controller.cpp:1126-1142`) and
`run_coupled` exits on relative ε ≤ 1e-7 with the absolute escape inert at its
0.0 default and `min iterations = 2`. So "converged in 2 Picard iterations"
certifies iterate x₁ — **the state actually written and time-stepped is x₂,
whose residual is never evaluated.** The watchdog is structurally blind at a
2-iterate exit (spared while iteration ≤ min). Anderson (magnetic depth 3,
γ-cap 100) can amplify a weakly-observed direction but at a 2-iterate exit has
one history column — a contributor, not the window-8 amplifier the brief
guessed.

**Important nuance:** this is *documented design*, adopted deliberately after
the greg3 flux-front freeze (the same comment block explains why a Picard trial
cannot be judged within one frozen assembly). Its safety assumption — the
Picard map is contractive, so x₂ is at least as good as x₁ — fails on a system
whose null directions cost nothing. Changing it is a design decision, not a
bug fix.

**Mechanism summary:** an exactly-singular, resistively-unpinned system is
advanced by an iteration that never checks the state it commits. The unchecked
second-map increment ratchets across steps (1.04×/frame in the quiet phase);
at t ≈ 2.125 s the compounding crossed into the regime where the committed
increments reorganize real currents, and the system jumped to a differently
magnetized state ~0.55 s ahead on its own curve. Messe et al. 2023 **§4**
(Eqs. 10–13 — not §2.7 as CLAUDE.md's routing table says) prescribes ε = 1e-11
precisely because loose tolerances let "current density per element oscillate…
between ~0 and ±2·jc" in this material class.

## 3. Why piecewise-vs-riva *appeared* causal

The comparison runs differ in build, backend (Blaze/Armadillo), executable,
algorithm keys, BDF order, and the save cadence that shapes the actual Δt
sequence. On a system living this close to singular, any of these moves the
onset time. The law is the *least* likely discriminator of the set (agreement
to 1e-28 in the visited regime). Grok's discriminator ranking for the A/B:
backend first, BDF order, Δt/save cadence, algorithm — resistivity type last.

## 4. Cleared / minor findings

- `EF` thin-shell rho wiring is symmetric across laws (dispatch,
  `cl_FEM_Calculator.cpp:267-331`); Picard/Newton kernel split verified.
- Ghost + Coulomb-gauge chi are **active by default** (`chi = 1e-4` since
  2026-08-27) although the deck's penalty blocks are commented out — the deck
  comments mislead; common-mode across both runs, so not a discriminator, but
  worth an explicit `chi : 0` if "off" was intended.
- Value/tangent mismatch of `rho_piecewise`/`drho_piecewise_dJ` below j1
  (raw vs parallel-combined): real, inert here, stays on the law debt list.
- The 2026-08-27 controller diff (thermal budget, DR-115 rank guard) does not
  fire at a 2-iterate exit — comparison confound only.
- The coarse 4.3 s "suspect" band was a different, already-fixed defect
  (DR-125, reset_timestep dof re-seed).

## 5. Recommended mitigation (ranked)

Physics/formulation decisions are Christian's; items marked ⚑ need his ruling.

1. **Clean A/B first (executable gate):** one build, one deck, 25 ms cadence,
   vary exactly one key per run: piecewise↔riva, Anderson on↔off,
   algorithm Newton↔Picard, Blaze↔Armadillo. Without this, every causal story
   stays unfalsifiable. (The coarse deck edited at 02:38 to Picard is the
   first leg of exactly this.)
2. ⚑ **Pin the gauge.** Repair/verify the air bearing (today the written air φ
   level floats despite `fixed = 1`), and anchor one φ dof per enclosed buffer
   patch (8 more). This removes the exact singularity, should drop the condest
   baseline by orders, and is independent of any physics choice — the
   constants are pure gauge. The free axial generator I₁ is legitimate and
   **stays free**.
3. ⚑ **Guard the accepted state.** Either evaluate the post-update residual
   before accepting a Picard/Anderson exit (one extra residual assembly per
   timestep), or raise `min iterations` / add an `absolute tolerance` so a
   2-iterate roundoff-floor exit cannot advance time. Interacts with the
   documented always-accept semantics and the greg3 history — design round
   required.
4. **Tighten the magnetic tolerance toward 1e-10/1e-11 only together with (3)
   or an absolute escape** — the Picard-2 floor is already −80 dB at onset and
   −69 dB at the first reject; a bare 1e-11 target with no floor escape will
   cut Δt forever.
5. ⚑ **Extend the sp-ap B axis below 10 mT** from defensible data or a
   documented low-field model — the whole operating range of this deck
   currently has no measured B-dependence, and post-jump the field sits at the
   clamp edge (kink in ∂ρ/∂B once Newton engages). Hygiene and lock-in
   insurance, not the onset fix.
6. **Doc fixes:** CLAUDE.md literature routing "§2.7" → "§4" for the
   checkerboarding/tolerance citation; deck-comment trap for the
   default-on chi/ghost penalties.

## 6. Open questions / follow-up gates

- ~~What does the single prescribed magnetic dof actually pin, and why does the
  written air φ still float?~~ **ANSWERED 2026-08-28 (parallel Fable session,
  reported by Christian): the air bearing was placed on a node lying in a
  PERIODIC PLANE, and the periodicity constraint eliminates/replaces that dof —
  the pin was silently consumed.** Christian's assessment: this and the
  buffer-island pinning are both straightforward fixes. Consequences for the
  automatic pinning are recorded in §8.3 A7.
- Instrument Gauss-point maxima (J/jc, B) and the accepted-state residual
  (including the I₁ component) for one window around 2.1 s — closes the
  "off-envelope evaluation" caveat and measures the ratchet directly.
- Does the fixed spatial pattern of the buffer amplitudes (envelope growth
  5e-3 → 4e7) correlate with the transport current level? (Cheap: re-run the
  frame probe on the surviving riva keyframes.)
- Register rows filed this session: DR-126 (unpinned φ constants / bearing),
  DR-127 (accepted-state guard ruling), DR-128 (sp-ap B-axis floor).

## 7. Provenance

| artifact | role |
|---|---|
| `cmake-build-debug/tapestack3d/out.txt`, `iv_results.csv`, keyframes 76–94 | noisy run (piecewise, efe46051, Armadillo, `belfem`) |
| `cmake-build-debug/tapestack3d_coarse/` keyframes 38–57 | clean comparison (riva era, c695af6c, Blaze, `hphiTrun`; deck not preserved) |
| `share/material/sp-ap.hdf5` (rebuilt 2026-08-27) | jc/n table, B-origin 10 mT |
| `tmp/ai_exchange/review_tapestack3d_jjc_noise.md` | frozen pre-registration, both audits, verification, reconciliation |
| Jury | Codex (medium ~80 %) and Grok (high on controller path) — convergent on §2c; every cited line re-verified this session |

---

## Addendum 2026-08-28 (night session, independent corroboration of §2a)

Added by the concurrent quench-front session, which reached §2a's hypothesis from the
code side (buffer = `DomainType::Buffer` -> phi kernels via the rho-less test in
`ThinShellFactory::create_buffers:2057`; the h-phi interface constrains only grad-phi, so a
patch constant is invisible to it) and then measured the effect of removing the layer.

**Matched-time condest, with vs without the buffer layer** (`tapestack3d` piecewise vs
`tapestack3d_nobuffer` riva, both from t = 0, byte-identical mesh, both BDF1; at t < 400 ms
the tape is deeply sub-critical — I(400 ms) = 0.42 A — so the two laws agree to ~1e-8 and the
buffer is the operative difference):

| t | with buffer | without buffer |
|---|---|---|
| 5 ms | 2.6e17 | **4.2e16** |
| cruise to 400 ms | climbs to 1.8e18, spikes 5.6e18 | **steady 1.8e17 - 3.4e17, no spikes** |

Removing 8 of the 9 floating constants drops the condest baseline by about an order of
magnitude and removes the spikes — consistent with §2a's reading that the unpinned constants
set that baseline, and it isolates the buffer patches as the dominant contributor while the
air-bearing defect (the 9th constant) remains. **This is a diagnostic, not a fix:** the buffer
is what electrically severs the tape halves, so the no-buffer deck is different physics. It
supports mitigation 2 (pin one phi dof per buffer patch + repair the air bearing) over any
deck-side workaround.

**Also observed in the coarse-dir 4.35 s band, same night, post-DR-125:** three independent
runs stalled with the magnetic Picard residual plateauing at -67.8 / -67.9 / -69.4 dB against
the deck's -70 dB target, and a promoted Newton actively degraded a -69.41 dB iterate to
-57.5 dB with relaxation collapsing to 0.07 (out8), while the same state with promotion
disabled advanced at 50x the timestep (out9). A fixed near-singular subspace is the natural
reading: the direct solve puts arbitrary content in the null directions, which harms a
Newton direction while a Picard map never inverts them.

**Caution this report should be read back onto that observation (§2c):** at a 2-iterate exit
the Picard path never evaluates the state it commits, so "Picard-only advances further" is
partly structural blindness rather than demonstrated robustness. The Picard-only
configuration should not be adopted as a recommendation until mitigation 3 (accepted-state
guard) is settled.

## Addendum (2026-08-28, later): why the bearing on point 9 does not pin the air

Question (Christian): how can the air float if the deck sets `bearing { nodes : 9 }`?
Answer, as far as the tree and the run data settle it:

**The single "prescribed (fixed) : 1" dof is NOT the bearing — it is λ₀,** the
transport-current generator, fixed to I(t) every step by `IWG_Maxwell::set_currents`
(`cl_IWG_Maxwell.cpp:116`; master-only, consistent with its own comment and with I_0
tracking the sigmoid exactly). The bearing pins nothing:

- φ(node 9) = −0.117 A at the very first keyframe (t = 25 ms) and −3.4688 in the raw
  t = 0.825 s memdump (matching the exodus to 6 digits — the written phi IS the solver
  state, no postprocessing shift). A working Dirichlet pin would read exactly 0.
- The coarse run is identical: −0.43 at frame 1, −754 at 2.1 s, **−2.1e9 at t = 4.35 s** —
  i.e. the same defect across `belfem`/`hphiTrun`, both builds. The 4.35 s coarse wall
  band was fought on top of φ-garbage nine orders above the physical scale.
- The raw solver φ field grows monotonically garbage-ward: range [−13.7, +57.2] A at
  0.825 s → [−1.6e7, 2.2e9] at 2.48 s (fine-run memdumps).

**Every static link of the bearing chain was verified individually and reads correct:**
the BC is parsed with domains [9]; gmsh 4.1 point-element blocks map vertex 9 → node 9
(an outer air corner, a sensible bearing point); vertex 9 survives enrichment
(`tapestack3d.bfm` `vertices/ids` contains 9); `Kernel::create_field` runs
`DofManager::set_equation`, which creates bearings and links them to the node's dofs
BEFORE `MaxwellFactory` calls `impose_dirichlet(0.0)` (`cl_MaxwellFactory.cpp:811`);
nothing frees the fix afterwards; `init_dof_values` respects `is_fixed`; the free/fixed
graph split happens on the master after all of this. Yet the composite provably fails —
which points at the two SILENT failure modes designed into the chain:

1. `Bearing::impose_dirichlet` returns quietly when the bearing carries no dofs
   (`cl_FEM_Bearing.cpp:90` — "empty bearing contains no dofs").
2. `BearingData::bearing(id)` returns the shared EmptyBearing for any unknown id
   (`cl_FEM_DofMgr_BearingData.hpp:97`) — a misses-the-map lookup is indistinguishable
   from success at the call site.

Plus one fragile, undocumented dependency: bearing linking needs the node→dof
back-links, which `connect_dofs_to_mesh` builds in exactly ONE place — inside the
hanging-dof T-matrix builder (`cl_FEM_DofMgr_DofData.cpp:3492`). A deck without hanging
dofs can never link any bearing. This deck HAS 145 141 hanging dofs, so the specific
break here is still open; static analysis is exhausted.

**Discriminating experiment (minutes, serial):** boot the same deck under gdb and watch
the impose happen — the factory phase suffices, kill after the first timestep prints:

    cd cmake-build-debug/tapestack3d
    gdb -q ../bin/belfem
      (gdb) break belfem::fem::Bearing::impose_dirichlet
      (gdb) run
      (gdb) print mID          # expect 9
      (gdb) print mNumDofs     # 0 = the silent no-op, mystery solved
      (gdb) print mNode->id()  # if >0: which node, then `print mDOFs[0]->id()`

If `mNumDofs` is 0 (or the breakpoint never fires), the failing link is upstream
(EmptyBearing lookup or empty link table); if it fires with a dof and fixes it, the
loss is downstream (graph split / a later re-creation) and the next breakpoint is
`DofData::reorder_dofs`. Either way: after the root cause, the fix wants a loud
failure mode — an empty bearing named in the deck should be a `BELFEM_ERROR`, not a
silent skip — plus the per-buffer-patch anchors of mitigation 2.

---

# 8. Implementation manuscript — automatic gauge pinning of the φ islands

**Author of the design:** Christian, 2026-08-28 (end of the night session), dictated as final
thoughts and written up here.
**Author of the machinery survey and audit:** the concurrent quench-front session, same night.
**Status:** design + audit only. No source modified. Sequenced behind the open bearing-imposition
question (DR-126, under investigation in a parallel session).

## 8.1 The design, as stated

> To pin the buffer nodes, we must first sweep over all blocks of the element and tag the nodes
> that are connected to a buffer block. Since these blocks are not necessarily connected, we must
> create a graph that contains all nodes that sit on buffers, and do a BFS to identify how many
> non-connected subgraphs we have. From there, we always fix the first node. We have the
> machinery for that in BELFEM.

In steps:

1. sweep the blocks, select those tagged `DomainType::Buffer`;
2. tag the nodes belonging to those blocks;
3. build a graph over the tagged nodes, edges = "share a buffer element";
4. traverse (BFS/DFS) to count the connected components — the φ islands;
5. fix one node per island.

The premise is correct and matches §2a's measurement: each tape's buffer φ is one independent
constant, and the h-φ interface pins only the gradient, never the constant.

## 8.2 The machinery that already exists

| need | what to use | where |
|---|---|---|
| connected components ("how many subgraphs") | `graph::find_connected_partitions( Graph & )` — runs `dfs`, returns the count, writes the partition id per vertex, orders partitions by size (0 = largest) | `src/math/graph/fn_Graph_find_connected_partitions.{hpp,cpp}` |
| traversal primitives | `graph::bfs( Graph &, Vertex * )`, `graph::bfs( Graph & )`, `graph::dfs( Graph & )` | `src/math/graph/fn_Graph_{bfs,dfs}.{hpp,cpp}` |
| the graph type | `typedef Cell< graph::Vertex * > Graph` | `cl_Graph_Vertex.hpp:329` |
| nodes are already graph vertices | `class Node : public Vertex` (mesh) and `class Vertex` (graph) | `src/mesh/cl_Node.hpp:29`, `src/math/graph/cl_Graph_Vertex.hpp:31` |
| node → element adjacency | `mesh::Vertex::number_of_elements()` / `add_element()` accessors. **Precondition:** the container is filled by `Mesh_ConnectivityCalculator` (`cl_Mesh_ConnectivityCalculator.cpp:81`, `tElement->node(k)->add_element(tElement)`), not unconditionally at mesh load — meshes loaded with connectivities disabled do not have it. Assert it is populated at the point of use, or flood-fill via `Element::node(k)` and a locally built node→element map instead | `src/mesh/cl_Vertex.hpp:204-208`, `cl_Mesh_ConnectivityCalculator.cpp:81` |
| buffer blocks already identified and walked | `ThinShellFactory::create_buffers()` — tags `DomainType::Buffer` by the **rho-less test** `if ( material->have( MaterialProperty::rho ) ) continue ;` and already loops the blocks' elements, edges and faces with flag hygiene | `src/fem/kernel/cl_ThinShellFactory.cpp:2030-2070` |
| marking / dedup without a map | the 8-slot entity flag bitset (`flag()`, `unflag()`, `is_flagged()`), house rule: clear before AND after the pass | `cl_Graph_Vertex.hpp`, used e.g. at `cl_ThinShellFactory.cpp:2039-2048` |
| the pin itself (existing, deck-facing) | `BoundaryConditionType::Bearing` — "a bare node constraint with no evaluated scalar"; deck path `bearing { nodes : ... }` parsed at `cl_MaxwellBoundaryConditionFactory.cpp:67-73` (shared with `Gauge`) | as cited |
| the pin's dof-side machinery | `class Bearing` (`node()`, `allocate_dof_container()`, `insert_dof()`) and `DofManager`'s `create_bearings()` | `src/fem/kernel/cl_FEM_Bearing.hpp`, `cl_FEM_DofMgr_BearingData.{hpp,cpp}` |

## 8.3 Audit — one blocking trap, and five design corrections

### A1 (BLOCKING). `find_connected_partitions` must NOT be called on live mesh nodes

It is the right algorithm but it writes into fields that mean something else on a `mesh::Node`:

- `find_connected_partitions` does `tVertex->set_index( tIndex++ )` — on a mesh node, `index()`
  is the **dense index used for dof and field-data mapping**;
- `dfs` does `set_owner( gNoOwner )` and then `set_owner( aNumSubGraphs++ )`
  (`fn_Graph_dfs.cpp:32,44`) — on a mesh node, `owner()` is the **MPI rank**, tested as
  `tNode->owner() == mCommRank` in `cl_FEM_DofMgr_FieldData.cpp:111,124,256`;
- it also overwrites `set_level()` and leaves flags set.

Calling it on the mesh's own nodes would silently destroy the parallel ownership map and the
index space, and the symptom would look like a distribution bug far from its cause. Three ways
out, in preference order:

**CHRISTIAN'S RULING (2026-08-28): mirror graph — "tag the nodes that sit on the buffer layers,
and create a duplicate graph. We have done this before."** This is the right call: it reuses
`find_connected_partitions` unchanged, and the throwaway vertices absorb the `set_index` /
`set_owner` writes harmlessly. **The idiom already exists in the tree, twice**, and a patch
should copy it rather than invent one:

- `CutFactory` (`src/homology/cl_CutFactory.cpp:1066-1100`) — the closest analogue: collect the
  flagged entities, `Graph tGraph( tCount, nullptr )`, one `new graph::Vertex()` each with
  **`set_id( tFace->id() )`** (carries the mesh id for mapping back) and
  **`set_index( tCount++ )`** (position in the collection, so `tFaces( tVertex->index() )`
  recovers the mesh entity), then wire adjacency through the shared-entity relation with
  temporary unflagging to exclude self-connection.
- `Mesh` (`src/mesh/cl_Mesh.cpp:2762-2781`) — same shape for facets, with the comment "one graph
  vertex per unique facet key, in sorted key order, so that `find_index_in_unique_cell()`
  replaces a map".

Applied here: flag the nodes of every `DomainType::Buffer` block (plus Air, per A2), collect
them into a `Cell< Node * >`, build the mirror `Graph` carrying `set_id( node->id() )` and
`set_index( position )`, wire two mirror vertices together when their nodes share a buffer
element, call `graph::find_connected_partitions( tGraph )`, read the component id off each
mirror vertex's `owner()`, choose each component's representative by **smallest `id()`** (A3),
then delete the mirror vertices. The mesh's own nodes are never written to — only flagged, and
the flags are cleared before and after.

Alternatives, recorded but not chosen: a hand-rolled flood fill over the buffer elements using
only the flag bitset (no allocation, but reimplements what the graph module already provides);
or save/restore of `index`/`owner`/`level` around a direct call — rejected outright, since it
breaks silently the day another field is added to `dfs`.

### A2. Pin per φ-component, not per buffer patch — and only if the component has no pin yet

A buffer patch that *does* touch the air φ region is already covered by the air's bearing.
Adding a second pin inside one component over-constrains it: it forces the potential difference
between two points to zero, which is wrong physics, not extra safety. The correct rule:

> For each connected component of the φ-carrying node set, if the component contains no
> already-prescribed dof, fix its representative node.

Running the component analysis over **all φ blocks (Air + Buffer)** rather than Buffer alone
gets this for free, and it also repairs the air region itself if its own bearing is ineffective
(§2a, DR-126).

### A3. "Fix the first node" must be a global rule, not an iteration-order rule

Under MPI each rank enumerates a subset in its own order, and `find_connected_partitions`
orders partitions by size, not by identity. The representative must be chosen by a
rank-invariant key — **the globally smallest node id in the component**. Otherwise different
ranks pin different nodes and the constraint set is inconsistent.

### A4. Sequencing: the imposition path must be proven before the automatic pins are added

§2a measured that the air φ level floats **despite** the deck's `bearing { nodes : 9 }` and
despite exactly one prescribed magnetic dof appearing in the dof table. If the imposition path
is itself broken (DR-126, under investigation), adding eight more pins will change nothing and
the negative result would wrongly discredit this design. **Verify the single existing bearing
bites first, then land the automatic pinning.**

### A5. Where the code belongs

`create_buffers()` runs during thin-shell construction, before the dof manager exists, while the
`Bearing` machinery lives in the dof manager. The clean split is: **identify components and
representative node ids** in the factory (it already owns the Buffer tagging), then **feed those
ids into the same path the deck's `bearing` uses**, so an automatic pin and a user-written one
are indistinguishable downstream — one code path, one behaviour, one place to debug.

### A7. Never pin on a periodic plane — and verify the pin survived

**This is the defect that produced the floating air constant** (§6): the deck's
`bearing { nodes : 9 }` sat on a node in a periodic plane, and the periodicity constraint
eliminated the very dof the bearing meant to prescribe. The pin was consumed silently, while
the dof table still reported one prescribed magnetic dof.

Two consequences for the automatic pinning of §8.1:

1. **A3's representative rule must exclude periodic-partner nodes.** "Smallest id in the
   component" is rank-invariant but blind: if that node lies in a periodic plane the auto-pin
   evaporates exactly as the manual one did — eight times over, and silently. Rule: smallest id
   **among nodes that are not periodic partners**; a component with no such node is a hard
   error worth reporting, not a case to skip quietly.
2. **Add a post-condition check.** After the dof manager is built, assert that every intended
   pin is present and prescribed in the final dof table. This whole failure class is a
   constraint that is requested, reported as applied, then quietly removed by a later stage;
   only an end-of-pipeline check catches it. Setup path, so `BELFEM_ERROR` per the error-tier
   rule — and it would have caught the original defect on its first run.

### A6. What must NOT be pinned

The free axial cohomology generator I₁ is legitimate physics and stays free (§5.2). This pass
touches φ **node** dofs only — never cut or generator dofs. Note also that thin-shell layer
nodes are duplicates of the surface nodes: geometric coincidence with an air node is not
connectivity, which is precisely why the patches are islands.

## 8.4 Verification gates for the eventual patch

1. **Count gate:** the pass reports the number of φ components found. Expect 8 buffer patches
   plus the air component on this deck; any other number is itself a finding.
2. **Conditioning gate:** the condest baseline should fall from the measured 1e17–1e18 to a
   non-singular value. Reference points already measured (addendum above): removing the buffer
   layer entirely — i.e. deleting 8 of the 9 constants — moved the baseline from 2.6e17 to
   4.2e16 at t = 5 ms and from a spiky 1.8e18–5.6e18 to a steady 1.8e17–3.4e17 in cruise.
   Pinning should achieve at least that **without** changing the physics.
3. **Constant gate:** the per-patch φ constants (measured 80768.5 ± 0.1 on tape 1, envelope
   growing to −3.9e7) must collapse to ≈ 0 ± the physical ripple, and must stop being redrawn
   every solve.
4. **Physics no-change gate:** with the pins in place, j/jc and B fields at matched keyframes
   must be unchanged within tolerance against a pre-patch run — the constants are pure gauge, so
   a physics change would mean the pin is wrong (A2 over-constraint being the first suspect).
5. `make check-fast` green; a `chi : 0` deck unchanged.

## Addendum 2 (2026-08-28, night): ROOT CAUSE FOUND AND VERIFIED — the bearing node is periodically condensed

Christian authorized autonomous runs; the gdb experiment and a follow-up reproducer ran in
`cmake-build-debug/bearing_probe/` (isolated copy of the deck; originals untouched).

**gdb result (serial boot of the unmodified deck):** the bearing chain WORKS at the dof level —
`impose_dirichlet` fires from `cl_MaxwellFactory.cpp:811` with the real bearing (mID = 9,
mNumDofs = 1), fixes the φ dof of node 9 (field index 104905 — the exact node probed in the
run data), and the flag survives to the first timestep (`mFixedFlag = true`,
`mDirichletValue = 0`). But the dof it pins carries `mIndex = gNoIndex` and
`mNumberOfSources = 1`: **it is a HANGING (condensed) dof.**

**Why node 9 hangs:** it sits at (0, −0.05, z = L) — on the **periodic target face**. The
deck's `topology { periodic { source : 11,13,53 ; target : 14,16,56 } }` condenses every
z = L node onto its z = 0 partner (weight 1). Mesh check: the hanging-node list contains
exactly the z = L twins of the low-id corner points — {7, 8, 9, 10, 14, 15, 16, 20} hanging,
their z = 0 partners {2, 3, 4, 5, 11, 12, 13, 17, 18, 19} free. A Dirichlet flag on a
condensed dof never reaches the assembled system: the free/fixed split skips hanging dofs
(they are the third, "condensed" category), and the T-matrix overwrites the value from the
source. The pin evaporates — **silently**, because `Bearing::impose_dirichlet` accepts a
hanging dof without complaint.

So: every individual link is correct; the defect is the *composition* bearing × periodicity.
The deck author could not have seen it — node 9 is a perfectly sensible bearing point on a
non-periodic deck, and nothing warns when periodicity demotes it.

**Reproducer (bearing moved to node 4, the z = 0 partner; 30 ms probe run, serial):** the
magnetic dof table flips from `prescribed (fixed) : 1` (= λ₀ only, the transport-current
generator — NOT the bearing, as established in Addendum 1) to **`prescribed (fixed) : 2`**
with the free count down by exactly one. φ-field and conditioning checks of the probe run
follow below.

**Fixes, in order:**
1. **Deck (immediate):** `bearing { nodes : 4 ; }` — any geometry point NOT on the periodic
   target face. This anchors the air constant. The 8 buffer-patch constants remain
   (DR-126 half b — per-patch anchors, Christian's call).
2. **Code (small, high value):** `Bearing::impose_dirichlet` must refuse silently-inert
   targets: `BELFEM_ERROR` when the dof is hanging (or, more graceful: chase the weight-1
   single-source chain and fix the SOURCE dof — mathematically equivalent for periodic
   pairs), and `BELFEM_ERROR` when a deck-named bearing resolves to the EmptyBearing or to
   zero dofs. Three silent no-ops become three loud errors.
3. The `prescribed (fixed)` count printed at init should list WHICH dofs (field + node id)
   when the count is small — one line that would have shortened this hunt by hours.

**Reproducer completed (verified by execution, keyframe t = 25 ms of the node-4 probe run):**
φ(node 4) = 0.000000 exactly (the pin holds), φ(node 9) = 0 as well (its condensed dof
follows the pinned source through the weight-1 periodic pair — the condensation mechanism
confirming itself), and the far-air φ is at physical scale (mean 4.8e-2 A, max 0.134 A) with
no floating constant. Step-1..4 conditioning stays 3.5e17–7.2e18: the 8 buffer-patch
constants still dominate the near-singularity, as §2a and the no-buffer diagnostic predicted.
The air-bearing half of DR-126 is closed at reproducer level; the buffer-anchor half and the
loud-failure code fix remain open.

---

## Addendum 3 (2026-08-29): INVESTIGATION CLOSED — the noise is gone on the fixed build

**Evidence level: run-verified.** The deck that filed this report
(`cmake-build-debug/tapestack3d`, sp-ap table, `resistivity type : piecewise`, same
`input.conf`) was rerun on a build carrying both fixes this investigation produced —
DR-126 (automatic gauge pinning of the φ islands, §8) and DR-127 (certified exit,
`devlog/dl20260829_dr127_certified_exit.md`). 137 timesteps, t = 0 → 3.3 s, past the
2.55 s that the Riva comparison deck had reached.

**The 2.125 s frame is no longer distinguishable from its neighbors.** Across
2.00–2.40 s, `max|B|` and `max j/jc` both grow by ×1.047 per frame — the same constant
factor at every single frame, tracking the transport-current ramp:

| t [s] | max\|B\| [T] | max j/jc | frame ratio |
|---|---|---|---|
| 2.100 | 3.2181e-03 | 4.5098e-02 | ×1.047 |
| **2.125** | **3.3690e-03** | **4.7212e-02** | **×1.047** |
| 2.150 | 3.5269e-03 | 4.9426e-02 | ×1.047 |

Measured defect signature for comparison (§1): `max|B|` ×3 **in a single step** at this
frame, tape-1 transverse current 5.8×/7.0× per frame, I₁ ×40. On the fixed build the
largest frame-to-frame I₁ ratio *anywhere in the run* is ×1.115 (at t = 2.25 s); I₁
grows ×8.0 smoothly over 1.17 s while I₀ ramps 9.5 → 17 A. The scatter measure
`std(j/jc)` over active nodes grows at the same ×1.047 and the active-node count climbs
monotonically 19306 → 19486 — a flux front, not noise.

Controller state through the window: every step exits certified at 2.6e-15…5.0e-15 after
exactly 2 Picard bodies; zero timestep cuts, resets, divergence strikes or stalls in 137
steps; conditioning flat at 1.26e9 (it was 1e17–1e21 when this report was written).

**Attribution.** The two fixes shipped together and this single run cannot separate them.
What the run does establish: with the gauge pinned and every exit measured, the deck is
smooth well past the onset. What it does *not* establish is that DR-127 alone would rescue
a near-null system — DR-126 removed the singularity that made the always-accept
contractivity assumption fail in the first place. Both remain justified on their own
arguments (§2a, §2b).

**§2c (material model) stands as written and is NOT cleared by this run.** DR-128 is live
in these very frames: `max|B|` stays between 2.7 mT and 5.6 mT across the window, entirely
below the sp-ap table's B-axis floor at 10 mT, so jc and n are B-independent (clamped)
throughout. The noise disappeared with that clamp still in place — which is evidence the
clamp was not the trigger, and no evidence at all that the clamp is acceptable. DR-128
stays open on its own merits.

**Remaining gate for DR-127** (not reachable from this deck): a segregated
`hphirun`/`hphiTrun` run that hits a thermal sub-step cut, the only path exercising the
`solve_thermal()` / `mResetThermal` return that the code-audit jury corrected.

**Status: this investigation is closed.** The reported symptom is resolved and reproduced
as resolved. Residual work lives in the register rows (DR-127 segregated gate, DR-128) and
no longer in this document.
