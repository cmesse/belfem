# The register sweep, the DR-18x cluster, and the DR-168 second opinion

**Date:** 2026-09-24
**Purpose:** Christian opened with "let's have a look at the debt register": three rows were done
in the tree but still open, two rows carried "ruling in progress" from a quiz whose answers had
survived but whose follow-ups had not run. One quiz settled eight rulings; a second session
landed five of them as the DR-18x cluster (DR-184, DR-182, DR-186, DR-187, DR-166), each its own
diff and gate, and ran the physics round DR-168 had been waiting for since 2026-09-22.
**Module:** `src/mesh` (`cl_Element_QUAD16`, `cl_Mesh_PeriodicityFactory`), `src/fem/maxwell`
(`cl_MaxwellFactory` thin-shell neighbour admission), `src/fem/thermal` (`cl_ThermalFactory`),
`src/math/graph` (`cl_Graph_Vertex`), `tests/mesh`, `tests/fem`; `todo/` registers
**AIs involved:** Claude Fable 5.1 (read, probe, plan, implementation, gates, reconciliation).
Plan jury and code jury: Codex `gpt-6-astra`/high + Grok `grok-4.7`/high; round 2 on the DR-186
hunks at astra/xhigh + grok-4.7/xhigh. DR-168 physics round: Codex `gpt-5.6-terra`/high + Grok
`grok-4.7`/high, blind.
**Status:** cluster landed in the working tree behind three jury rounds; DR-168 ruled; commit and
the shared-tree `make check` owed
**Verification:** executable — scratch builds against a patched copy of the Debug tree's library
(`tests/mesh` facet-table oracle 8/8 and periodicity 10/10, thermal admission 2/2; the QUAD16 and
the periodic refusal RED on the pristine library); four DR-186 deck legs (two positive decks clean,
the B-H-conductor deck refused in the NDEBUG and the assert-active build alike); syntax gate on
every changed unit under both trees' flags. Reviewed, not verified: the parallel shape of the
DR-186 lookup (rank 0's dof manager holds only its own blocks); `make check` in the shared trees.
**Exchange:** `tmp/ai_exchange/review_dr18x_cluster.md` (+ `.diff`, `_r2.diff`),
`tmp/ai_exchange/review_dr168_ghost_operator.md` (+ `dr168_ghost_operator_brief.md`)

---

## 1. The register sweep

Currentness against the tree, then the quiz. Three rows were closed in code and still open:
DR-185 (`25dbc53`, Tier 2 slice 1), DR-160 (ruled O9 "warn" in the Tier 2 quiz, landed `9cb58e3`),
DR-188 (fixed 2026-09-23, closed in status, never struck). Christian struck all three. DR-168 and
DR-172 read "ruling in progress" because the 2026-09-22 quiz answers were recorded in the closed
Tier 1 plan but the second opinion Christian asked for on DR-168 never ran and DR-172's answer was
tentative. Eight rulings today, all by quiz: DR-168 physics round first; DR-187 narrowed refusal;
DR-182 refuse the overlap; DR-186 restore the guards; DR-172 narrow refusal confirmed; DR-166
keep the assert, fix the limit; DR-175 widened with the Tier 2 leaks; strikes as above. Open
population after the sweep: `[F]` DR-168, DR-171; `[P]` 14; `[W]` 3.

A peer's commit (`2f970b6`) swept this session's edits to `debt_register_closed.md` and
`test_gap_followups.md` along with its own todo tidy-up, so HEAD briefly had the three rows
struck in the archive and still listed in the open register; the open register's commit closes
that.

## 2. The cluster plan, and what the reads changed

Every claim in the plan was read from the tree, and DR-184's scope was executed: a scratchpad
oracle placed every Lagrange element type's nodes at their parameter coordinates and compared
each facet table with the convex hull of its corners and each edge table with the segment of its
first two nodes. Over 32 types, one defect — QUAD16 facet 1 (and edge 1, which forwards to it).
HEX64 is clean. That oracle became `tests/mesh/test_FacetTables.cpp`.

DR-187's ruling premise was wrong in both directions. The thermal kernel selects exactly
`Conductor`, `ThinShell`, `Ferro`, `Buffer`, and the dof manager builds blocks from that
selection, so the dispatch's `Air` arm was dead and its `default:` arm unreachable. The thermal
equation selects no sidesets at all, and the dof manager builds sideset groups only for a
non-empty selection, so the sideset loop iterated an empty container; its `Dirichlet` /
`Neumann` arms tested mesh types that nothing sets, and the enriched integration orders 5 / 13
it computed were never applied. Consequence recorded on the row: thermal sideset boundary
conditions are unsupported — a Dirichlet section dereferences the dof manager's null empty
sideset, Neumann is refused at parse. No shipped deck declares one. Christian ruled the dead
loop deleted and the gap recorded.

The plan jury (12 findings, all confirmed by source) corrected three things before code was
written: `mesh::Block` has no entity flag bitset (dedupe by a setup-scope `Set< id_t >`); a
builtin ferromagnet keeps a constant μ₀ unless a B-H curve is loaded, so the negative deck must
use `curve : RoxieIron`; and a third `uint32_t` limit sits in `init_vertex_container`.

## 3. What landed

| row | change | gate |
|---|---|---|
| DR-184 | `mNodes[ 8 ]` → `mNodes[ 7 ]` on QUAD16 facet 1 | `FacetTables.Quads` RED on the pristine library at facet 1 and edge 1, GREEN after; 8/8 over all families |
| DR-182 | `tag_periodic_sidesets` refuses a plane sideset whose type is neither `Default` nor `Periodic`, naming sideset, plane and type | three `Periodicity` cases on `PeriodicBrick` with a new lateral sideset; the refusal RED before, GREEN after |
| DR-166 | three limits via `decltype( mVertexCounter )`, one comment | syntax; no behaviour change |
| DR-187 | `ThermalFactory::block_is_supported` (public static) read by selection and dispatch; `default:` → `BELFEM_ERROR`; `Air` arm, dead sideset loop, `tOrder`, `mUseEnrichment` deleted | `ThermalBlockAdmission` over every enumerator 2/2 |
| DR-186 | see §4 | four deck legs |

Syntax gate: every changed unit, `-Werror` flags of both trees (the `build` tree is Debug today,
`cmake-build-debug` Release; the Release tree has no built library, so the executable gate ran
on the Debug tree's flags only).

## 4. DR-186: the code jury found the Lawrencium mechanism

The first landing restored both guards as Christian had ruled in the morning: a warning in every
build plus `BELFEM_ASSERT`. The deck legs passed (racetrack and tapestack3d clean; the racetrack
with an iron B-H conductor warned once per conductor block, 36 for 36, and continued under
NDEBUG; aborted on the assert with assertions active). Then the code jury (astra/high +
grok-4.7/high) found the P1 both the old commented-out error and the new code shared: the
material was read through `mMagneticField->block( id )->material()`, and the dof manager on rank
0 of a parallel run holds only the blocks rank 0 owns or borders — a miss returns the empty block,
whose material is null, so a valid copper block reads as "magnetic". That is the most plausible
"throws on Lawrencium without cause" (medium-high; the deck is not on record). Both auditors and
the coding philosophy also held that a setup-path wrong-answer check belongs on `BELFEM_ERROR`,
and that a rank-0-only assert leaves the other ranks running. Christian re-ruled: material by
label from the factory's own map, both predicates always active, the warning gone. The negative
deck now aborts with the named refusal in the NDEBUG unit; the positive decks are unchanged.

The other round-1 findings, all confirmed and applied: the facet test swallowed a
`std::runtime_error` from any edge index (now only edge 0 may throw); the thermal refusal text
said Air "has no thermal kernel" although `T_phi` serves Air; the lateral-sideset expectation
could pass vacuously (a node-count check added); and the O3 by-catch — naming QUAD9TS `quad9`
as the linear TS arms name their parents — was reverted, because it would corrupt the
edge-function factory's "implement EF_%s" message and the `element_type()` round trip. One
refuted: a mixed-type block judged by one element (a `mesh::Block` holds one element type).

Round 2 (astra/xhigh + grok-4.7/xhigh) on the DR-186 hunks alone: both approve. Codex's one nit — the comment called the ε test MaxwellData's predicate, which is exact equality — fixed. Grok's one surviving point: the μ refusal's single text could not tell a null material from a B-H curve or a μ function from an off-μ₀ constant, the undiagnosable shape the Lawrencium comment was avoiding; the message now names the cause, and the negative deck prints "its permeability is a B-H curve or a function, not a constant". Recorded as pre-existing residue: a sideset with two linear neighbour types of different edge counts, a `Ferro` neighbour that this check does not see while the postprocessor assumes μ₀ on both sides, and the order refusal being the `.bfm`-reload gate only.

## 5. DR-168: three derivations, one answer

The brief carried the kernel's four block lines and the two `D` operators. My own expansion:
`Dm·q = +∂z h_m`, `Ds·q = −∂z h_s`, so the coded consistency term is `ρ_harm (∂z h_m − ∂z h_s)`,
a slope jump, where the SWIP average is `(ρ_harm/2)(∂z h_m + ∂z h_s)`. For an exact solution with
continuous flux `φ = ρ_m ∂z h_m = ρ_s ∂z h_s` the coded term gives `2φ (ρ_s − ρ_m)/(ρ_m + ρ_s)`
where consistency needs `φ`: right only on the ray `ρ_s = 3ρ_m`, zero at equal ρ, and an O(1)
interface error that in-plane refinement cannot remove because α scales with the physical layer
thickness. Codex and Grok reached the same factor independently, Grok by Green's formula. No
reading makes the jump right: in SWIP the contrast lives in the flux weights and the penalty
scale, not in the sign pattern, and the 2026-03-19 devlog's own transcription of Ern eq. 2.16
(`δ⁺ρ⁺ = ρ⁺ρ⁻/(ρ⁺+ρ⁻)`) is that half weight. Ghost facets exist only between two conductive
layers of different materials, never beside a buffer layer, so the magnesia argument does not
reach this kernel. The penalty floor is settled by algebra: the code (`η k_reg`;
`2η(k_m + k_reg)` one-sided) is right, the module doc and the kernel comment are wrong. Grok
corrected the brief on one point: the existing `GhostConsistencyCoefficient` test does pin the
coded slave sign through the adjoint rows and will go red under the average, as its own comment
says. Christian ruled: adopt the SWIP average as the operator of record — one code round, its own
plan and jury, consumer inventory of the ghost-on decks first, since every ghost-on result
changes. Not verified: the volume side's curl scaling of the oracle (Grok's open risk), which the
manufactured-flux fixture sidesteps by asserting the kernel alone.

Christian struck all five rows at the close on the DR-42/49 pattern; the open population is `[F]` DR-168, DR-171; `[P]` 9; `[W]` 3.

## 6. Owed

- Christian: `make check` in both trees; the commit of the register, plan, devlog and the five
  diffs.
- DR-186 in parallel: a 2-rank run of a thin-shell deck beside a conductor would observe the
  fixed lookup; not run today.
- DR-168 code round (its own plan). Follow-ups plan R1, R9, R10 unchanged.
