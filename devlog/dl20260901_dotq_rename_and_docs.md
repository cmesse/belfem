# `heatloss` Renamed to `dotQ`, and Documented

**Date:** 2026-09-01
**Purpose:** Record the rename of the ohmic-dissipation mesh global from
`heatloss` to `dotQ`, and the documentation that now says what it actually is.
**Module:** `fem/maxwell`, `fem/kernel`

## Why

`7323d528` added a mesh global accumulating the ohmic dissipation of the
h-domain and wrote it to Exodus as `heatloss`. The name is ambiguous in the one
way that matters: it does not say whether the number is a **power** or an
**energy**. It is a power. An AC loss per cycle — the thing a user reading
"heatloss" would most likely assume they were getting — is its time integral,
which nothing in the tree takes.

`dotQ` fixes that at the name: a dot is a rate. It also matches the local
variable the kernels already used (`dotQ`) and the argument name (`aDotQ`).

## The rename

18 edits across five files, applied as count-asserted exact matches. Nothing in
`src/` says `heatloss` any more.

| From | To | Sites |
|---|---|---|
| `"heatloss"` (the Exodus global variable name) | `"dotQ"` | 7 |
| `'heatloss'` in the duplicate-global error message | `'dotQ'` | 1 |
| `save_heatloss()` | `save_dotQ()` | 6 |
| `Controller::reset_heatloss()` | `reset_dotQ()` | 5 |
| `Controller::collect_heatloss()` | `collect_dotQ()` | 5 |

The helper names went with the string deliberately. Renaming only the string
would have left `save_heatloss()` writing into a global called `dotQ`, which is
the exact ambiguity the rename exists to remove.

A units note was added at `save_dotQ()` in `mt_maxwell_h.hpp` — watts, not
energy, and the global is zeroed before every assembly so what reaches Exodus is
the value from that timestep's last assembly.

**This is a breaking change for any downstream script reading the Exodus global
by name.** Acceptable here: the field is one day old, unpushed, documented
nowhere, present in no deck and in no schema, and `grep` finds no consumer in
the repository.

## Name collision, examined and accepted

`dotQ` was already in use: `cl_IWG_StaticHeatConduction.cpp:63` declares
`mFluxFields = { "dotQ" }` and reads it at `:181` through
`mCalc->node_data( "dotQ" )`. That one is a **nodal field** on a static heat
conduction problem; the new one is a **mesh global** on a magnetic problem.

Not a technical collision. `Mesh` keeps the two in separate namespaces
(`field_exists()` vs `global_variable_exists()`), and Exodus does too —
`EX_NODAL` versus `EX_GLOBAL` in `cl_Mesh_ExodusWriter.cpp`. The two also live
on different meshes in a coupled run. Both quantities are genuinely heat rates,
so the symbol is honest in both places; they differ in scope and normalization,
not in kind. Recorded here because a reader meeting `dotQ` twice in BELFEM
should know it is deliberate.

## Documentation

New `§9.3 Ohmic Dissipation Global (dotQ)` in
`src/fem/maxwell/doc/maxwell_usage_guide.md`, and a pointer to it from the
H-kernel table in `src/fem/maxwell/doc/README.md`. It records:

- the definition and the **units — watts**, with the explicit statement that AC
  loss needs a downstream time integral that BELFEM does not perform
- the scope: one scalar for the whole mesh, no per-block decomposition
- the lifecycle table (created / zeroed per assembly / accumulated per element /
  MPI-summed) and the consequence that the stored value is from the timestep's
  **last** assembly, not an average over its iterations
- which five of the six h-kernels contribute, and why `h_ghost()` does not — its
  `rho` enters a Nitsche stabilization coefficient, not a dissipation term
- a note on the resistivity clamp, and why it does not bite (see the correction
  below); `rho_clamped()` records when it fires but is not propagated to `dotQ`

No input-contract update was needed: `dotQ` is an output, created
unconditionally, and reads no `input.conf` key. `check_doc_claims.py` 37/37.

## Codex language sweep

Run over §9.3 (`gpt-5.6-luna`, medium, slug `dotq_docs`). It returned a
readability pass and flagged two of my claims, **both correct and both applied**:

1. "set once in `cl_Communicator.cpp`" was too strong. The comment at
   `cl_Communicator.cpp:62-64` says executables *may* narrow the window. No
   executable in this tree does — but the design intent is that they can, so the
   text now says defaults are set in `Communicator::set_globals()`, that an
   executable may narrow them and none here does, and that a deck cannot.
2. "`gRhoMin` is zero, so the low side never bites" was wrong in the corner.
   `compute_rho()` does test `tRho < gRhoMin`; with `gRhoMin == 0` only a
   negative — nonphysical — resistivity trips it. Text corrected to say exactly
   that.

## The clamp claim, wrong twice and then measured

The first version of §9.3 said `dotQ` is a lower bound "wherever the power law
drives `rho` above `gRhoMax`". **That is false**, and the way it was caught is
worth recording because two rounds of plausible reasoning missed it.

A peer session challenged the framing on physical grounds: a quenched HTS sits
around 1e-6 to 1e-2 Ohm*m, nowhere near the 1e10 cap, so the caveat could not
mean "the quench regime". Correct. Their replacement was that an E-J power law
driven deep over-critical has `rho ~ (J/Jc)^(n-1)`, which at high `n` exceeds
1e10 for `J` a modest multiple of `Jc` — a transient Newton excursion rather
than a physical state. That is also false, and for a reason neither of us had
looked up: **the material layer caps `rho` structurally.**

All three HTS law families put the power-law channel *in parallel* with the
normal-state channel:

- `rho_powerlaw()`, all eight overloads — `1.0/((1.0/rhon) + (1.0/rhoPL))`
  (`powerlaws.hpp:193, 214, 242, 269, 290, 309, 331, 348`)
- `rho_riva()` — the same, written branch-stably *because* the author
  anticipated this case: `// parallel combination, branch-stable against a huge`
  `rhoPL`
- `rho_piecewise()` — returns the unbounded `rhoPL` branch only below its knee
  `j1 = jc * 10^( 2.5 / n )`, where `rhoPL <= ~10^2.5 * ec / jc`; above the knee
  it returns `rhoFF`, then `rhon`

So `rho <= rhon` always, however far over-critical an iterate drives `|j|`. The
`(J/Jc)^(n-1)` scaling is right about `rhoPL` and irrelevant to the returned
value — `rhoPL` is exactly what the parallel combination discards when it grows.
The upper clamp can only fire if a material's own `rho( T )` exceeds 1e10 Ohm*m
(no shipped material does; a user plugin could) or if an executable narrows the
window (none does). Which is precisely what the source comment at
`cl_Communicator.cpp:62-64` says: *the defaults make the clamp a no-op.*

Both of us had read that comment. Codex quoted its second clause to me
("executables may narrow the window") and I applied that without noticing the
first clause was the load-bearing one. The peer noticed the first clause and
substituted a new mechanism instead of checking what the law returns. Nobody
opened `powerlaws.hpp` until the third round.

§9.3 now leads with why the clamp does not bite, names the parallel combination
as the reason, states the two ways it could still fire, and demotes
`dotQ`-as-lower-bound to real-but-latent: *with stock materials and stock
bounds, treat `dotQ` as unclamped.*

The peer then supplied the piece of evidence that beats both of our arguments,
independently verified here: the author states the design intent in the source.
`powerlaws.hpp:2466-2469` documents that "the power-law channel is evaluated in
log10 space; past a cap the parallel combination is ρn to machine precision, so
the residual returns ρn and the tangents return the matching normal-branch
values instead of the raw inf/inf", and the implementation at `:2503` carries
the same statement at the cap itself — `if ( ! ( lg <= 250.0 ) ) return false ;`,
NaN-aware by construction, after which `rho_riva()` returns `rhon` outright.
`rho_powerlaw()` arrives at the same place by arithmetic rather than by an
early-out: an overflowing `rhoPL` sends `1/(1/rhon + 1/rhoPL)` to `rhon`.
So the no-op is deliberate. §9.3 now cites it, and a reader who has that line
needs neither of our arguments. Cited by greppable sentence rather than by line
number, per the anchor rule.

The durable lesson is not about resistivity. Two AI sessions produced two
confident, physically-reasoned, mutually-inconsistent claims about a clamp whose
behaviour is fixed by four lines of arithmetic in a file neither had opened. The
correct move at round one was to read the law, not to reason about the regime.

## Posture

**Reviewed, not verified.** Nothing was built or run. The rename is mechanical
and the residual grep is clean, but a compile is owed before this is trusted:
the three renamed methods are declared in `cl_FEM_Controller.hpp` and defined in
`cl_FEM_Controller.cpp`, and `save_dotQ` is a header inline with five call
sites in one translation unit.

Still open on this feature, unchanged by this session and inherited from
`7323d528`'s own jury round: no time integration, no per-block decomposition,
no clamp-contamination flag on the reported value, and a two-rank run still
owed. A peer session is rewriting `todo/ac_loss_postprocessing.md` around
exactly that remainder.
