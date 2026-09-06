# DR-131's last residual discharged: a deck ran to completion on Darwin

**Date:** 2026-09-01
**Topic:** Close out DR-131's "no deck has been run to completion on Darwin" residual, unblocked by
the DR-155 fix; DR-155's own R5 gate run in the same session
**Module:** `examples/`, `src/sparse` (gate execution only — no source change in this session)
**Machine:** macOS 15.7.9 x86_64, gcc 16.1.0 / Apple clang, 8 physical / 16 logical cores

---

## What was owed

DR-131's gate had run on 2026-08-28/09-01 and refuted the row's own suspicion, but it ended with one
honest residual: *"no deck has been run to completion on Darwin, because that same run then died in
`matvec_csc` — a defect unrelated to plugins, filed as DR-155."* With DR-155 fixed, the block was
gone and the residual became executable.

## Setup

A sandbox tree, so nothing touched the shared build trees or the repository working copy:

- `/tmp/build` — CMake tree, **Release** (`USE_DEBUG=OFF`), **`USE_SHARED_LIBS=ON`**,
  `USE_BELFEM_OPENMP=OFF`, `USE_MKL=OFF`, Blaze. Shared is the configuration that actually
  exercises the Darwin `dlopen` / `dynamic_lookup` path this row is about.
- `/tmp/belfem-dr131` — install prefix, left separate from the earlier session's `/tmp/belfem`.
- `/tmp/build/deck`, `/tmp/build/disk` — copies of the two decks, so `gmsh` and Exodus output
  landed outside the repository.

Both decks' plugins were built against the installed prefix with the CMakeLists DR-131 fixed.
Every run used a **default environment with `OMP_STACKSIZE` unset**, verified absent before launch.

## `tape_quench_usermat` — cleared the crash, then stopped deliberately

The step that used to kill it now converges: `BDF1 1 | t : 0.5000 ms`,
`timestep succeeded | magnetic : 1.033e-10 | thermal : 6.684e-12`, 58 s. That is DR-155's R5 bar.

It then ran well past the bar — **54 accepted timesteps, 6 rejected, t = 10.98 ms, 6 h 50 m wall,
zero crash signatures across 2003 log lines** — with five plugin-supplied materials across two
`dlopen`ed libraries, including the `JcFunctionUserDefined` vtable case.

**It was stopped, not failed.** The deck is `simulation time : 0.3 s` of a stiff coupled quench at
295,315 dofs; per-step cost ran 5-13 minutes through the stiff patch, and completion projects to
**~114 hours**. An earlier estimate of ~375 h was too pessimistic — it extrapolated while the
adaptive controller was collapsing Δt to 0.166 ms; Δt later recovered to 0.982 ms. Neither number is
viable, so deck completion was moved to a cheaper in-scope deck.

## `disk_pulse` — ran to completion

`examples/disk_pulse` is one of the four decks DR-131 fixed, and much cheaper: 7,565 nodes /
47,974 elements against `tape_quench`'s 21,303 / 146,265, magnetic-only with cohomology cuts rather
than a coupled quench.

**Result: 106 accepted timesteps, 0 rejected, t = 500.0000 ms — the deck's full `simulation time` —
0 crash signatures in 3085 log lines, 100 Exodus frames written, and the process exited on its
own.** Sub-second timesteps throughout, Δt pinned at its 5 ms ceiling.

**The plugin was proven on the run path, not merely present on disk:** `lsof` on the live pid showed
`/private/tmp/build/disk/src/build/bgpulse.so` mapped, alongside the sandbox
`/private/tmp/belfem-dr131/lib/libbelfem.0.9.0.dylib` it resolves against.

That completion is what 54 `tape_quench` steps could not give: the **teardown path** — output
finalisation and `dlclose` at exit. The plugin *contract* was already better evidenced by
`tape_quench`'s richer surface; what was missing was the exit, and this supplies it.

## By-catch

`hatch_turtle()` fired live at startup on the blind default — *"threads requested by this rank : 16,
PHYSICAL cores available here : 8"*, with its own banner line *"the runtime defaulted to the affinity
mask, which counts every hyperthread. It is not a recommendation."* This is executed confirmation of
a design decision taken earlier the same day: the oversubscription warning was deliberately left on
the plain `OMP` define rather than moved to `BELFEM_OMP`, and under the blanket sweep first proposed
it would have vanished on every default build.

`lsof` also showed `libgomp.1.dylib` mapped with `BELFEM_OMP` off — consistent with the corrected
documentation: the process is still multi-threaded through the third-party solvers and the matrix
backend, just not through BELFEM's own pragmas.

## Status

**DR-131's residual is discharged; the strike is Christian's ruling and is not taken here.**
DR-155's R5 gate is green. Both rows' remaining residue is named in the register:
DR-155 still owes a Linux `-DUSE_MKL=ON` build, the only configuration that can compile the
forwarded MKL branch in `SpMatrix::multiply` at all.

No `src/` change was made in this session — this was gate execution against the code that landed
earlier the same day.
