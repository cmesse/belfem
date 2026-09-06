# DR-131's Darwin gate ran: the templates were right, the example decks were broken

**Date:** 2026-09-01
**Purpose:** Record the execution of DR-131's never-run Darwin gate, the four-deck defect it
found instead, and the unrelated `matvec_csc` crash the end-to-end run surfaced (filed as DR-155).
**Module:** `examples/*/src`, `src/physics/materials` (templates), `src/sparse` (by-catch)

---

## What was owed

DR-131 had stood since 2026-08-28 as the round's one residue standing on reasoning alone: both
plugin templates had been given `-Wl,-undefined,dynamic_lookup` on `APPLE` and switched
`SHARED` → `MODULE`, on a source trace by Codex and Grok, and **neither auditor ran a Mac and
neither did Claude.** The row named the `MODULE` change as the higher-risk half — untested on
Darwin, and it changes the artifact name there, so a Mac user's existing
`create_material("./libmyalloy.dylib", …)` would break.

`devlog/dl20260831_plugin_build_against_installed_prefix.md` had closed the installed-prefix half
the day before and closed with the same admission: *"Darwin is untouched and unmeasured."*

This session was run on Christian's Mac — macOS 15.7.9, x86_64, gcc 15 and Apple clang — so the
gate was finally available.

## Build and install

Full rebuild of the shared `build/` tree (exit 0, zero errors; CMake had re-run since 2026-08-27,
so this was effectively cold), then `cmake --install build --prefix /tmp/belfem`, which avoids
touching the tree's cached `CMAKE_INSTALL_PREFIX`. The prefix carries
`bin/{belfem,db2exo,electricalCircuit,gas,material}`, `lib/libbelfem.0.9.0.dylib` + symlinks,
629 headers, and `share/belfem/{examples,templates,material,fluid,python}`.
`scripts/check_doc_claims.py`: 37/37.

One stale entry noticed on the way: `config/globals.cmake:19` still lists `hphirun hphiTrun` in
`BELFEM_INSTALL_EXECUTABLES`, but `src/executables/CMakeLists.txt:39-45` has both commented out as
retired. Install skips them silently; the `build/bin` copies are 2026-08-27 leftovers that no
longer rebuild. Not fixed here — noted so the list is not read as truth.

## The row's own suspicion: refuted

| Gate | Result |
|---|---|
| `UserMaterialTemplate.cmake` + `example_user_material.cpp`, source tree, gcc | `libmyalloy.so`, Mach-O bundle, rc=0 |
| same, Apple `clang++` | rc=0 |
| `UserLibraryTemplate.cmake` + `example_user_defect.cpp` + `example_user_source.cpp` | `libcustom.so`, rc=0 |
| `UserMaterialTemplate.cmake` against `/tmp/belfem` | `-- Using installed BELFEM: /tmp/belfem`, rc=0 |

`nm -gU` shows `_MyAlloy_init` and `_MySuperconductor_init` exported; `nm -gu` shows 11 undefined
`belfem::` symbols left for the host; `otool -L` shows nothing linked but libstdc++ and libSystem.
The artifact is `lib<name>.so`, which is the name the docs and `example_user_material.cpp:15`
actually use.

So `MODULE` plus `dynamic_lookup` is correct on Darwin, and the row's stated risk does not exist:
**no Mac user's `.dylib` load path breaks, because on this tree no Mac user ever had a working one
to break.** Which the next section is about.

## What the gate found instead

The 2026-08-28 fix went into the two templates. It never went into the four example decks, which
are what a user actually runs. Building `examples/tape_quench_usermat/src` unchanged on this Mac:

```
ld: symbol(s) not found for architecture x86_64
make[2]: *** [CMakeFiles/usermat.dir/build.make:100: usermat.dylib] Error 1
```

Two independent fatals, both present in all four decks:

1. **No `dynamic_lookup`.** The decks link `SHARED` and link nothing, so Mach-O's two-level
   namespace rejects the undefined `belfem::` symbols at **link** time. Hard build failure — the
   exact mechanism DR-131 described for the templates.
2. **`if(APPLE) SUFFIX ".dylib"`.** Even had it linked, the artifact is `usermat.dylib` while the
   deck's own `examples/tape_quench_usermat/input.conf:71` names `src/build/usermat.so`. The deck
   could not have loaded its own plugin.

Every shipped deck was dead on arrival on macOS. That is a stronger and more user-visible finding
than the row it came from, so DR-131 was **widened** rather than branched, per the 2026-08-29
extend-don't-branch policy — same contract, same platform, same round.

**Fixed** in all four (`disk_pulse`, `block3d`, `tape_quench_usermat`, `undulator2d`):
`SHARED` → `MODULE`, the `.dylib`/`.dll` suffix overrides dropped so the name matches the deck's
`file :` on every platform, `target_link_options(… "-Wl,-undefined,dynamic_lookup")` under `APPLE`,
and `cmake_minimum_required` 3.10 → 3.13 because `target_link_options` needs it. 4/4 rebuild green
and produce `.so`.

## The dlopen half, which had never been run anywhere

`dl20260831` closed with "nothing was `dlopen`ed by an installed solver binary". Three levels of
evidence, weakest first:

1. A minimal host linking `/tmp/belfem/lib/libbelfem.dylib` loaded the template plugin and
   evaluated it — `rho(77K)=2.301800e-09`, `cp`, `lambda` all live.
2. The same host loaded all five `tape_quench_usermat` materials from the deck's own `usermat.so`,
   including `hts`, which crosses the `JcFunctionUserDefined` vtable — the very symbols that had
   failed to link before the fix.
3. **The real gate:** `gmsh -3 tape.geo`, then the installed `belfem` on the deck. It reached
   295,315 free dofs and entered the first BDF1 step, which is well past material creation.

(3) alone is positive evidence but circumstantial, so it was falsified: hiding `usermat.so` makes
the same run abort at `Could not load library src/build/usermat.so`, naming dyld's full search.
The load path is on the run path, and it passes.

**DR-131's gate is executed.** Retagged `[RUN-BLOCKED]` → `[RUN]`. The strike is Christian's ruling
and was not taken here.

## By-catch: the run then died in the sparse kernel

That end-to-end run did not finish. It crashed in the first timestep with `Bus error: 10`, and the
crash report's faulting thread is `matvec_csc._omp_fn.0` on `gomp_thread_start` — nothing to do
with plugins.

`splinalg.f90:103`'s `!$omp reduction(+:y)` places each worker thread's private copy of the whole
result vector on that thread's **stack**. Isolated in a standalone driver so `n` was the only
variable:

| n | private `y` per thread | 16 threads |
|---|---|---|
| 262,000 | 2046.9 KiB | ok |
| **262,144** | **2048.0 KiB = 2 MiB** | **crash** |
| 295,315 (the deck) | 2307.1 KiB | crash |

Controls at the deck's size: `OMP_STACKSIZE=64M` survives, `OMP_NUM_THREADS=1` survives. The
boundary landing on 2 MiB / 8 B **to the byte** is what identifies libgomp's Darwin worker stack as
the mechanism — measurement, not inference.

Not Darwin-only: Linux's 8 MiB default moves the same cliff to 1,048,576 rows rather than removing
it. That is unmeasured (no Linux box this session) and is logged as O2, because it decides whether
this is a porting note or a live limit on large 3D runs everywhere.

Filed as **DR-155**, plan in `todo/matvec_csc_omp_stack_overflow.md`. **No `src/sparse` change was
made** — O1 (Fortran `allocatable` vs. a `cl_SpMatrix` member scratch) is a hot-path design fork
and is Christian's call. `OMP_STACKSIZE=64M` is a working stopgap today.

## Evidence ladder

**Verified by execution:** the template builds (4 configurations), the four deck builds before and
after the fix, the three dlopen levels plus the falsification control, the install, and the
`matvec_csc` boundary sweep with both controls.

**Reviewed, not verified:** nothing claimed here rests on review alone.

**Not established:** no deck has been run to *completion* on Darwin — DR-155 stops it — so the
plugin machinery is proven to load and evaluate, not to carry a full simulation. The `libc++` /
`libstdc++` ABI question below is reasoned, not measured.

## Residue

- No deck runs to completion on Darwin until DR-155 is fixed.
- The clang template build produces a plugin against `libc++` while the host here is
  gcc/`libstdc++`. Each links cleanly alone; mixed across the `dlopen` boundary that is a silent
  `std::string` ABI mismatch. Neither template nor docs warn about it, and the templates' own
  "keep these in agreement" note covers `NDEBUG` and `BELFEM_INT64` but not the standard library.
  Not filed — raised for a ruling on whether it belongs in the template comment block.
- `config/globals.cmake:19` names two retired executables.
