# Material Data Search Path and the share/ Directory Rename

**Date:** 2026-08-12
**Purpose:** Record the `$BELFEM_DATA/material` fallback for material data files
and the `fluidprop` → `fluid` / `matprop` → `material` rename of the shared data
directories
**Module:** `src/physics/materials`, `src/physics/gastables`

## Motivation

The undulator example in `cmake-build-debug/undulator` names two material data
files in its `input.conf`:

```
ybco { builtin : ybco ; file : sp-ap.hdf5 ; ... }
iron { builtin iron ; bhfile: MatData/bhdata.hdf5 ; curve: RoxieIron ; }
```

Neither file was in the run directory: the databases live in the shared data
tree, and the run directory's `MatData` is empty. Both paths therefore had to be
copied into every run directory, or spelled as machine-specific absolute paths.
The gas tables already solved this — `gastables::data_path()` resolves them
through `gBelfemDataPath`, set from `$BELFEM_DATA` by
`Communicator::set_globals()` — but the materials module had no equivalent and
handed the input string straight to `HDF5` or `dlopen`.

## Directory rename

The shared data directories were renamed, so that the subdirectory names read as
what they hold rather than as an abbreviation:

| before | after |
|---|---|
| `$BELFEM_DATA/fluidprop` | `$BELFEM_DATA/fluid` |
| `$BELFEM_DATA/matprop` | `$BELFEM_DATA/material` |

`share/material` also absorbs the HDF5 databases that used to sit in the
repository's `database/` directory (`bhdata.hdf5`, `sst-1.hdf5`, `sp-ap.hdf5`).

`fluidprop` survived only inside the gastables module and the table-building
scripts, and was converted in place:

- `fn_GT_data_path.cpp` — `gGastablesSubdir` and the four relative fallbacks
- `fn_GT_data_path.hpp`, `cl_GT_RefGasFactory.{hpp,cpp}` — doc and error text
- `src/physics/CMakeLists.txt` — module comment
- `scripts/fluidprop/{build_tables.py,README.md}` — `--out-dir` targets

The script directory keeps its name; it is not part of `$BELFEM_DATA`. The
`.gitignore` that documents the provenance of the four generated `.inp` files
moved with them into `share/fluid`.

## The resolver

New: `src/physics/materials/fn_material_data_path.{hpp,cpp}`.

```cpp
string material::data_path();                    // $BELFEM_DATA/material, or ""
string material::data_file( const string & aFile );
```

`data_file()` tries three locations in order:

1. the path as written, relative to the run directory, or absolute;
2. the same relative path below `data_path()`;
3. the **file name alone** below `data_path()`.

Three design points, each of which was a choice rather than a default:

**The run directory wins.** Step 1 comes first so that a local copy of a
database always overrides the shared one. A search path that consulted
`$BELFEM_DATA` first would make a run's own data silently inert.

**Step 3 exists for `MatData/bhdata.hdf5`.** An input file carries the directory
layout of the machine it was written on. Falling back to the bare file name lets
the same `input.conf` run anywhere the file is in the shared tree, which is what
made the undulator example work without editing it. `examples/sidecoating`
(`MatData/sp-ap.hdf5`) benefits identically.

**An unresolved path is returned unchanged, not turned into an error.** Two
reasons. The opener already reports it — `HDF5` in `OPEN_RDONLY` mode checks
`file_exists( mPath )` and names the file — and reporting it there names the
file the user actually wrote. More importantly, two of the four call sites are
`dlopen`, which searches `$LD_LIBRARY_PATH` for a name that is not a path at
all; a `BELFEM_ERROR` in the resolver would have removed a lookup that works
today.

## Call sites

All four paths that an input file can name are routed through the resolver:

| key | site |
|---|---|
| `bhfile` | `MaterialFactory::create_bh_curve` |
| HTS `file` | `MaterialFactory::create_jc_function( path, label )` |
| `custom { file }` | `MaterialFactory::create_material( libraryPath, label )` |
| `defect { file }` | `Material::read_defect` |

The B-H branch of the `MaterialFactory` constructor built a
`material::BhSplineCurve` directly, bypassing its own factory method; it now
calls `this->create_bh_curve()`, which is where the resolution happens. That was
the only place in the tree that constructed one of these loaders outside the
materials module — checked by grep for `BhSplineCurve`, `JcFunctionDatabase` and
`UserDefinedMaterial(`.

## Evidence

Reviewed, not verified end-to-end: no solver run was made.

- The three touched translation units syntax-check clean under
  `g++ -std=gnu++17 -fsyntax-only` with the `build/` tree's own flags.
- A standalone scratchpad probe linked `fn_material_data_path.cpp` against
  `filetools.cpp` and exercised the resolver from two working directories.
  From `cmake-build-debug/undulator`, with `BELFEM_DATA` pointing at `share`:

  ```
  sp-ap.hdf5                 ->  .../share/material/sp-ap.hdf5     (step 2)
  MatData/bhdata.hdf5        ->  .../share/material/bhdata.hdf5    (step 3)
  share/material/sst-1.hdf5  ->  .../share/material/sst-1.hdf5     (step 3)
  /no/such/file.hdf5         ->  /no/such/file.hdf5                (unresolved)
  ```

  Run from the repository root, `share/material/sst-1.hdf5` resolves to itself
  at step 1 — the local-copy-wins case.
- With `gBelfemDataPath` empty, every input is returned unchanged, so a tree
  without `$BELFEM_DATA` behaves exactly as before this change.
- `scripts/check_doc_claims.py`: 32/32 claims hold.

Not exercised: the two `dlopen` call sites, and any MPI rank count above one.
Resolution is a pure function of the filesystem, which every rank sees alike, so
ranks cannot disagree — but that is an argument, not a run.

## Documentation

Per the input contract rule, both artifacts were updated in this session:

- `doc/input_file_reference.md` §5 — a "Where a material file is looked for"
  paragraph covering all four keys
- `doc/input_schema.yaml` — a `path_resolution` block on the `materials`
  section, anchored on `material::data_file`, with `resolution:
  material_data_path` tagged onto each of the four `type: path` keys

Also updated: `src/physics/materials/doc/README.md`, and the Doxygen parameter
docs on the three factory methods and `Material::read_defect`.
