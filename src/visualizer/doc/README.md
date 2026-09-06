# Visualizer Module Documentation {#visualizer_index}

**Module:** `src/visualizer`
**Purpose:** Index of documentation for BELFEM's optional VTK-based visualization

---

## Overview

The `visualizer` module uses VTK to render meshes and curves and provides the
`visualize` executable.

**This module is opt-in and OFF by default.** Enable it with `-DUSE_VTK=ON`
(`CMakeLists.txt`: `option( USE_VTK "Use Visualization Toolkit" OFF )`).
`src/CMakeLists.txt` adds the directory only under that guard. A default build contains
neither the library nor the executable, and nothing else in the tree depends on it. The
usual path is to let BELFEM write Exodus, VTK and HDF5 files and open them directly in
external viewers such as ParaView or VisIt.

---

## Key Classes

| Class | Files | Role |
|-------|-------|------|
| **`vtk::MeshView`** | cl_VTK_MeshView.{hpp,cpp} | Renders a mesh |
| **`vtk::Curve`** | cl_VTK_Curve.{hpp,cpp} | Renders a curve |
| `vtk::BlockActor` | cl_VTK_MeshView.{hpp,cpp} | VTK actor for one mesh block |
| `vtk::SideSetActor` | cl_VTK_MeshView.{hpp,cpp} | VTK actor for one sideset |

The classes live in namespace `belfem::vtk`; the `VTK_` prefix is on the file names, not
the class names.

`vtktypes.hpp` provides the `belfem::vtk` smart-pointer aliases (`Actor`, `Points`, `Renderer`, ...);
the element-type map `vtk_type()` lives in the mesh module (`src/mesh/vtktools.hpp`).
`visualize.cpp` provides `main()` for the `visualize` executable.

---

## See Also

- **Source code** - header files in `src/visualizer/`
- [Mesh module](../../mesh/doc/README.md) - mesh I/O, including the VTK writer that the
  external-viewer path uses
- [Getting Started](../../../doc/getting_started.md) - build options

---

**Status:** this index was written from the module's headers and its CMake guard. Because
the module is off in a default build it is the least exercised part of the tree; treat
anything here as unverified until you have built with `-DUSE_VTK=ON`.
