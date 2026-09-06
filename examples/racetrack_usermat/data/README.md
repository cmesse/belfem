# racetrack_usermat plugin data

**Date:** 2026-08-31
**Purpose:** Measured tables read by the plugins in `../src`.

This directory is on the plugin include path as a compiled-in absolute path:
`src/CMakeLists.txt` passes it to the compiler as `BELFEM_USER_DATA_DIR`, so a
plugin that reads a table here works no matter which directory the solver is
launched from.

The deck currently ships no tables — `src/matlib.cpp` defines the bulk
conductor from constants alone. To add a measured property, drop a two-column
text file here and read it with the `usermat::Table` helper (see
`../../tape_quench_usermat/src/user_table.hpp`, which is copied per example):

```cpp
real bulk_rho( const Material* mat, real T )
{
    static usermat::Table tTable( "rho_bulk.txt" ) ;
    return tTable( T ) ;
}
```

No build change is needed — the data directory is already wired.
