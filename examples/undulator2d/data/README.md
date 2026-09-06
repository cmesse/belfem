# undulator2d plugin data

**Date:** 2026-08-31
**Purpose:** Measured tables read by the plugins in `../src`.

This directory is on the plugin include path as a compiled-in absolute path:
`src/CMakeLists.txt` passes it to the compiler as `BELFEM_USER_DATA_DIR`, so a
plugin that reads a table here works no matter which directory the solver is
launched from.

The deck currently ships no tables — `src/current.cpp` defines its transport
current analytically as a trapezoid. To drive the coils from a measured
waveform instead, drop a two-column `t I` text file here and read it with the
`usermat::Table` helper (see `../../tape_quench_usermat/src/user_table.hpp`,
which is copied per example):

```cpp
real my_current( const real t )
{
    static usermat::Table tTable( "I_vs_t.txt" ) ;
    return tTable( t ) ;
}
```

No build change is needed — the data directory is already wired.
