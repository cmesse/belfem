# Summarise the Doxygen problem log after 'make doc'.
#
# Doxygen writes to a file (WARN_LOGFILE) so the output does not scroll past in the
# build. It logs two severities and they are not interchangeable: `warning:` for
# documentation defects (a bad @param, an unresolved @ref) and `error:` for output
# it failed to produce at all. A run can emit zero warnings while silently dropping
# every diagram in the site, so both are counted here.
#
# History worth keeping: until 2026-08-18 this script matched `warning:` only. The
# tree was carrying 1471 `error:` lines -- one per graph that Doxygen requested and
# Graphviz never delivered -- and none of them ever reached the build output.
#
# Note that Doxygen also writes CONFIGURATION complaints (obsolete tags, options not
# compiled in) to stderr rather than to this log, so a clean log is not by itself proof
# of a clean run -- watch the console too.
#
# THE WARNING BASELINE IS VERSION-BOUND. It was 0/0 on 2026-08-18 against Doxygen
# 1.18.0. Measured 2026-08-28 on Doxygen 1.9.1 it is 90/0, and every one of the 90 is
# in a single file, doc/lessons_learned_evidence.md. Nothing else in the tree warns.
#
# Those 90 are a 1.9.1 limitation rather than bad markup: 80 are "found </tt> tag
# without matching <tt>", and backtick parity is even on every line of that file, a
# synthetic table reproducing the same constructs emits nothing, and the warnings only
# begin ~40 rows into a single 110-row, 11-column table. The remaining 10 are \ref and
# \includedoc complaints from the same tables. Clearing them would mean rewriting a
# file that is deliberately dense and is exempt from prose sweeps by CLAUDE.md, so the
# file is quarantined in this number instead.
#
# What is NOT hidden here: as of 2026-08-28 every warning outside that one file has
# been fixed, not baselined. If this number rises, something real broke.
#
# Expect a DIFFERENT number on a different Doxygen. On an upgrade, re-measure and reset
# rather than assuming a regression.
#
# Invoked as: cmake -DLOG=<path> -P report_doxygen_warnings.cmake

set(BELFEM_DOXYGEN_WARNING_BASELINE 90)

# Errors mean output is MISSING, never merely undocumented, so this stays at zero.
# It was 1 until 2026-08-28 (a version-bound layout file the installed Doxygen
# rejected); the layout is now generated for the installed version, so 0 is real.
set(BELFEM_DOXYGEN_ERROR_BASELINE 0)

if(NOT EXISTS "${LOG}")
    message(STATUS "Doxygen: no problem log at ${LOG}")
    return()
endif()

file(STRINGS "${LOG}" _warn_lines REGEX "warning:")
file(STRINGS "${LOG}" _error_lines REGEX "error:")
list(LENGTH _warn_lines _warnings)
list(LENGTH _error_lines _errors)

# Errors mean output is missing, so they are reported first and never tolerated
# above the baseline.
if(_errors GREATER BELFEM_DOXYGEN_ERROR_BASELINE)
    message(WARNING
        "Doxygen reported ${_errors} errors -- output is MISSING from the site, not "
        "merely undocumented. The usual cause is a graph backend that could not run; "
        "check that Graphviz 'dot' is on PATH, then read ${LOG}.")
endif()

if(_warnings GREATER BELFEM_DOXYGEN_WARNING_BASELINE)
    math(EXPR _new "${_warnings} - ${BELFEM_DOXYGEN_WARNING_BASELINE}")
    message(WARNING
        "Doxygen produced ${_warnings} warnings, ${_new} more than the known baseline "
        "of ${BELFEM_DOXYGEN_WARNING_BASELINE}. Read ${LOG} and fix the new ones.")
elseif(_warnings LESS BELFEM_DOXYGEN_WARNING_BASELINE)
    message(STATUS
        "Doxygen: ${_warnings} warnings, below the baseline of "
        "${BELFEM_DOXYGEN_WARNING_BASELINE}. Lower BELFEM_DOXYGEN_WARNING_BASELINE in "
        "config/doc/report_doxygen_warnings.cmake to lock the improvement in.")
else()
    message(STATUS
        "Doxygen: ${_warnings} warnings (known baseline), ${_errors} errors. Log: ${LOG}")
endif()
