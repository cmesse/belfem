# Make the project logo a link to the BELFEM project site.
#
# Doxygen has no configuration option for this: PROJECT_LOGO is emitted as a bare
# <img> inside <td id="projectlogo">, and the only way to wrap it is a custom
# HTML_HEADER.
#
# A custom header is NOT checked into the repository, deliberately. Doxygen does not
# warn when a custom header is stale -- verified 2026-08-18 on 1.18.0, where a header
# missing $treeview and $mathjax produced a silent build with both features simply
# absent. A header file checked in against one Doxygen version would therefore lose
# the navigation tree, search, dark mode or MathJax on a future upgrade, with no
# warning and no error. Instead the `doc` target regenerates the header from the
# installed Doxygen with `doxygen -w html` and this script applies the one edit.
#
# The markup is matched by REGEX, not as a literal, because it is version-dependent:
# Doxygen 1.18 emits an extra $logosize marker inside the <img> tag that 1.9.1 does
# not. A literal needle pinned to one version turns every other version into a build
# failure -- which is exactly what happened on 1.9.1, where the 1.18-shaped needle
# aborted `make doc`. The pattern below accepts any <img> attribute set; only a
# change to the surrounding <td id="projectlogo"> structure is treated as a genuine
# template change worth failing on.
#
# Invoked as: cmake -DHEADER=<file> -DURL=<url> -P patch_doxygen_header.cmake

if(NOT EXISTS "${HEADER}")
    message(FATAL_ERROR "patch_doxygen_header: no header at ${HEADER}")
endif()

file(READ "${HEADER}" _html)

# One pattern matches BOTH the unpatched template and an already-wrapped one, so the
# rewrite is idempotent by construction and a stale or wrong URL is CORRECTED rather
# than accepted. An earlier version returned early on any "<a href=" it found, which
# silently kept whatever link was already there.
#
# The <img> attributes vary by Doxygen version -- 1.18 adds a $logosize marker that
# 1.9.1 does not emit -- so they are matched loosely. The <td> wrapper is what has to
# hold; if that changes, Doxygen has restructured the header and we want to know.
string(REGEX MATCH
       "<td id=\"projectlogo\">(<a[^>]*>)?(<img[^>]*/>)(</a>)?</td>"
       _needle "${_html}")

if(NOT _needle)
    # Fail loudly. Silently shipping an unlinked logo is how this kind of change rots.
    message(FATAL_ERROR
        "patch_doxygen_header: the project-logo markup in the generated header no longer "
        "matches what this script expects. Doxygen has changed its header template.\n"
        "Re-derive the expected markup with:\n"
        "    cd $(mktemp -d) && doxygen -w html header.html footer.html style.css\n"
        "(run it in an EMPTY directory -- in a build tree Doxygen picks up ./Doxyfile and\n"
        "aborts on HTML_HEADER before writing anything), then update the regex in\n"
        "config/doc/patch_doxygen_header.cmake.")
endif()

set(_img "${CMAKE_MATCH_2}")
set(_linked "<td id=\"projectlogo\"><a href=\"${URL}\" title=\"BELFEM project site\">${_img}</a></td>")
string(REPLACE "${_needle}" "${_linked}" _html "${_html}")
file(WRITE "${HEADER}" "${_html}")
