# Write a copy of a Doxygen configuration with HTML_HEADER blanked.
#
# Why this exists: `doxygen -w html <header> <footer> <css> <config>` validates the
# config's HTML_HEADER for EXISTENCE before it writes anything. Pointing that config
# at the very file the command is about to create is a chicken-and-egg -- Doxygen
# exits 1 with "header file '...' does not exist" and the `doc` target dies on its
# second command. Passing no config at all avoids the check, but Doxygen's own
# commentary (Doxyfile.in, HTML_HEADER section) says the generated template "is
# dependent on the configuration options used (e.g. the setting GENERATE_TREEVIEW)"
# and recommends passing the project config. So we pass the project config with the
# one self-referential tag removed, which is both config-aware and unable to trip
# the existence check.
#
# Only HTML_HEADER is blanked. HTML_FOOTER and HTML_EXTRA_STYLESHEET are already
# empty in this project; if either is ever set to a generated file, it needs the
# same treatment here.
#
# Invoked as: cmake -DIN=<config> -DOUT=<config copy> -P strip_doxygen_html_header.cmake

if(NOT EXISTS "${IN}")
    message(FATAL_ERROR "strip_doxygen_html_header: no configuration at ${IN}")
endif()

file(READ "${IN}" _cfg)

# Match only a real assignment, so the commented-out occurrences in Doxygen's own
# explanatory block above the tag are left alone: a '#' after the newline cannot match.
# Leading whitespace is tolerated because nothing guarantees a future `doxygen -u`
# keeps assignments left-aligned.
string(REGEX REPLACE "(^|\n)[ \t]*HTML_HEADER[ \t]*=[^\n]*"
       "\\1HTML_HEADER            =" _cfg "${_cfg}")

file(WRITE "${OUT}" "${_cfg}")

# Verify, rather than assume. If the tag were left pointing at the header this config
# is used to generate, `doxygen -w html` would abort on the existence check and take
# the whole `doc` target with it -- the original defect this script was written for.
file(READ "${OUT}" _check)
if(NOT _check MATCHES "(^|\n)[ \t]*HTML_HEADER[ \t]*=[ \t]*(\n|$)")
    message(FATAL_ERROR
        "strip_doxygen_html_header: HTML_HEADER is not blank in ${OUT}. Doxygen would "
        "validate it for existence and abort before writing the header template.")
endif()
