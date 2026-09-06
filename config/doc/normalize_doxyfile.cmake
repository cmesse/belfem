# Produce the working Doxygen configuration: copy the pristine configure_file output,
# then normalise the copy to the INSTALLED Doxygen's tag set with `doxygen -u`.
#
# The copy is not optional. `make doc` does not re-run CMake, and the Doxygen version
# is not a configure dependency, so normalising the configure_file output in place
# would leave a permanently downgraded file in the build tree -- and after a later
# Doxygen upgrade `-u` would re-add the newer tags at their DEFAULTS instead of the
# values configured from Doxyfile.in. Copy-then-normalise keeps the result a pure
# function of (pristine configuration, installed Doxygen).
#
# `-u` must parse the pristine file to rewrite it, so it reports every tag the
# installed Doxygen does not know -- 46 of them when a 1.18-maintained Doxyfile.in
# meets Doxygen 1.9.1. That is expected migration chatter, not a defect, and printing
# it on every `make doc` buries anything that matters. It is summarised to one line
# here. Any OTHER diagnostic `-u` produces is passed through untouched.
#
# Invoked as: cmake -DDOXYGEN=<exe> -DPRISTINE=<config> -DWORKING=<config>
#                   -P normalize_doxyfile.cmake

if(NOT EXISTS "${PRISTINE}")
    message(FATAL_ERROR "normalize_doxyfile: no configuration at ${PRISTINE}")
endif()

configure_file("${PRISTINE}" "${WORKING}" COPYONLY)

execute_process(
    COMMAND "${DOXYGEN}" -u "${WORKING}"
    RESULT_VARIABLE _rc
    OUTPUT_QUIET
    ERROR_VARIABLE _err)

if(NOT _rc EQUAL 0)
    message(FATAL_ERROR "normalize_doxyfile: 'doxygen -u' failed (${_rc}):\n${_err}")
endif()

# Split the expected tag chatter from anything else, and keep the rest visible.
set(_ignored 0)
set(_other "")
string(REPLACE "\n" ";" _lines "${_err}")
foreach(_line IN LISTS _lines)
    if(_line MATCHES "ignoring unsupported tag")
        math(EXPR _ignored "${_ignored} + 1")
    elseif(_line MATCHES "[^ \t]")
        string(APPEND _other "${_line}\n")
    endif()
endforeach()

if(_ignored GREATER 0)
    message(STATUS
        "Doxygen: dropped ${_ignored} configuration tag(s) this Doxygen does not "
        "support while normalising Doxyfile.in to the installed version. The features "
        "they configure are unavailable until Doxygen is upgraded.")
endif()

if(_other)
    message(STATUS "Doxygen reported while normalising the configuration:\n${_other}")
endif()
