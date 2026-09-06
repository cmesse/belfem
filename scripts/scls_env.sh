# scls_env.sh — idempotent SCLS toolchain pin for BELFEM tooling. SOURCE this file.
# SCLS_PY_PREFIX names the pinned Python install; override it in the environment
# if yours lives elsewhere. It is only prepended when the directory exists.
# Layout contract: SCLS root /opt/scls with debug/, gcc/, mkl/ subtrees. Only
# /opt/scls/gcc/bin belongs on PATH — debug/ and mkl/ are library stacks that
# CMake discovers; Python is NOT part of SCLS and has its own pinned prefix.

# SCLS is a Linux-only stack: the pinned prefixes below do not exist on macOS
# (which carries its own toolchain, see the g++ note further down). Skip the
# pin entirely on any other OS instead of failing every non-interactive script.
case "$(uname -s)" in
    Linux) : ;;
    *) return 0 2>/dev/null || exit 0 ;;
esac

SCLS_GCC_BIN="/opt/scls/gcc/bin"
SCLS_PY_PREFIX="${SCLS_PY_PREFIX:-$HOME/Applications/python}"

# idempotent prepend: never duplicate, never reorder an already-correct shell
case ":$PATH:" in *":$SCLS_GCC_BIN:"*) : ;; *) PATH="$SCLS_GCC_BIN:$PATH" ;; esac
if [ -d "$SCLS_PY_PREFIX/bin" ]; then
    case ":$PATH:" in *":$SCLS_PY_PREFIX/bin:"*) : ;; *) PATH="$SCLS_PY_PREFIX/bin:$PATH" ;; esac
fi
export PATH

# guards — a wrong toolchain must not run silently from scripts/hooks;
# interactive shells get a warning only.
# g++ is only required to EXIST: on Linux the system compiler is preferred by
# policy — SCLS ships no compiler binary; its mpicxx wraps the bare system g++
# (`mpicxx -show` => "g++ -I/opt/scls/gcc/include -L/opt/scls/gcc/lib ...").
# Exceptions carry their own toolchains and would need their own pin:
# Rocky 8 (system gcc too old) and macOS (clang unless a real gcc is installed).
scls_env_fail=0
for scls_tool in cmake mpicxx ; do
    case "$(command -v "$scls_tool" 2>/dev/null)" in
        /opt/scls/*) : ;;
        *) echo "scls_env.sh: $scls_tool resolves to '$(command -v "$scls_tool" 2>/dev/null || echo MISSING)', not under /opt/scls/" >&2
           scls_env_fail=1 ;;
    esac
done
unset scls_tool
if ! command -v g++ >/dev/null 2>&1 ; then
    echo "scls_env.sh: g++ not found on PATH (SCLS mpicxx wraps the system g++)" >&2
    scls_env_fail=1
fi
case "$(command -v python3 2>/dev/null)" in
    "$SCLS_PY_PREFIX"/*) : ;;
    *) echo "scls_env.sh: python3 resolves to '$(command -v python3 2>/dev/null || echo MISSING)', not under pinned prefix $SCLS_PY_PREFIX" >&2
       scls_env_fail=1 ;;
esac

if [ "$scls_env_fail" -ne 0 ]; then
    case "$-" in
        *i*) echo "scls_env.sh: WARNING — wrong toolchain, continuing because this shell is interactive" >&2
             unset scls_env_fail ;;
        *)   unset scls_env_fail
             return 1 2>/dev/null || exit 1 ;;
    esac
else
    unset scls_env_fail
fi
