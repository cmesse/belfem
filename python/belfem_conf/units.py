"""The unit table, derived from `unit_to_si` rather than transcribed.

`stringtools.cpp` carries ~800 lines of

    else if ( tUnit == "V" )
    {
        tScale = 1.0;
        tMass += tPower;
        tLength += tPower * 2;
        tCurrent -= tPower;
        tTime -= tPower * 3;
    }

which is a machine-readable table: a token, a scale, and the exponent it adds to
each of the seven SI dimensions. Copying that into Python would produce yet
another artifact to fall out of step — the same mistake the enum lists nearly
made. So it is parsed from the source at run time.

Only the DIMENSION matters for validation. BELFEM's own unit checks compare the
dimension and not the token, which is why `kA` is accepted where `A` is expected
and `mum` where `m` is: same dimension, different scale.
"""

from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path

# `aValue.second` in the C++, in order.
DIMENSIONS = ("length", "mass", "time", "current",
              "temperature", "substance", "brightness")

_VAR_TO_DIM = {
    "tLength": "length",
    "tMass": "mass",
    "tTime": "time",
    "tCurrent": "current",
    "tTemperature": "temperature",
    "tSubstance": "substance",
    "tBrightness": "brightness",
}

# A branch may test several spellings at once —
#   else if ( tUnit == "K" || tUnit == "°K" || tUnit == "C" || tUnit == "°C" )
# so the whole condition is captured and every quoted token taken from it.
# Matching only a single `== "X" )` silently dropped every multi-token branch,
# temperature and resistance among them.
_BRANCH = re.compile(
    r'if\s*\(\s*((?:tUnit\s*==\s*"[^"]+"\s*(?:\|\|\s*)?)+)\)\s*\{(.*?)\n(\s*)\}',
    re.S,
)
_TOKENS = re.compile(r'tUnit\s*==\s*"([^"]+)"')
_BUMP = re.compile(
    r"(tLength|tMass|tTime|tCurrent|tTemperature|tSubstance|tBrightness)"
    r"\s*([+-])=\s*tPower(?:\s*\*\s*([0-9.]+))?"
)

# Dimension signatures for the names the schema uses, so a schema saying
# `dimension: length` can be compared with a parsed token.
NAMED: dict[str, dict[str, float]] = {
    "dimensionless": {},
    "length": {"length": 1},
    "mass": {"mass": 1},
    "time": {"time": 1},
    "current": {"current": 1},
    "temperature": {"temperature": 1},
    "frequency": {"time": -1},
    "angle": {},                       # radians are dimensionless in this scheme
    "voltage": {"mass": 1, "length": 2, "current": -1, "time": -3},
    "current_density": {"current": 1, "length": -2},
    "electric_field": {"mass": 1, "length": 1, "current": -1, "time": -3},
    "magnetic_field": {"current": 1, "length": -1},
    "heat_flux": {"mass": 1, "time": -3},
}

# Schema spellings that name a unit rather than a dimension.
UNIT_ALIASES = {
    "a": "current", "v": "voltage", "k": "temperature", "s": "time",
    "m": "length", "hz": "frequency", "a/m": "magnetic_field",
    "a/m^2": "current_density", "v/m": "electric_field",
    "w/m^2": "heat_flux", "ohm": "resistance",
}


@lru_cache(maxsize=4)
def table(root: str) -> dict[str, dict[str, float]]:
    """{unit token: {dimension: exponent}} straight out of unit_to_si."""
    src = (Path(root) / "src/core/stringtools.cpp").read_text(errors="replace")
    try:
        body = src[src.index("unit_to_si( const string"):]
    except ValueError:
        return {}

    out: dict[str, dict[str, float]] = {}
    for match in _BRANCH.finditer(body):
        condition, block = match.group(1), match.group(2)
        if "tScale" not in block:
            continue                       # not a unit branch
        dims: dict[str, float] = {}
        for var, sign, factor in _BUMP.findall(block):
            step = float(factor) if factor else 1.0
            dim = _VAR_TO_DIM[var]
            dims[dim] = dims.get(dim, 0.0) + (step if sign == "+" else -step)
        signature = {d: e for d, e in dims.items() if e != 0.0}
        for token in _TOKENS.findall(condition):
            out[token] = signature
    return out


_ATOM = re.compile(r"^([^\^]+?)(?:\^([+-]?[0-9.]+))?$")


# Temperature spellings that `Section::create_key` converts BEFORE unit_to_si
# ever sees them (°F and Rankine never reach the table), so they must be known
# here or a legal `x : 70 °F ;` reports as an unknown unit.
_CREATE_KEY_TEMPERATURES = {"°F", "R", "°R"}


def _preprocess(token: str) -> str:
    """The rewrites `unit_to_si` applies before matching tokens.

    Mirrored from stringtools.cpp: spaces and brackets removed, `*` handled by
    the caller's split, `µ` -> `mu`, `²` -> `^2`, `³` -> `^3`. Skipping these
    rejected `µm` and `mm²`, both legal to the C++.
    """
    token = token.replace(" ", "").replace("(", "").replace(")", "")
    token = token.replace("µ", "mu")
    return token.replace("²", "^2").replace("³", "^3")


def dimension_of(token: str, root: str) -> dict[str, float] | None:
    """Dimension signature of a unit token, or None if unrecognised.

    Compound units are composed, not looked up: `unit_to_si` splits on `/` and
    `*` and reads a trailing `^N` as a power, so `A/m^2` never appears in the
    table as a whole. Treating the table as a flat token list rejected `A/m^2`
    and `V/m` — both of which sit in four of the six example decks.
    """
    token = _preprocess(token.strip())
    if token in ("", "-"):
        return {}
    if token in _CREATE_KEY_TEMPERATURES:
        return {"temperature": 1}

    known = table(root)
    if token in known:                    # atomic, and possibly irregular
        return known[token]

    # unit_to_si splits at the FIRST '/' only; a second '/' stays glued to a
    # denominator word, matches no branch, and aborts. `W/m/K` is illegal to
    # BELFEM, so accepting it here would be a false pass.
    if token.count("/") > 1:
        return None

    total: dict[str, float] = {}
    sign = 1.0
    resolved = 0
    # split keeping the separators, so everything after the '/' inverts
    for piece in re.split(r"([*/])", token):
        if piece == "/":
            sign = -1.0
            continue
        if piece in ("*", ""):
            continue
        m = _ATOM.match(piece.strip())
        if not m:
            return None
        base, power = m.group(1), m.group(2)
        dims = known.get(base)
        if dims is None:
            return None
        resolved += 1
        exponent = float(power) if power else 1.0
        for dim, value in dims.items():
            total[dim] = total.get(dim, 0.0) + sign * exponent * value

    if not resolved:
        return None                       # nothing but separators
    return {d: e for d, e in total.items() if e != 0.0}


def expected(name: str) -> dict[str, float] | None:
    """Resolve a schema `dimension:` spelling to a signature."""
    key = str(name).strip().lower()
    key = UNIT_ALIASES.get(key, key)
    if key == "resistance":
        return {"mass": 1, "length": 2, "current": -2, "time": -3}
    return NAMED.get(key)


def describe(sig: dict[str, float]) -> str:
    if not sig:
        return "dimensionless"
    return " ".join(
        f"{d}^{e:g}" if e != 1 else d for d, e in sorted(sig.items())
    )
