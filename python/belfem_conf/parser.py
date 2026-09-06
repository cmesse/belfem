"""A lossless parser for `input.conf`.

Design premise: **the original text is the source of truth.** Every node holds
the byte span it came from, and nothing is normalised on the way in. That buys
two things:

  * round-tripping an unedited deck returns the original bytes, trivially and
    provably — no formatter can drift;
  * an edit later rewrites only its own span, so hand-written comments, tabs and
    layout survive. The example decks are hand-commented and git-tracked, and a
    tool that reflows them turns a one-key change into a whole-file diff.

The C++ parser normalises aggressively (lowercases keys and section types,
collapses whitespace, drops comments) because it only has to *read*. This one
has to read AND write, so it keeps the raw and normalises on access instead.

Structure recognised, following `cl_InputFile.cpp` / `cl_Input_Section.cpp`:

    // comment, or # comment          -> trivia, preserved verbatim
    type [ : label ] { ... }          -> section; header is the text before `{`
    key [ : value ] ;                 -> statement; a bare `key ;` is a flag

Comments are scanned first, because a `{`, `}` or `;` inside one is not
structure. Statement values may themselves contain `:` — `sidesets : 4:9 ;` is
one key and one value — so a statement splits on its FIRST colon only.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class Statement:
    """One `key : value ;` (or bare `key ;` flag)."""

    span: tuple[int, int]          # whole statement including the ';'
    key_span: tuple[int, int]
    value_span: tuple[int, int] | None
    text: str                      # the document's full text, for slicing
    # comment mask over the whole text. A statement span runs from the previous
    # terminator, so it can contain whole comment lines — the mask keeps them
    # out of `key`/`value` exactly as the C++ strips them before parsing.
    # Slicing raw would fold `// note` into the next key name.
    mask: bytearray | None = None

    def _clean(self, start: int, end: int) -> str:
        if self.mask is None:
            return self.text[start:end]
        return "".join(self.text[i] for i in range(start, end)
                       if not self.mask[i])

    @property
    def key(self) -> str:
        return self._clean(*self.key_span).strip()

    @property
    def name(self) -> str:
        """Lowercased key, which is how the C++ stores and looks it up."""
        return " ".join(self.key.lower().split())

    @property
    def value(self) -> str:
        if self.value_span is None:
            return "true"          # a bare flag is stored as the string "true"
        return " ".join(self._clean(*self.value_span).split())

    @property
    def key_offset(self) -> int:
        """Offset of the key's first real character, for diagnostics.

        The span starts at the previous terminator, often a line (or a comment
        block) above — pointing a message there sends the reader to the wrong
        line.
        """
        for i in range(*self.key_span):
            if not self.text[i].isspace() and (self.mask is None or not self.mask[i]):
                return i
        return self.span[0]

    @property
    def raw(self) -> str:
        return self.text[slice(*self.span)]


@dataclass
class Section:
    span: tuple[int, int]          # header through closing brace
    header_span: tuple[int, int]
    body_span: tuple[int, int]     # inside the braces, exclusive
    text: str
    statements: list[Statement] = field(default_factory=list)
    sections: list["Section"] = field(default_factory=list)
    parent: "Section | None" = None
    # The header with comments removed, reduced to the last line before the
    # brace as the C++ does. Stored rather than derived, because the span runs
    # from wherever the previous statement ended: for the first section in a
    # file it swallows the licence banner. Keeping the span intact matters for
    # editing; keeping the TEXT right matters for reading.
    header_text: str = ""
    # Non-comment text in the header span that the C++ would discard. Almost
    # always empty; when it is not, the deck contains something inert.
    header_discarded: str = ""
    # Non-comment text between the last terminator and this section's closing
    # brace (or, on the root, the end of file). A statement missing its `;` is
    # the usual cause, and the C++ drops such a line without a word — see
    # cl_Input_Section.cpp, which keys only on lines containing ';'.
    trailing_discarded: str = ""
    # offset of the discarded trailing text, for diagnostics
    trailing_offset: int = 0
    # offset of the header line itself, for diagnostics
    header_offset: int = 0

    @property
    def header(self) -> str:
        return self.header_text.strip()

    @property
    def header_raw(self) -> str:
        """The span verbatim, comments and all — for editing, not for reading."""
        return self.text[slice(*self.header_span)]

    @property
    def type(self) -> str:
        head = self.header.split(":", 1)[0]
        return " ".join(head.lower().split())

    @property
    def label(self) -> str:
        parts = self.header.split(":", 1)
        if len(parts) < 2:
            return ""
        return " ".join(parts[1].lower().split())

    @property
    def body(self) -> str:
        return self.text[slice(*self.body_span)]

    # -- lookup, mirroring the C++ accessors ------------------------------

    def key_exists(self, name: str) -> bool:
        target = " ".join(name.lower().split())
        return any(s.name == target for s in self.statements)

    def get(self, name: str) -> str | None:
        """Last wins, matching the C++ map: a duplicate key overwrites."""
        target = " ".join(name.lower().split())
        found = None
        for s in self.statements:
            if s.name == target:
                found = s.value
        return found

    def all_named(self, name: str) -> list[Statement]:
        """Every statement with this key, in order.

        Needed because `layers` blocks legitimately repeat a material name and
        the order is the physical stack; `get()` would silently collapse them.
        """
        target = " ".join(name.lower().split())
        return [s for s in self.statements if s.name == target]

    def section(self, type_: str, label: str | None = None) -> "Section | None":
        """LAST match wins. The C++ keeps named sections in a map keyed on
        type:label, so a duplicate section silently overwrites the earlier
        entry — a first-match here would inspect a section BELFEM ignores.
        (Index-based iteration, `find_all`, still sees every duplicate, which
        is also what the C++ per-index accessors do.)"""
        found = None
        for s in self.sections:
            if s.type == type_ and (label is None or s.label == label):
                found = s
        return found

    def find_all(self, type_: str) -> list["Section"]:
        return [s for s in self.sections if s.type == type_]

    @property
    def path(self) -> str:
        bits = []
        node = self
        while node is not None and node.header_span != (0, 0):
            bits.append(f"{node.type}:{node.label}" if node.label else node.type)
            node = node.parent
        return "/".join(reversed(bits))


class ParseError(Exception):
    def __init__(self, message: str, text: str, offset: int):
        line = text.count("\n", 0, offset) + 1
        col = offset - (text.rfind("\n", 0, offset) + 1) + 1
        super().__init__(f"{message} (line {line}, column {col})")
        self.line, self.column = line, col


def _comment_mask(text: str) -> bytearray:
    """1 where a character is inside a comment, so structure scanning skips it.

    `//` and `#` both truncate to end of line — the second via `clean_string`,
    which `InputFile::remove_comments` calls on every line.
    """
    mask = bytearray(len(text))
    i, n = 0, len(text)
    while i < n:
        ch = text[i]
        if ch == "#" or (ch == "/" and i + 1 < n and text[i + 1] == "/"):
            end = text.find("\n", i)
            end = n if end == -1 else end
            for j in range(i, end):
                mask[j] = 1
            i = end
        else:
            i += 1
    return mask


def parse(text: str) -> Section:
    """Parse a deck into a tree of spans. Raises ParseError on unbalanced input."""
    mask = _comment_mask(text)
    root = Section(span=(0, len(text)), header_span=(0, 0),
                   body_span=(0, len(text)), text=text)

    stack = [root]
    cursor = 0          # start of the pending header / statement text

    for i, ch in enumerate(text):
        if mask[i]:
            continue

        if ch == "{":
            header_start, header_end = cursor, i
            head, dropped, at = header_of(text, mask, header_start, header_end)
            node = Section(
                span=(header_start, i),          # closed out at '}'
                header_span=(header_start, header_end),
                body_span=(i + 1, i + 1),        # closed out at '}'
                text=text,
                parent=stack[-1],
                header_text=head,
                header_discarded=dropped,
                header_offset=at,
            )
            stack[-1].sections.append(node)
            stack.append(node)
            cursor = i + 1

        elif ch == "}":
            if len(stack) == 1:
                raise ParseError("unmatched '}'", text, i)
            node = stack.pop()
            node.body_span = (node.body_span[0], i)
            node.span = (node.span[0], i + 1)
            # Text since the last terminator dies here without a ';' — usually
            # a statement missing its semicolon. The C++ drops it silently;
            # record it so the validator can say so instead.
            dropped = _uncommented(text, mask, cursor, i).strip()
            if dropped:
                node.trailing_discarded = " ".join(dropped.split())
                node.trailing_offset = cursor
            cursor = i + 1

        elif ch == ";":
            stmt_text = text[cursor:i]
            if stmt_text.strip():
                key_end = i
                value_span = None
                colon = _first_colon(text, cursor, i, mask)
                if colon is not None:
                    key_end = colon
                    value_span = (colon + 1, i)
                stack[-1].statements.append(
                    Statement(span=(cursor, i + 1),
                              key_span=(cursor, key_end),
                              value_span=value_span,
                              text=text,
                              mask=mask)
                )
            cursor = i + 1

    if len(stack) != 1:
        raise ParseError(
            f"{len(stack) - 1} unclosed section(s); expected '}}'",
            text, len(text))

    # same rule at end of file: root-level text with no ';' is dropped by the
    # C++ — record it rather than swallow it
    dropped = _uncommented(text, mask, cursor, len(text)).strip()
    if dropped:
        root.trailing_discarded = " ".join(dropped.split())
        root.trailing_offset = cursor

    return root


def _uncommented(text: str, mask: bytearray, start: int, end: int) -> str:
    """The characters in [start, end) that are not inside a comment."""
    return "".join(text[i] for i in range(start, end) if not mask[i])


def header_of(text: str, mask: bytearray, start: int, end: int) -> tuple[str, str, int]:
    """The section header, and whatever preceded it on earlier lines.

    A header span runs from the previous terminator to the `{`, so it can pick
    up a licence banner, a blank line, or stray text. The C++ works on a
    line-split buffer and takes the LAST line before the brace — as the
    reference puts it, "the section header is the line before `{`". This
    matches that, and hands back the discarded remainder so a caller can
    mention it rather than swallow it: `examples/sidecoating` carries a lone
    `/` there, which the parser ignores and which nobody has noticed.
    """
    cleaned = _uncommented(text, mask, start, end)
    lines = [ln.strip() for ln in cleaned.splitlines()]
    for i in range(len(lines) - 1, -1, -1):
        if lines[i]:
            discarded = " ".join(x for x in lines[:i] if x)
            # Offset of the kept line in the ORIGINAL text, so a diagnostic can
            # point at the header a reader sees rather than at wherever the
            # span happens to begin — which is the end of the previous
            # statement, often a line or more earlier.
            hit = text.rfind(lines[i], start, end)
            return lines[i], discarded, (hit if hit != -1 else start)
    return "", "", start


def _first_colon(text: str, start: int, end: int, mask: bytearray) -> int | None:
    """Index of the first non-comment ':' in [start, end), or None.

    First only: `sidesets : 4:9 ;` is one key and one value, and splitting on
    every colon would turn the range into a second key.
    """
    for i in range(start, end):
        if text[i] == ":" and not mask[i]:
            return i
    return None
