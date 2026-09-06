"""A deck you can read, and later edit without reflowing it.

`Document` owns the original bytes and a span tree over them. Serialising an
untouched document returns the original bytes unchanged — not "formatted the
same way", but the same bytes. Edits are recorded as span replacements and
applied at serialise time, so everything outside an edited span is carried
through verbatim.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from .parser import ParseError, Section, parse


def read_verbatim(path) -> str:
    """Read without newline translation.

    `newline=""` matters: Python's universal-newline mode turns CRLF into LF on
    read, so a CRLF deck would come back rewritten on save and the round trip
    would not be byte-identical. Passed to open() rather than Path.read_text()
    because that keyword only exists from Python 3.13 and this must run on the
    3.9 in the SCLS toolchain.
    """
    with open(path, encoding="utf-8", newline="") as handle:
        return handle.read()


def write_verbatim(path, text: str) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        handle.write(text)


@dataclass(order=True)
class _Edit:
    start: int
    end: int
    replacement: str


@dataclass
class Document:
    text: str
    root: Section
    path: Path | None = None
    _edits: list[_Edit] = field(default_factory=list)

    @classmethod
    def load(cls, path: str | Path) -> "Document":
        p = Path(path)
        return cls.loads(read_verbatim(p), path=p)

    @classmethod
    def loads(cls, text: str, path: Path | None = None) -> "Document":
        return cls(text=text, root=parse(text), path=path)

    # -- writing ----------------------------------------------------------

    @property
    def modified(self) -> bool:
        return bool(self._edits)

    def replace_value(self, statement, new_value: str) -> None:
        """Rewrite one statement's value, leaving its key and layout alone."""
        if statement.value_span is None:
            raise ValueError(
                f"{statement.key!r} is a flag and has no value span; "
                "rewriting it would have to invent syntax"
            )
        self._edits.append(_Edit(*statement.value_span, replacement=new_value))

    def dumps(self) -> str:
        if not self._edits:
            return self.text                      # the byte-identical path

        edits = sorted(self._edits)
        for a, b in zip(edits, edits[1:]):
            if b.start < a.end:
                raise ValueError("overlapping edits")

        out, cursor = [], 0
        for e in edits:
            out.append(self.text[cursor:e.start])
            out.append(e.replacement)
            cursor = e.end
        out.append(self.text[cursor:])
        return "".join(out)

    def save(self, path: str | Path | None = None) -> Path:
        target = Path(path) if path else self.path
        if target is None:
            raise ValueError("no path to save to")
        write_verbatim(target, self.dumps())
        return target

    # -- reading ----------------------------------------------------------

    def section(self, type_: str, label: str | None = None) -> Section | None:
        return self.root.section(type_, label)

    def sections(self, type_: str) -> list[Section]:
        return self.root.find_all(type_)

    def walk(self):
        """Every section in the document, depth first."""
        stack = list(self.root.sections)
        while stack:
            node = stack.pop(0)
            yield node
            stack = list(node.sections) + stack

    def layer_stack(self, tape: str) -> list[tuple[str, str]]:
        """The `layers : <tape>` stack as ordered (material, thickness) pairs.

        Read positionally and in order, with duplicates kept, because that is
        what the C++ does: `read_thin_shell_data` word-splits the raw lines,
        repeated material names are legal, and the order IS the physical stack
        bottom to top. Going through a dict here would silently collapse a
        deck like copper / ybco / copper into two layers.
        """
        section = self.root.section("layers", tape.lower())
        if section is None:
            return []
        return [(s.key, s.value) for s in section.statements]


__all__ = ["Document", "ParseError", "Section"]
