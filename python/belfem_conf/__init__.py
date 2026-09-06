"""Tooling for the BELFEM `input.conf` contract.

The contract has two artifacts that must agree: `doc/input_file_reference.md`
(for people) and `doc/input_schema.yaml` (for programs). This package is what
checks that they still agree with the C++ that actually parses the deck.

Three commands: `drift` (schema against the C++ parse sites), `check` (deck
validation against the schema), and `roundtrip` (the lossless-parser gate).
The configurator GUI is planned on the same schema — see
`todo/input_conf_configurator_plan.md`.
"""

__version__ = "0.2.0"
