# DR-100 fixed: the orphaned dummy send that poisoned the comm fabric

**Date:** 2026-08-27
**Topic:** `examples/RLC_Circuit` np=4 and the gantry np=10 `MPI_ERR_TRUNCATE` crashes —
one two-line defect, one month old, found by a probe ladder and killed by a deletion.

## The defect

`Postprocessor::recover_fields()` has an early return for workers that own no
postprocessor nodes. Since `b88d9d37` ( 2026-07-19 ) that path sent a dummy `Matrix`
**twice** — the second send tagged `[DIAG]`, claiming to pair with a "patch_count gather
below" that never existed anywhere in the tree. Root collects exactly once, so every save
step left one unmatched 3-`index_t` size header ( 12 bytes in default builds ) queued on
the worker's base tag.

The comm fabric gives each rank pair exactly two tags and relies purely on execution
order, so the stray never errored where it was born: it sat in the queue until the next
1-element size collect on the same ( source, tag ) matched it — `MPI_ERR_TRUNCATE` in
whatever innocent collect ran next. RLC np=4 died in `collect_matrices` under `save_IV`;
gantry np=10 died one exchange earlier, in `synch_node_indices`. Same producer, different
victims — which is exactly why the crash site misled for so long.

## Why some rank counts never crash

The dummy is only sent by ranks that own **zero** postprocessor nodes, so the trigger is
partition layout, not load or timing: np=2 RLC and np=8 gantry draw layouts where every
rank owns nodes, np=4 / np=10 do not. The 2026-08-26 "latent timing race" reclassification
in the register was wrong and is retracted — the defect is fully deterministic per layout.

## What died in audit before the fix ( recorded so nobody re-walks it )

- My "field-count handshake" producer theory: the early 4-byte counts from
  `FieldData::distribute` are benign — workers block immediately after sending them
  ( Codex, by source order ).
- My "rank-divergent postprocessor list" theory: a brace-scan misread; the Maxwell
  factory creates postprocessors on all ranks, decisions max-reduced and redistributed.
- The pre-registered C1 ( barrier at the top of `DofManager::postprocess` ): dropped as
  pairing the benign handshake only. Both auditors independently converged on the dummy
  as the true producer, matching the probe-ladder amendment.

## The fix ( commit this session )

- **`cl_FEM_Postprocessor.cpp`**: deleted the orphan second send and its barrier; the
  early return now mirrors the normal path message-for-message. The site comment records
  the history.
- **`commtools.cpp`**: the two-tags-per-pair ordering contract is documented on
  `comm_tag()` — the load-bearing rule is that any conditional early-out must reproduce
  the communication pattern of the path it replaces, message for message.
- All DR-100 probes stripped, including the committed per-request Waitall decode in
  `collect( Vector )` ( restored to the pre-probe one-liner ). The `[DIAG]` `_Volumes`
  weight guard stays — it belongs to the slave-surface seam campaign.

## Gates

- **G1** — RLC np=4: red baseline aborts signal 6 at the *first* save; fixed binary ran
  the full transient, 210+ saves, zero MPI errors. The drain probes were still in that
  binary and reported clean queues at every exchange boundary — the fabric is now
  actually empty between exchanges, not merely surviving.
- **G2** — gantry np=10 ( probe-stripped binary, private build tree ): first exodus save,
  `save_IV`, and onward iteration all clean at the layout that used to die.
- **G3** — `make check-fast` 9/9.
- Audits: Codex **accept** ( one comment-accuracy finding — the header is 3-`index_t`,
  12 bytes only in default builds — applied ), Grok **accept**.

## Residue

- DR-113 filed ( P3, deferred ): the fabric has no structural defense against this defect
  class; options ( debug drain check at chokepoints, per-exchange tag separation ) are
  recorded in the register. Wrong week for a fabric change.
- Round record: `tmp/ai_exchange/dr100_handshake_ordering.md`.
