# Devlog 2026-04-21 — Disabling CLion Auto-Import for C++

**Date:** 2026-04-21
**Topic:** Disabling automatic `#include` additions in CLion
**AIs involved:** Claude
**Claude Confidence:** high

## Summary

The user requested help to disable CLion's feature that automatically adds `#include` statements from third-party libraries. Investigation of the project's `.idea` configuration confirmed that these settings are managed at the IDE level rather than within the project files.

## Key Findings

- CLion auto-import settings for C++ are not found in the `.idea` directory, implying they are IDE-wide.
- The relevant settings are located under `Editor | General | Auto Import | C/C++`.

## Changes Made / Proposed

- Provided instructions to the user to disable `Show import popup` and `Insert imports on completion`.
- Identified the "Exclude from auto-import and completion" list as an alternative for targeting specific libraries.

## Files Updated

- `todo/ai_exchange.md` (documented finding)
