# Citation Index Completion

**Date:** 2026-09-21
**Purpose:** Close the remaining gap between the works cited in `src/` and `doc/` and the citation
key index in `doc/literature_references.md`. The index and `scripts/check_citations.py` landed on
2026-09-21; 24 cited works still resolve to nothing, so a reader cannot reach the source behind
them. The mechanism is transcription, not research: most of these works are named with a partial
reference in the module document that cites them, and need the full reference plus a DOI.
**Module:** `doc` (+ `src/physics/gasmodels`, `src/physics/materials`, `src/fem`, `src/sparse`)
**AIs involved:** Claude (index + checker)
**Status:** BLOCKED ON O3 — R1–R5 and R7 done 2026-09-22. **505 citations resolve, 4 do not**, and all four are the O3 works that cannot be identified from the tree. R6 (whether the checker becomes a gate) is now answerable and is the only step left besides O2. O1 resolved 2026-09-22

> **Scope guards:**
> - Out of scope: adding the works to `literature/`. The index records the reference so a reader
>   can find the work independently; holding a copy is a separate decision.
> - Out of scope: rewriting the citations in `src/`. Where a citation is a spelling variant of an
>   indexed work, the **index** gains the alias — the source comment is left alone (minimal rename).
> - A DOI is transcribed from a document that carries it, or left out. It is never composed.

---

## 1. Gap table

| Gap | Where it shows | Closed by |
|---|---|---|
| 24 cited works absent from both authoritative tables | `scripts/check_citations.py` exits 1 | R1–R4 |
| ~~`Gerritsen 1952` in code vs `Gerritsen 1956` in the materials bibliography~~ | ~~`cl_Material_Silver.cpp:37,64`~~ | **closed 2026-09-22** — wrong year label, right DOI |
| Module docs carry full references that the index does not | `src/fem/doc/timestepping_strategy.md` §7 and peers | R5 |

---

## 2. Steps

- [x] **R1 — Gas models (7 works).** **Done 2026-09-22** — five index rows; every DOI was already in the module's own headers, and `Leachman 2009` and `Penoncello & Lemmon 2009` turned out to be the same paper. Original text: `Assael 2011`, `Friend 1989`, `Leachman 2009`,
      `Lemmon & Jacobsen 2004`, `Muzny 2013`, `Penoncello & Lemmon 2009`, `Robinson 1978`.
      Cited in `src/physics/gasmodels/doc/README.md:600-605` and in the EoS headers. These are
      NIST-lineage correlations; the reference-quality papers carry DOIs.
- [◐] **R2 — Maxwell and FEM theory.** **7 of 10 done 2026-09-22** — Nédélec 1980 and 1986, Juntunen & Stenberg, Bíró & Preis, Shahbazi, Hayasaka & Ito (Crossref-confirmed), Zaglmayr and Schoof & Yarberry (transcribed from `nedelec_derivation.md:487-489`). Herrmann 1978, Kuczmann & Iványi 2008 and Meunier 2008 moved to O3. Original text: `Bíró & Preis 1989`, `Hayasaka & Ito 2019`,
      `Herrmann 1978`, `Juntunen & Stenberg 2009`, `Kuczmann & Iványi 2008`, `Meunier 2008`,
      `Nédélec 1980`, `Shahbazi 2005`, `Zaglmayr 2006`, `Schoof & Yarberry 1994`.
      Nédélec 1980 is the foundational edge-element paper and should arguably also be routed in
      CLAUDE.md's literature table. (after: none)
- [x] **R3 — Solvers and time stepping (3 works).** **Done 2026-09-22** — all three Crossref-confirmed; the MUMPS and STRUMPACK DOIs cross-check against references [24] and [26] of Messe et al. 2023. Original text: `Duff & L'Excellent 2000` and `Ghysels 2016`
      (`src/sparse/doc/solver_memory_and_compression.md:281`), `Grigorieff 1983`
      (`src/fem/doc/timestepping_strategy.md:46`, `cl_FEM_Controller.cpp:2890`).
- [x] **R4 — Materials (3 works).** **Done 2026-09-22** — three rows in the materials bibliography, all Crossref-confirmed. Original text: `Giebe & Blechschmidt 1931`, `Koester & Franz 1961`,
      `Zou & Balandin 2001`. These belong in the **materials bibliography**, not the general
      index. ~~and the `Gerritsen` year conflict in O1~~ — closed 2026-09-22.
- [x] **R5 — Sweep the module reference sections.** **Done 2026-09-22, absorbed into R7** — every module reference list that holds an unresolved work was read and harvested. Original text: Several module docs keep their own reference
      list (`timestepping_strategy.md` §7, `anderson_acceleration_theory.md` §6,
      `ghost_penalty_stabilization.md`, `thin_cut_nonunit_rectification.md`,
      `cohomology_algorithms.md`). Those already-transcribed references were harvested into the
      index on 2026-09-21; re-run the sweep after R1–R4 so no module list holds a work the index
      does not. (after: R1, R2, R3, R4)
- [x] **R7 — Harvest the parenthesized-citation population.** **Done 2026-09-22** — 22 index rows plus two aliases (`Sousa 2022b` → Alves 2022b, `Zunino 2009` → Ern et al. 2009). Six DOIs came from Crossref, five were already in the gasmodels and Maxwell headers and each was corroborated by an independent bibliographic entry elsewhere in the tree; nine are books with no DOI, carried with edition and publisher. Three matcher defects fixed in passing: a possessive (`Nédélec's 1980`), an initials fragment (`J.-C 1980`), and citations inside fenced blocks, which quote data files rather than cite works (`SVEHLA (1994)` in a `trans.inp` excerpt). Original text: Widening the matcher to accept
      `Author (YYYY)` exposed ~26 further works, nearly all in module reference lists
      (`gasmodels/doc/gasmodels_usage_guide.md:1139-1162`, `math/graph/doc/README.md:203-209`,
      `fem/interpolation/doc/interpolation_usage_guide.md:1877-1878`,
      `fem/iwg/doc/README.md:260-261`). Those lists already carry title, journal, volume and pages;
      several DOIs are in the gasmodels headers but in different files from the list that names the
      work (`10.1016/0009-2509(72)80096-4` for Soave 1972, `10.1063/1.1349047` for Span 2000,
      `10.1016/0378-3812(85)87016-3` for Schmidt & Wagner 1985, `10.1063/1.555898` for
      Setzmann & Wagner 1991). Mostly transcription; the residue needs Crossref. (after: O3)
- [ ] **R6 — Wire the checker into the routine.** *Now answerable: the checker is at 4, all blocked on O3.* Decide whether `check_citations.py` joins
      `check_doc_claims.py` in the session checklist, in the post-commit auto-review, or stays
      on-demand. (after: R5)

---

## 3. Open questions

- [x] **O1 — `Gerritsen 1952` or `Gerritsen 1956`?** **RESOLVED 2026-09-22 → 1956; the code's
      year label was wrong, its DOI was right** (Christian). There are three Gerritsen papers in
      *Physica*, all confirmed against Crossref: the 1952 electrical resistances of noble/transition
      alloys with Linde (10.1016/S0031-8914(52)80222-8), the 1953 magnetoresistances
      (10.1016/S0031-8914(53)80007-8), and the 1956 *Thermal conductivity of some dilute silver
      alloys* with Linde (10.1016/S0031-8914(56)90035-0). The last is what
      `set_lambda_coefficients` fits, and `sample AG2T, RRR 407` is one of its dilute-alloy
      samples. The two year labels in `cl_Material_Silver.cpp:37,64` were corrected; the
      bibliography row gained Linde as second author. The 1952 and 1953 papers are **not** indexed
      — nothing in the tree cites them (silver resistivity comes from Matula, magnetoresistance
      from Lüthi, Strom-Olsen, Stout and Iwasa).
- [ ] **O3 — Four works cannot be identified from the tree.** Each is cited by surname and year
      with no fuller reference anywhere, and composing one from recollection is exactly what the
      `literature` skill forbids:
      - `Herrmann 1978` (`postprocessor_recovery_theory.md:598`) — carries a superconvergence
        theorem for derivatives at Gauss points. The result is usually attributed to Barlow 1976,
        so the year and the name together may be a mis-citation rather than a findable work.
      - `Kuczmann & Iványi 2008` and `Meunier 2008` (`maxwell_weak_forms.md:367`) — the document
        calls them "general Maxwell FEM references **from the notes**", so the source is somebody's
        lecture notes rather than the tree.
      - `Robinson 1978` (`cl_GM_EoS_AlphaFunctionFactory.hpp:40`) — the 1978 revision of the
        Peng-Robinson m(ω) correlation, unambiguous as physics (`AlphaFunctionFactory.cpp:53-57`
        implements the ω ≤ 0.491 switch), but the bibliographic form is not in the tree. It is
        usually a GPA research report, which has no DOI.
- [ ] **O2 — Does the index want an `unheld: true` marker?** A reader cannot currently tell an
      entry whose `.txt` is in `literature/` from one that is only a reference. Cheap to add,
      and it would let the skill answer "can I read this one" without an `ls`.

---

## 4. Residue as of 2026-09-22 (widened matcher)

Live list: `scripts/check_citations.py`. Snapshot:

```
  Herrmann 1978                      src/fem/maxwell/doc/postprocessor_recovery_theory.md:598
  Kuczmann & Iványi 2008             src/fem/maxwell/doc/maxwell_weak_forms.md:367
  Meunier 2008                       src/fem/maxwell/doc/maxwell_weak_forms.md:367
  Robinson 1978                      src/physics/gasmodels/cl_GM_EoS_AlphaFunctionFactory.hpp:40
```
