# Material property literature sources {#physics_materials_material_property_sources}

**Date:** 2026-08-28
**Purpose:** One table per property family naming the literature behind every material fit, with the full bibliography
**Module:** src/physics/materials

The material classes identify their data sources in comments within their
`create_*()` routines. This document collects those statements in one place
and presents them in a standardized form. The `file:line` positions in the
class sources are authoritative. When a comment gives only a vague attribution,
such as "literature value" or "experimental data", the table reports that
wording instead of inventing a source. Model-level references for the E-J power
law are in `resistivity_laws.md`; sources for homogenization are in
`alloy_homogenization.md`; and sources for the Callaway model are in
`callaway_thermal_conductivity.md`.

Property symbols: `alpha` denotes thermal expansion, `cp` specific heat, `rho`
electrical resistivity (the quantity stored by the `rho` group of a generated
database), `E` Young's modulus, and `nu` Poisson ratio.

## Main table

| Material | alpha | cp | rho | E | nu |
|---|---|---|---|---|---|
| Copper | Touloukian TPRC | Touloukian TPRC | Matula 1979 (Bloch-Grueneisen, 273.15 K anchor); Kohler: De Launay 1959, Benz 1969, Arentz 1982, Clausecker 1969, Strom-Olsen 1967 | Wachtman fit vs Blanke 1989 | copper.org / Wolfram Cloud (RT value) |
| Aluminum | Touloukian TPRC | Touloukian TPRC | Hust-Lankford 1984 anchor; Debye curve: Desai 1984, Cook 1975; Kohler: Luethi 1960, Fickett 1972 | Wachtman fit vs Blanke 1989 | Grueneisen-consistent, RT value Wolfram Cloud |
| Chromium | e3sconf 2024 (primary), Touloukian TPRC (secondary) | unnamed literature spread (sanity: cp(293 K) = 451 vs 449) | White-Woods 1959 anchor; Debye: Anderson 1970; Kohler: Kozlova-Kondorskii 1963, Arajs-Dunmyre 1965 | Armstrong-Brown 1964 (scaled to 279 GPa at RT) | Wolfram Cloud (RT value) |
| Silver | Touloukian TPRC | Touloukian TPRC | Matula 1979 (273.15 K anchor); Kohler: Luethi 1960, Strom-Olsen 1967, Stout 1939, Iwasa 1993, Li 2015 | Wachtman fit vs Blanke 1989 | Wolfram Cloud (RT value) |
| Indium | Touloukian TPRC (6..374 K, rms 0.005%) | Touloukian TPRC; low-T gamma/beta: TPRC curve 3 | White-Sondheimer via PhysRevB 23 3845 anchor; Kohler: Luethi 1960, PhysRev 120 1167 | Kim-Ledbetter 1998 | Kim-Ledbetter 1998 |
| White tin | Touloukian TPRC | Touloukian TPRC; Sommerfeld/Debye: O'Neil 1965 | White 1968 via Hariharan 1979 anchor; Debye high-T: Meaden 1965 via Hall 1968; Kohler: Luethi 1960 | Wachtman fit vs Blanke 1989 | Wolfram Cloud (RT value) |
| Lead | Touloukian TPRC | Touloukian TPRC | White 1968 via Hariharan 1979 anchor; Debye high-T: Meaden 1965 via Hall 1968; Kohler: Luethi 1960 | Wachtman fit vs Blanke 1989 | Wolfram Cloud (RT value) |
| Iron | Touloukian TPRC | Touloukian TPRC (Curie anomaly smoothed) | White-Woods 1959 (T < 295 K, RRR 40..104), Arajs-Colvin 1964 (300..1300 K), magnetization: Crangle-Goodman 1971; Kohler: Klaffky-Coleman 1974, Luethi 1960; valid to 860 K | Wachtman fit vs Blanke 1989 | Wolfram Cloud (RT value) |
| Nickel | Touloukian TPRC | Touloukian TPRC; gamma/beta: TPPM | White-Woods 1959, Farrell-Greig 1968; magnetization: Crangle-Goodman 1971; Kohler: Luethi 1960 | Wachtman fit vs Blanke 1989 (two curves across the Curie point) | Wolfram Cloud (RT value) |
| Magnesia (MgO) | Simon 1994, Durand 1936 | Touloukian TPRC curve 4 + Simon 1994 | n/a (insulator) | Blanke 1989 (Simon 1994 judged too low for HTS tapes) | Simon 1994 |
| Hastelloy C-276 | Lu 2008 | Lu 2008 (15..120 K); RT Debye tail anchored at unnamed literature cp = 427 | Lu 2008 (log-log above 10.73 K) | Cryogenics 2006 + MATWEB | Cryogenics 2006 + MATWEB |
| YBCO | Salomons 1987 | Baak 1989 (T < 20 K), Lang 1988 (50..350 K) | Sommerfeld 2003 (Bloch-Grueneisen, RRR ~ 50) | Lei-Ledbetter 1991 (normalized curve, user-scaled) | Lei-Ledbetter 1991 |
| Alloy (homogenized) | see alloy_homogenization.md | see alloy_homogenization.md | see alloy_homogenization.md (Nordheim, Bruggeman, Hust-Lankford) | see alloy_homogenization.md (Eshelby, Hill, Voigt/Reuss) | see alloy_homogenization.md |

The supporting property thermal conductivity, `lambda`, is included here
because users often ask for it next. Copper, Aluminum, and Iron follow
Hust-Lankford 1984. Silver is fitted to Mendelssohn 1952, Gerritsen 1956, and
Li 2015; the class header notes that the RRR data are sparse. Indium, White
tin, Lead, Magnesia, and Nickel follow Touloukian TPRC. YBCO uses the Callaway
model fitted to Sommerfeld 2003, which the code flags as a degenerate interim
fit. Hastelloy has no lambda source and uses only the Wiedemann-Franz
construction.

## Bibliography

| Key | Authors | Year | Journal / Publisher | DOI or identifier |
|---|---|---|---|---|
| Touloukian TPRC | Y. S. Touloukian et al. | 1970 ff. | Thermophysical Properties of Matter (TPRC data series), IFI/Plenum | series; cited per curve in-code |
| TPPM | Thermophysical Properties of Matter, vol. 1 | 1970 | IFI/Plenum | volume/page cited in-code |
| Matula 1979 | R. A. Matula | 1979 | J. Phys. Chem. Ref. Data 8, 1147 | 10.1063/1.555614 |
| Blanke 1989 | W. Blanke (ed.) | 1989 | Thermophysikalische Stoffgroessen, Springer | ISBN book |
| Hust-Lankford 1984 | J. G. Hust, A. B. Lankford | 1984 | NBSIR 84-3007 / NBS SP 260-90 | 10.6028/nbs.ir.84-3007 |
| White-Woods 1959 | G. K. White, S. B. Woods | 1959 | Phil. Trans. R. Soc. A 251, 273 | 10.1098/rsta.1959.0004 |
| De Launay 1959 | J. De Launay et al. | 1959 | J. Phys. Chem. Solids | 10.1016/0022-3697(59)90038-1 |
| Benz 1969 | M. G. Benz | 1969 | J. Appl. Phys. 40 | 10.1063/1.1657896 |
| Arentz 1982 | R. Arentz et al. | 1982 | Phys. Rev. B 26, 2727 | 10.1103/PhysRevB.26.2727 |
| Clausecker 1969 | K. Clausecker | 1969 | Z. Physik | 10.1007/BF02422537 |
| Strom-Olsen 1967 | J. O. Strom-Olsen | 1967 | Proc. R. Soc. A | jstor.org/stable/2415889 |
| Desai 1984 | P. D. Desai et al. | 1984 | J. Phys. Chem. Ref. Data | 10.1063/1.555725 |
| Cook 1975 | J. G. Cook et al. | 1975 | ORNL-5079 report | ORNL-5079 |
| Luethi 1960 | B. Luethi | 1960 | Dissertation, ETH Zuerich (Widerstandsaenderung von Metallen in hohen Magnetfeldern) | thesis |
| Fickett 1972 | F. R. Fickett | 1972 | FNAL report (Magnetoresistivity of Copper and Aluminum at Cryogenic Temperatures) | report |
| Armstrong-Brown 1964 | P. E. Armstrong, H. L. Brown | 1964 | Trans. AIME | as cited in-code |
| Kozlova-Kondorskii 1963 | N. Kozlova, E. Kondorskii | 1963 | Soviet Physics JETP | as cited in-code |
| Arajs-Dunmyre 1965 | S. Arajs, G. R. Dunmyre | 1965 | J. Appl. Phys. | 10.1063/1.1703039 |
| Anderson 1970 | Anderson, Stewart, Ramsay | 1970 | phys. status solidi (b) | 10.1002/pssb.19700370137 |
| e3sconf 2024 | (chromium expansion dataset) | 2024 | E3S Web of Conferences | 10.1051/e3sconf/202459202009 |
| Smith-Fickett 1995 | D. R. Smith, F. R. Fickett | 1995 | J. Res. NIST 100 | 10.6028/jres.100.012 |
| Mendelssohn 1952 | K. Mendelssohn et al. | 1952 | Proc. Phys. Soc. A 65 | 10.1088/0370-1298/65/6/301 |
| Gerritsen 1956 | A. N. Gerritsen | 1956 | Physica | 10.1016/S0031-8914(56)90035-0 |
| Li 2015 | Li et al. | 2015 | IOP Conf. Ser.: Mater. Sci. Eng. 102 | 10.1088/1757-899X/102/1/012027 |
| Stout 1939 | J. W. Stout et al. | 1939 | J. Am. Chem. Soc. | 10.1021/ja01871a006 |
| Iwasa 1993 | Y. Iwasa et al. | 1993 | Cryogenics | 10.1016/0011-2275(93)90199-X |
| Kim-Ledbetter 1998 | S. Kim, H. Ledbetter | 1998 | Mater. Sci. Eng. A | 10.1016/S0921-5093(98)00490-0 |
| O'Neil 1965 | H. R. O'Neal et al. | 1965 | Phys. Rev. 137, A748 (in-code cited as 1964) | 10.1103/PhysRev.137.A748 |
| White 1968 | G. K. White | 1968 | Experimental Techniques in Low Temperature Physics (book) | via Hariharan 1979 |
| Hariharan 1979 | Hariharan et al. | 1979 | Pramana 13, 117 | 10.1007/BF02872130 |
| Meaden 1965 | G. T. Meaden | 1965 | Electrical Resistance of Metals (book) | via Hall, NBS TN 365 (1968) |
| Crangle-Goodman 1971 | J. Crangle, G. M. Goodman | 1971 | Proc. R. Soc. A | 10.1098/rspa.1971.0044 |
| Arajs-Colvin 1964 | S. Arajs, R. V. Colvin | 1964 | phys. status solidi | 10.1002/pssb.19640060317 |
| Klaffky-Coleman 1974 | R. W. Klaffky, R. V. Coleman | 1974 | Phys. Rev. B 10, 2915 | 10.1103/PhysRevB.10.2915 |
| Farrell-Greig 1968 | T. Farrell, D. Greig | 1968 | J. Phys. C 1 | 10.1088/0022-3719/1/5/326 |
| Simon 1994 | N. J. Simon | 1994 | NISTIR 5030 | report |
| Durand 1936 | M. A. Durand | 1936 | Physics 7 | 10.1063/1.1745396 |
| Lu 2008 | J. Lu, E. S. Choi, H. D. Zhou | 2008 | J. Appl. Phys. 103 (Hastelloy C-276 at cryogenic temperatures) | 10.1063/1.2899058 |
| Cryogenics 2006 | (Hastelloy mechanical dataset) | 2006 | Cryogenics 46 | 10.1016/j.cryogenics.2006.01.014 |
| Cryogenics 2023 | (Hastelloy composition) | 2023 | Cryogenics | 10.1016/j.cryogenics.2023.103776 |
| Salomons 1987 | E. Salomons et al. | 1987 | Physica B+C | 10.1016/0378-4363(87)90092-1 |
| Baak 1989 | J. Baak et al. | 1989 | Physica C | 10.1016/0921-4534(89)91125-8 |
| Lang 1988 | M. Lang et al. | 1988 | Z. Phys. B | 10.1007/BF01312506 |
| Lei-Ledbetter 1991 | M. Lei, H. Ledbetter | 1991 | NISTIR 3980, Fig. 4.3 | report |
| Sommerfeld 2003 | Sommerfeld et al. | 2003 | Phys. Rev. B 67, 174520 | 10.1103/PhysRevB.67.174520 |
| CRC Handbook | CRC Handbook of Chemistry and Physics | n/a | CRC Press (isotope masses and abundances) | edition not stated in-code |

## Known gaps (as found in the sources, 2026-08-28)

- The mass density `ref_density` is an uncited literal for every material
  except Silver, which cites Smith-Fickett 1995.
- Hastelloy C-276 has no data source for `lambda`; it uses only the
  Wiedemann-Franz construction.
- The `Alloy` class contains no in-code citations; its sources are documented
  exclusively in `alloy_homogenization.md`.
- The iron and nickel spin-disorder resistivity amplitudes are described as
  "typical values" without a source.
- The Debye curve for Indium and the Debye temperature of 275 K for YBCO are
  bare fits.
- "Wolfram Cloud" is the room-temperature Poisson source for seven metals,
  although it is a lookup service rather than a primary reference.
- Base-class model names such as Bloch-Grueneisen, Matthiessen, Kohler,
  Wachtman, Pippard, and Wiedemann-Franz are used throughout, but none has a
  bibliographic entry anywhere in the tree.

## Citation defects in the sources (typos worth a cleanup pass)

- In `cl_Material_Iron.cpp:168` and
  `cl_Material_Nickel.cpp:145,149`, the DOI lacks its leading `1`
  (`0.1098/rsta.1959.0004`).
- In `cl_Material_Silver.cpp:99`, the DOI `0.6028/jres.100.012` lacks its
  leading `1`.
- In `cl_Material_Magnesia.cpp:148`, "imon, 1994" should read "Simon"; the
  comment at `:114` says "poisson" above code that builds E.
- In `cl_Material_Iron.cpp:100`, "Arays" should read "Arajs".
- In `cl_Material_Lead.cpp:93` and `cl_Material_WhiteTin.cpp:84`,
  "referene" is misspelled.
- The statement "fitted against experimental data" in
  `cl_Material_Silver.cpp:33` contradicts the Blanke 1989 / Wolfram Cloud
  attribution at `:80`.
- The citation in `powerlaws.hpp:44` gives the Messe et al. 2023 paper a title
  that does not match the one in `doc/literature_references.md`; one of the
  two is wrong.
- The in-code citation "O'Neil et al, 1964" for Phys. Rev. 137, A748 gives the
  wrong year; the volume was published in 1965.
- The in-code citation "Gerritsen 1952" is paired with a 1956 Physica DOI.

Line numbers above are advisory; grep for the citation text when in doubt.
