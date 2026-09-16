# Material property literature sources {#physics_materials_material_property_sources}

**Date:** 2026-09-15 (first version 2026-08-28)
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
database), `E` Young's modulus, and `nu` Poisson's ratio.

The nine pure metals represent isotropic polycrystals. Their E and nu are obtained from
the quasi-harmonic closure of `Metal::create_mech` (bulk and shear modulus as
exponentials of the logarithmic thermal strain; Garai-Laugier 2007 derive the
form for the bulk modulus, BELFEM applies it to the shear modulus as well), fit
over 0..300 K to the elastic data named in the E column, which span 4 K to room
temperature; data from single-crystal sources are Hill-averaged, ultrasonic (adiabatic) bulk moduli
are converted to isothermal values with the material's own alpha, cp and rho
before the fit. The resulting served moduli are isothermal (dynamic
moduli, not static ones). Blanke 1989 is a compilation without sources or
data points; it is not used for any level or slope, only as the citable
reference for the temperature ceiling: `T_max` of copper, silver, chromium
and white tin is set where the quasi-harmonic E(T) departs from Blanke's curve
by more than 5 % in shape. In the overlap below room temperature Blanke lies
2..4 % above the ultrasonic literature for Cu, Al, Ag and within 1 % for Fe;
its lead value (16 GPa at 300 K) is a static-type modulus, 33 % below the
dynamic Hill average, and is the one exception where Blanke sets a level.

## Main table

| Material | alpha | cp | rho | E | nu |
|---|---|---|---|---|---|
| Copper | Touloukian TPRC | Touloukian TPRC | Matula 1979 (Bloch-Grueneisen, 273.15 K anchor); Kohler: De Launay 1959, Benz 1969, Arentz 1982, Clausecker 1969, Strom-Olsen 1967 | Ledbetter 1981 (polycrystalline material, 5..295 K); quasi-harmonic fit | Ledbetter 1981 |
| Aluminum | Touloukian TPRC | Touloukian TPRC | Hust-Lankford 1984 anchor; Debye curve: Desai 1984, Cook 1975; Kohler: Luethi 1960, Fickett 1972 | Kamm-Alers 1964 (single-crystal data, Hill average); quasi-harmonic fit | Kamm-Alers 1964 |
| Chromium | e3sconf 2024 (primary), Touloukian TPRC (secondary) | unnamed literature spread (sanity: cp(293 K) = 451 vs 449) | White-Woods 1959 anchor; Debye: Anderson 1970; Kohler: Kozlova-Kondorskii 1963, Arajs-Dunmyre 1965 | Palmer-Lee 1971 cryogenic plateau (single-crystal data, Hill average); a common softening constant; the spin-density-wave anomalies are not represented | Palmer-Lee 1971, held at 0.2371 |
| Silver | Touloukian TPRC | Touloukian TPRC | Matula 1979 (273.15 K anchor); Kohler: Luethi 1960, Strom-Olsen 1967, Stout 1939, Iwasa 1993, Li 2015 | Neighbours-Alers 1958 (single-crystal data, Hill average); quasi-harmonic fit | Neighbours-Alers 1958 |
| Indium | Touloukian TPRC (6..374 K, rms 0.005%) | Touloukian TPRC; low-T gamma/beta: TPRC curve 3 | White-Sondheimer via PhysRevB 23 3845 anchor; Kohler: Luethi 1960, PhysRev 120 1167 | Kim-Ledbetter 1998 (polycrystalline Varshni-fit data); quasi-harmonic fit | Kim-Ledbetter 1998 |
| White tin | Touloukian TPRC | Touloukian TPRC; Sommerfeld/Debye: O'Neil 1965 | White 1968 via Hariharan 1979 anchor; Debye high-T: Meaden 1965 via Hall 1968; Kohler: Luethi 1960 | Rayne-Chandrasekhar 1960 (tetragonal single-crystal data at 4.2, 77, 300 K, Voigt-Reuss-Hill); quasi-harmonic fit | Rayne-Chandrasekhar 1960 |
| Lead | Touloukian TPRC | Touloukian TPRC | White 1968 via Hariharan 1979 anchor; Debye high-T: Meaden 1965 via Hall 1968; Kohler: Luethi 1960 | K from Waldorf-Alers 1962 (single-crystal data); E level and shape from Blanke 1989 (static-type curve, 95..300 K) - static moduli by decision, see the class header | derived from K and E: 0.435 at 300 K (handbook 0.44) |
| Iron | Touloukian TPRC | Touloukian TPRC (Curie anomaly smoothed) | White-Woods 1959 (T < 295 K, RRR 40..104), Arajs-Colvin 1964 (300..1300 K), magnetization: Crangle-Goodman 1971; Kohler: Klaffky-Coleman 1974, Luethi 1960; valid to 860 K | Rayne-Chandrasekhar 1961 (single-crystal data, zero field, Hill average); quasi-harmonic fit | Rayne-Chandrasekhar 1961 (isothermal nu flat) |
| Nickel | Touloukian TPRC | Touloukian TPRC; gamma/beta: TPPM | White-Woods 1959, Farrell-Greig 1968; magnetization: Crangle-Goodman 1971; Kohler: Luethi 1960 | Alers-Neighbours-Sato 1960 (single-crystal data at 10 kOe, Hill average); quasi-harmonic fit; the demagnetized Delta-E dip is not represented | Alers-Neighbours-Sato 1960 |
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
| Ledbetter 1981 | H. M. Ledbetter | 1981 | phys. stat. sol. (a) 66, 477 (polycrystalline copper, 5..295 K) | 10.1002/pssa.2210660209 |
| Ledbetter-Naimon 1974 | H. M. Ledbetter, E. R. Naimon | 1974 | J. Phys. Chem. Ref. Data 3, 897 (copper review) | 10.1063/1.3253150 |
| Ledbetter-Reed 1973 | H. M. Ledbetter, R. P. Reed | 1973 | J. Phys. Chem. Ref. Data 2, 531 (iron, nickel, Fe-Ni; section 15 on the Delta-E effect) | 10.1063/1.3253127 |
| Kamm-Alers 1964 | G. N. Kamm, G. A. Alers | 1964 | J. Appl. Phys. 35, 327 (aluminum single crystal, 0..300 K) | 10.1063/1.1713309 |
| Neighbours-Alers 1958 | J. R. Neighbours, G. A. Alers | 1958 | Phys. Rev. 111, 707 (silver and gold single crystals, 0..300 K) | 10.1103/PhysRev.111.707 |
| Waldorf-Alers 1962 | D. L. Waldorf, G. A. Alers | 1962 | J. Appl. Phys. 33, 3266 (lead single crystal, 0..300 K) | 10.1063/1.1931149 |
| Rayne-Chandrasekhar 1960 | J. A. Rayne, B. S. Chandrasekhar | 1960 | Phys. Rev. 120, 1658 (beta tin single crystal, 4.2..300 K) | 10.1103/PhysRev.120.1658 |
| Rayne-Chandrasekhar 1961 | J. A. Rayne, B. S. Chandrasekhar | 1961 | Phys. Rev. 122, 1714 (iron single crystal, 4.2..300 K) | 10.1103/PhysRev.122.1714 |
| Alers-Neighbours-Sato 1960 | G. A. Alers, J. R. Neighbours, H. Sato | 1960 | J. Phys. Chem. Solids 13, 40 (nickel single crystal at 10 kOe, 0..760 K) | 10.1016/0022-3697(60)90125-6 |
| Palmer-Lee 1971 | S. B. Palmer, E. W. Lee | 1971 | Phil. Mag. 24, 311 (chromium single crystal, 4.2..345 K, zero field) | 10.1080/14786437108227390 |
| Garai-Laugier 2007 | J. Garai, A. Laugier | 2007 | J. Appl. Phys. 101, 023514 (isothermal bulk modulus as an exponential of the integrated thermal expansion) | 10.1063/1.2424535 |
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
- Lead is the one metal served with static-level moduli (E, G, nu from
  Blanke's static-type curve with the crystal's bulk modulus); the dynamic
  Hill average of the same crystal is 24 GPa against the static 16 GPa at
  300 K, because lead's shear modulus relaxes strongly at low frequency.
- White tin's two softening constants are based on three temperatures.
- Chromium's spin-density-wave anomalies (bulk modulus and Poisson ratio) and
  nickel's demagnetized Delta-E dip are deliberately not represented; the
  class headers document these omissions.
- Base-class model names such as Bloch-Grueneisen, Matthiessen, Kohler,
  Pippard, and Wiedemann-Franz are used throughout, but none has a
  bibliographic entry in the tree; the quasi-harmonic elastic closure cites
  Garai-Laugier 2007.

## Citation defects in the sources

The typos listed here on 2026-08-28 (DOIs missing their leading `1` in the iron, nickel and silver
sources, "imon" for Simon in the magnesia source, "Arays" for Arajs in the iron source, "referene"
in the lead and tin sources, and "O'Neil 1964" for O'Neal 1965 in the tin source) were corrected in
the code comments on 2026-09-16, as was the Messe et al. 2023 title in `powerlaws.hpp` (now the
published one, with DOI 10.1088/1361-6668/acf7f9). Silver's file header names Neighbours and Alers
1958 for the elastic data since the same date. Nothing from that list is open.
