# BELFEM Literature References {#doc_literature_references}

**Purpose:** Reference guide to the literature used in BELFEM development.

**Note:** The `./literature/` directory contains proprietary content and is maintained as an independent repository (in `.gitignore`). This document provides references so users can locate the cited works independently.

**Last Updated:** 2026-08-11

---

## Paper Alias Index (retired — decoder for older records)

**The `paperN` aliases were retired on 2026-08-11.** New citations use author and year
("Messe et al. 2023, §2.7"); see `CLAUDE.md`. They were converted out of the source comments and
the module documentation in the same pass, so the aliases now survive only in dated devlog entries,
which are kept as written.

This table remains the authoritative alias → source file → citation map, so those older records stay
readable and auditable. Source `.txt` files live in the proprietary `literature/` repository; the
citations below let anyone locate the works independently.

| Alias | Source file | Citation |
|---|---|---|
| paper0 | `dular2021.txt` | Dular et al. 2021, "Stability of Mixed FE Formulations for HTS", IEEE TASC, DOI:10.1109/TASC.2021.3098724 |
| paper1 | `messe2023.txt` | Messe et al. 2023, "BELFEM: a special purpose FE code for magnetodynamic modeling of HTS tapes", SUST, DOI:10.1088/1361-6668/acf7f9 |
| paper2 | `messe2022.txt` | Messe 2022, "A Special Purpose Finite-Element Framework for HTS Applications", HTS 2022 Conference, HAL:hal-03791404 |
| paper3 | `arsenault2023.txt` | Arsenault et al. 2023, "Magnetodynamic H-φ Formulation", IEEE TASC, DOI:10.1109/TASC.2023.3293449 |
| paper4 | `riva2023.txt` | Riva et al. 2023, "H-φ in Sparselizard with DDM", IEEE TASC, DOI:10.1109/TASC.2023.3240389 |
| paper5 | `alves2022a.txt` | Alves et al. 2022a, "3D FE Thin-Shell Model for HTS Tapes", IEEE TASC, DOI:10.1109/TASC.2022.3143076 |
| paper6 | `alves2022b.txt` | Alves et al. 2022b, "Thin-shell approach for modeling superconducting tapes in H-φ formulation", SUST, DOI:10.1088/1361-6668/ac3f9e |
| paper7 | `arsenault2021.txt` | Arsenault et al. 2021, "Implementation of H-φ Formulation in COMSOL Multiphysics", IEEE TASC, DOI:10.1109/TASC.2020.3033998 |
| paper8 | `alves2024.txt` | Alves et al. 2024, "2-D Thin-Shell Model Based on H-φ-Formulation in COMSOL", IEEE TASC, DOI:10.1109/TASC.2024.3473850 |
| paper9 | `schnaubelt2023.txt` | Schnaubelt et al. 2023, "Electromagnetic Simulation of No-Insulation Coils Using H-φ TSA", IEEE TASC, DOI:10.1109/TASC.2023.3258905 |
| paperA | `schnaubelt2023.txt` | Same work as paper9 — July 2026 devlogs use `paperA` for it (alias drift); both resolve here |

Newer papers cited by author-year (no `paperN` alias assigned):

| Source file | Citation |
|---|---|
| `arsenault2026.txt` | Arsenault et al. 2026, Erratum to "Magnetodynamic H-φ Formulation" (corrects Arsenault et al. 2023's air-domain coupling), IEEE TASC, DOI:10.1109/TASC.2026.3686487 |
| `denis2026.txt` | Denis et al. 2026, "Simultaneous Multi-Scale Homogeneous H-Phi Thin-Shell Model for Stacked HTS Coils", IEEE TASC, DOI:10.1109/TASC.2026.3652981 |
| `dular1997.txt` | Dular et al. 1997, "A Generalized Source Magnetic Field Calculation Method for Inductors of Any Shape", IEEE Trans. Magn., DOI:10.1109/20.582518 |
| `dular1999.txt` | Dular et al. 1999, "A Natural Method for Coupling Magnetodynamic H-Formulations and Circuit Equations", IEEE Trans. Magn., DOI:10.1109/20.767308 |
| `luccini2025.txt` | Lucchini 2025, "Evaluating Magnetization Losses in 3-D CORC Tapes With Integral and Finite-Element Methods", IEEE TASC, DOI:10.1109/TASC.2025.3544512 (filename `luccini` — cite as Lucchini) |
| `wozniak2025.txt` | Wozniak et al. 2025, "Influence of Critical Current Defect on Operation, Quench Detection and Protection of a Conduction-Cooled Pancake REBCO Coil", IEEE TASC 35(5):4604006, DOI:10.1109/TASC.2025.3532246 |
| `badel2021.txt` | Badel et al. 2019, "Modeling of 'quench' or the occurrence and propagation of dissipative zones in REBCO high temperature superconducting coils", SUST 32(9):094001, DOI:10.1088/1361-6668/ab181f (filename `badel2021` is the HAL deposit year, HAL:hal-02509494 — cite as Badel et al. 2019) |

FVM papers (aliases F0–F4):

| Alias | Source file | Citation |
|---|---|---|
| F0 | `aavatsmark2002.txt` | Aavatsmark 2002, "An Introduction to Multipoint Flux Approximations for Quadrilateral Grids", Comput. Geosci. 6:405-432, DOI:10.1023/A:1021291114475 |
| F1 | `agelas2008.txt` | Agélas, Di Pietro, Masson 2008, "A Symmetric and Coercive Finite Volume Scheme for Multiphase Porous Media Flow", FVCA V, DOI:10.1515/JNUM.2008.006 |
| F2 | `klausen2006.txt` | Klausen & Winther 2006, "Robust Convergence of Multi Point Flux Approximation on Rough Grids", Numer. Math. 104:317-337, DOI:10.1007/s00211-006-0023-4 |
| F3 | `ingram2010.txt` | Ingram, Wheeler, Yotov 2010, "A Multipoint Flux Mixed Finite Element Method on Hexahedra", SIAM J. Numer. Anal. 48(4):1281-1312, DOI:10.1137/090766176 |
| F4 | `wheeler2011.txt` | Wheeler, Xue, Yotov 2011, "A Family of Multipoint Flux Mixed Finite Element Methods for Elliptic Problems on General Grids", Procedia CS 4:918-927, DOI:10.1016/j.procs.2011.04.097 |

## Citation Key Index (machine-checked)

**This table is the source of record for full references, DOIs, and ISBNs.**
`scripts/check_citations.py` reads it between the markers below, scans `src/` and `doc/` for
author-year citations, and reports any that do not resolve. Keep the table shape: the key cell may
carry several accepted spellings separated by `;`, and the third cell is the full reference.

Property-data sources are **not listed here**. They live in
`src/physics/materials/doc/material_property_sources.md`, whose bibliography the checker reads as a
second authoritative table.

<!-- BEGIN CITATION INDEX -->

| Key | Cite as | Reference | DOI / ISBN |
|---|---|---|---|
| Messe 2023 | Messe et al. 2023, §2.7 | C. Messe, N. Riva, S. Viarengo, G. Giard, F. Sirois, "BELFEM: a special purpose finite element code for the magnetodynamic modeling of high-temperature superconducting tapes", Supercond. Sci. Technol. 36:114001, 2023 | 10.1088/1361-6668/acf7f9 |
| Messe 2022 | Messe 2022 | C. Messe, "A Special Purpose Finite-Element Framework for HTS Applications", 8th Int. Workshop on Numerical Modelling of HTS, 2022 | HAL:hal-03791404 |
| Alves 2022a | Alves et al. 2022a, §III | B. de Sousa Alves, M. Laforest, F. Sirois, "A 3-D Finite-Element Thin-Shell Model for HTS Tapes", IEEE Trans. Appl. Supercond. 32:7500411, 2022 | 10.1109/TASC.2022.3143076 |
| Alves 2022b; Sousa 2022b; de Sousa Alves 2022b | Alves et al. 2022b, Eq. 19 | B. de Sousa Alves, V. Lahtinen, M. Laforest, F. Sirois, "Thin-shell approach for modeling superconducting tapes in the H-φ formulation", Supercond. Sci. Technol. 35:024001, 2022 | 10.1088/1361-6668/ac3f9e |
| Alves 2024 | Alves et al. 2024 | B. de Sousa Alves, A. Arsenault, F. Sirois, "2-D Thin-Shell Model Based on the H-φ-Formulation in COMSOL Multiphysics", IEEE Trans. Appl. Supercond. 34:7500910, 2024 | 10.1109/TASC.2024.3473850 |
| Arsenault 2021 | Arsenault et al. 2021 | A. Arsenault, F. Sirois, F. Grilli, "Implementation of the H-φ Formulation in COMSOL Multiphysics for Simulating the Magnetization of Bulk Superconductors", IEEE Trans. Appl. Supercond. 31:6800111, 2021 | 10.1109/TASC.2020.3033998 |
| Arsenault 2023 | Arsenault et al. 2023, §II | A. Arsenault, B. de Sousa Alves, G. Giard, F. Sirois, "Magnetodynamic H-φ Formulation", IEEE Trans. Appl. Supercond. 33:6801306, 2023 — **read with the 2026 erratum** | 10.1109/TASC.2023.3293449 |
| Arsenault 2026 | Arsenault et al. 2026 (erratum) | A. Arsenault, B. de Sousa Alves, G. Giard, F. Sirois, Erratum to "Magnetodynamic H-φ Formulation", IEEE Trans. Appl. Supercond. 36, 2026 — corrects the air-domain coupling | 10.1109/TASC.2026.3686487 |
| Denis 2026 | Denis et al. 2026 | L. Denis, B. Vanderheyden, C. Geuzaine, "Simultaneous Multi-Scale Homogeneous H-Phi Thin-Shell Model for Stacked HTS Coils", IEEE Trans. Appl. Supercond. 36, 2026 | 10.1109/TASC.2026.3652981 |
| Dular 1997 | Dular et al. 1997 | P. Dular, F. Henrotte, F. Robert, A. Genon, W. Legros, "A Generalized Source Magnetic Field Calculation Method for Inductors of Any Shape", IEEE Trans. Magn. 33:1398-1401, 1997 | 10.1109/20.582518 |
| Dular 1999 | Dular et al. 1999 | P. Dular, C. Geuzaine, W. Legros, "A Natural Method for Coupling Magnetodynamic H-Formulations and Circuit Equations", IEEE Trans. Magn. 35:1626-1629, 1999 | 10.1109/20.767308 |
| Dular 2021 | Dular et al. 2021 | J. Dular, M. Harutyunyan, L. Bortot, S. Schöps, B. Vanderheyden, C. Geuzaine, "On the Stability of Mixed Finite-Element Formulations for High-Temperature Superconductors", IEEE Trans. Appl. Supercond. 31:8200412, 2021 | 10.1109/TASC.2021.3098724 |
| Riva 2023 | Riva et al. 2023 | N. Riva, A. Halbach, M. Lyly, C. Messe, J. Ruuskanen, V. Lahtinen, "H-φ Formulation in Sparselizard Combined With Domain Decomposition Methods", IEEE Trans. Appl. Supercond. 33:4900405, 2023 | 10.1109/TASC.2023.3240389 |
| Riva 2021 | Riva 2021, Ch. 5 | N. Riva, "Quench behavior of high-temperature superconductor tapes for power applications", PhD thesis, EPFL, 2021 | EPFL thesis 8754 |
| Schnaubelt 2023 | Schnaubelt et al. 2023 | E. Schnaubelt, M. Wozniak, S. Schöps, A. Verweij et al., "Electromagnetic Simulation of No-Insulation Coils Using H-φ Thin Shell Approximation", IEEE Trans. Appl. Supercond. 33:4900906, 2023 | 10.1109/TASC.2023.3258905 |
| Wozniak 2025 | Wozniak et al. 2025 | M. Wozniak, E. Schnaubelt, S. Atalay, B. Bordini, J. Dular et al., "Influence of Critical Current Defect on Operation, Quench Detection and Protection of a Conduction-Cooled Pancake REBCO Coil", IEEE Trans. Appl. Supercond. 35:4604006, 2025 | 10.1109/TASC.2025.3532246 |
| Badel 2019 | Badel et al. 2019, §4.3 | A. Badel, B. Rozier, B. Ramdane, G. Meunier, P. Tixador, "Modelling of 'quench' or the occurrence and propagation of dissipative zones in REBCO high temperature superconducting coils", Supercond. Sci. Technol. 32:094001, 2019 — local file `badel2021` is the HAL deposit year | 10.1088/1361-6668/ab181f |
| Lucchini 2025 | Lucchini 2025 | F. Lucchini, "Evaluating Magnetization Losses in 3-D CORC Tapes With Integral and Finite-Element Methods", IEEE Trans. Appl. Supercond. 35, 2025 — local file is misspelled `luccini2025` | 10.1109/TASC.2025.3544512 |
| Aavatsmark 2002 | Aavatsmark 2002, §3 | I. Aavatsmark, "An Introduction to Multipoint Flux Approximations for Quadrilateral Grids", Comput. Geosci. 6:405-432, 2002 | 10.1023/A:1021291114475 |
| Agelas & Masson 2008; Agélas & Masson 2008; Agélas 2008; Agelas 2008; Masson 2008 | Agelas & Masson 2008 | L. Agélas, R. Masson, "Convergence of the finite volume MPFA O scheme for heterogeneous anisotropic diffusion problems on general meshes", C. R. Acad. Sci. Paris 346:1007-1012, 2008 — **see the note below on the FVCA V paper** | 10.1016/j.crma.2008.07.015 |
| Ingram 2010 | Ingram et al. 2010 | R. Ingram, M. F. Wheeler, I. Yotov, "A Multipoint Flux Mixed Finite Element Method on Hexahedra", SIAM J. Numer. Anal. 48:1281-1312, 2010 | 10.1137/090766176 |
| Klausen & Winther 2006; Klausen 2006 | Klausen & Winther 2006 | R. A. Klausen, R. Winther, "Robust Convergence of Multi Point Flux Approximation on Rough Grids", Numer. Math. 104:317-337, 2006 | 10.1007/s00211-006-0023-4 |
| Wheeler 2011 | Wheeler et al. 2011 | M. F. Wheeler, G. Xue, I. Yotov, "A Family of Multipoint Flux Mixed Finite Element Methods for Elliptic Problems on General Grids", Procedia Comput. Sci. 4:918-927, 2011 | 10.1016/j.procs.2011.04.097 |
| Pellikka 2013 | Pellikka et al. 2013, Eq. 10 | M. Pellikka, S. Suuriniemi, L. Kettunen, C. Geuzaine, "Homology and Cohomology Computation in Finite Element Modeling", SIAM J. Sci. Comput. 35:B1195-B1214, 2013 | 10.1137/130906556 |
| Chen & Freedman 2010; Chen 2010; Chen & Freedman 2011 | Chen & Freedman 2010 | C. Chen, D. Freedman, "Hardness Results for Homology Localization", Discrete Comput. Geom. 45:425-448 — conference version 2010, journal version 2011 | 10.1007/s00454-010-9322-8 |
| Dey 2011 | Dey et al. 2011 | T. K. Dey, A. N. Hirani, B. Krishnamoorthy, "Optimal Homologous Cycles, Total Unimodularity, and Linear Programming", SIAM J. Comput. 40:1026-1044, 2011 | 10.1137/100800245 |
| Dunfield & Hirani 2011; Dunfield 2011 | Dunfield & Hirani 2011 | N. M. Dunfield, A. N. Hirani, "The least spanning area of a knot and the optimal bounding chain problem", Proc. 27th Ann. Symp. Computational Geometry, 135-144, 2011 | 10.1145/1998196.1998218 |
| Mrozek & Batko 2009; Mrozek 2009 | Mrozek & Batko 2009 | M. Mrozek, B. Batko, "Coreduction Homology Algorithm", Discrete Comput. Geom. 41:96-118, 2009 | 10.1007/s00454-008-9073-y |
| Haken 1961 | Haken 1961 | W. Haken, "Theorie der Normalflächen", Acta Math. 105:245-375, 1961 | 10.1007/BF02559591 |
| Costantini 1998 | Costantini 1998 | M. Costantini, "A Novel Phase Unwrapping Method Based on Network Programming", IEEE Trans. Geosci. Remote Sens. 36:813-821, 1998 | 10.1109/36.673674 |
| Bathe 2016; Bathe 2014 | Bathe, §8.6 | K.-J. Bathe, *Finite Element Procedures*, 2nd ed., published by the author, Watertown MA; copyright 2014, the copy here is the 2016 4th printing. Prentice Hall published the **first** edition only | ISBN 978-0-9790049-5-7 |
| Hughes 2000; Hughes 1987 | Hughes, Ch. 4 | T. J. R. Hughes, *The Finite Element Method: Linear Static and Dynamic Finite Element Analysis*. Dover, 2000 — unabridged reprint of Prentice-Hall, 1987; the local extraction is the 1987 printing | ISBN 0-13-317025-X (1987) |
| Zienkiewicz 2013 | Zienkiewicz et al. 2013 (Vol. 1) | O. C. Zienkiewicz, R. L. Taylor, J. Z. Zhu, *The Finite Element Method: Its Basis and Fundamentals*, 7th ed. Butterworth-Heinemann (Elsevier), 2013 | ISBN 978-1-85617-633-0 |
| Zienkiewicz 2014 | Zienkiewicz et al. 2014 (Vol. 2) | O. C. Zienkiewicz, R. L. Taylor, D. D. Fox, *The Finite Element Method for Solid and Structural Mechanics*, 7th ed. Butterworth-Heinemann (Elsevier), 2014 | ISBN 978-1-85617-634-7 |
| Belytschko 2014 | Belytschko et al. 2014 | T. Belytschko, W. K. Liu, B. Moran, K. I. Elkhodary, *Nonlinear Finite Elements for Continua and Structures*, 2nd ed. Wiley, 2014 | ISBN 978-1-118-63270-3 |
| Brenner & Scott 2008; Brenner 2008 | Brenner & Scott 2008, Ch. 5 | S. C. Brenner, L. R. Scott, *The Mathematical Theory of Finite Element Methods*, 3rd ed. Springer, 2008 | ISBN 978-0-387-75933-3 |
| Boffi 2013 | Boffi et al. 2013, Ch. 3 | D. Boffi, F. Brezzi, M. Fortin, *Mixed Finite Element Methods and Applications*, Springer Series in Computational Mathematics 44. Springer, 2013 | 10.1007/978-3-642-36519-5 |
| Monk 2003 | Monk 2003, §6.1 | P. Monk, *Finite Element Methods for Maxwell's Equations*. Oxford University Press, 2003 | ISBN not captured in the local extraction |
| Arnold 2018 | Arnold 2018 | D. N. Arnold, *Finite Element Exterior Calculus*, CBMS-NSF Regional Conference Series in Applied Mathematics 93. SIAM, 2018 — **single author**; Falk and Winther are co-authors of the Acta Numerica articles, not of this book | ISBN 978-1-61197-553-6 |
| Bronshtein 2015 | Bronshtein | I. N. Bronshtein, K. A. Semendyayev, G. Musiol, H. Mühlig, *Handbook of Mathematics*, 6th ed. Springer, 2015 | ISBN 978-3-662-46220-1 |
| Cormen 2009 | Cormen et al. (CLRS), Ch. 26 | T. H. Cormen, C. E. Leiserson, R. L. Rivest, C. Stein, *Introduction to Algorithms*, 3rd ed. MIT Press, 2009 | ISBN 978-0-262-03384-8 |
| Gross & Kotiuga 2004; Gross 2004 | Gross & Kotiuga, §3 | P. W. Gross, P. R. Kotiuga, *Electromagnetic Theory and Computation: A Topological Approach*, MSRI Publications 48. Cambridge University Press, 2004 | ISBN 0-521-80160-5 |
| Iwasa 2009 | Iwasa 2009, §8.2 | Y. Iwasa, *Case Studies in Superconducting Magnets: Design and Operational Issues*, 2nd ed. Springer, 2009 | ISBN 978-0-387-09799-2 |
| Russenschuck 2010 | Russenschuck 2010, Ch. 7 | S. Russenschuck, *Field Computation for Accelerator Magnets: Analytical and Numerical Methods for Electromagnetic Design and Optimization*. Wiley-VCH, 2010 | ISBN 978-3-527-40769-9 |
| Glukhikh 2018 | Glukhikh et al. 2018 | V. Glukhikh, O. Filatov, B. Kolbasov (eds.), *Fundamentals of Magnetic Thermonuclear Reactor Design*. Woodhead Publishing (Elsevier), 2018 | ISBN 978-0-08-102470-6 |
| Barth 2013 | Barth 2013, Ch. 3 | C. Barth, *High Temperature Superconductor Cable Concepts for Fusion Magnets*, dissertation, Karlsruhe Institute of Technology. KIT Scientific Publishing, 2013 | ISBN 978-3-7315-0065-0 |
| Martin 2025 | Martin 2025, Ch. 5 | R. C. Martin, *Clean Code: A Handbook of Agile Software Craftsmanship*, 2nd ed. Pearson, 2025 — the copyright page reads 2026; BELFEM documents cite it as 2025 | ISBN 978-0-13-539857-9 |
| Oliveira & Stewart 2006; Oliveira 2006 | Oliveira & Stewart 2006, §7.9 | S. Oliveira, D. E. Stewart, *Writing Scientific Software: A Guide to Good Style*. Cambridge University Press, 2006 | 10.1017/CBO9780511617973 |
| Rouson 2011 | Rouson et al. 2011, §1.7 | D. Rouson, J. Xia, X. Xu, *Scientific Software Design: The Object-Oriented Way*. Cambridge University Press, 2011 | 10.1017/CBO9780511977381 |
| Thomas & Hunt 2020; Thomas 2020; Hunt 2020 | Thomas & Hunt 2020, Topic 25 | D. Thomas, A. Hunt, *The Pragmatic Programmer: Your Journey to Mastery*, 20th anniversary ed. Addison-Wesley, 2020 — cite by Topic and Tip, never by page | ISBN 978-0-13-595705-9 |
| Evans 2017 | Evans 2017, §3.2 | J. A. Evans, *Mathematical Foundations of the Finite Element Method*, lecture series, University of Colorado Boulder, 2017 | lecture notes |
| Felippa IFEM | Felippa IFEM, Ch. 11 | C. A. Felippa, *Introduction to Finite Element Methods (IFEM)*, course notes, University of Colorado Boulder | course notes |
| Rhyner 1993 | Rhyner 1993 | J. Rhyner, "Magnetic properties and AC-losses of superconductors with power-law current-voltage characteristics", Physica C 212:292-300, 1993 | DOI not confirmed locally |
| Duron 2004 | Duron et al. 2004 | J. Duron, F. Grilli, B. Dutoit, S. Stavrev, "Modeling the E-J relation of high-Tc superconductors in an arbitrary current range", Physica C 401:231-235, 2004 | DOI not confirmed locally |
| Plummer 1987 | Plummer & Evetts 1987 | C. P. Plummer, J. E. Evetts, "Dependence of the shape of the resistive transition on composite inhomogeneity in multifilamentary wires", IEEE Trans. Magn. 23(2):1179-1182, 1987 | DOI not confirmed locally |
| Kim 1962 | Kim 1962 | Y. B. Kim, C. F. Hempstead, A. R. Strnad, "Critical persistent currents in hard superconductors", Phys. Rev. Lett. 9, 1962 — the Kim model, reached here through Messe et al. 2023 §2.6 | DOI not confirmed locally |
| Arnold 2001 | Arnold et al. 2001 | D. N. Arnold, D. Boffi, R. S. Falk, L. Gastaldi, "Finite element approximation on quadrilateral meshes", Commun. Numer. Methods Eng. 17(11):805-812, 2001 | 10.1002/cnm.450 |
| Arnold 2002 | Arnold et al. 2002 | D. N. Arnold, D. Boffi, R. S. Falk, "Approximation by quadrilateral finite elements", Math. Comp. 71(239):909-922, 2002 | 10.1090/S0025-5718-02-01439-4 |
| Arnold 2005 | Arnold et al. 2005 | D. N. Arnold, D. Boffi, R. S. Falk, "Quadrilateral H(div) finite elements", SIAM J. Numer. Anal. 42(6):2429-2451, 2005 | 10.1137/S0036142903431924 |
| Falk 2011 | Falk et al. 2011 | R. S. Falk, P. Gatto, P. Monk, "Hexahedral H(div) and H(curl) finite elements", ESAIM: M2AN 45(1):115-143, 2011 | 10.1051/m2an/2010034 |
| Lüthi 1960; Luethi 1960 | Luethi 1960 | B. Luethi, *Widerstandsänderung von Metallen in hohen Magnetfeldern*, dissertation, ETH Zürich, 1960 — spelled `Luethi` in the materials bibliography | thesis |

| Anderson 1965 | Anderson 1965 | D. G. Anderson, "Iterative procedures for nonlinear integral equations", J. ACM 12, 1965 | 10.1145/321296.321305 |
| Walker & Ni 2011; Walker 2011; Ni 2011 | Walker & Ni 2011 | H. F. Walker, P. Ni, "Anderson acceleration for fixed-point iterations", SIAM J. Numer. Anal. 49, 2011 — the type-II form implemented here | 10.1137/10078356X |
| Irons & Tuck 1969; Irons 1969 | Irons & Tuck 1969 | B. Irons, R. C. Tuck, "A version of the Aitken accelerator for computer iteration", Int. J. Numer. Meth. Eng. 1, 1969 | DOI not confirmed locally |
| Küttler & Wall 2008; Kuettler & Wall 2008; Küttler 2008 | Küttler & Wall 2008 | U. Küttler, W. A. Wall, dynamic-relaxation fixed-point FSI coupling, Comput. Mech., 2008 | 10.1007/s00466-008-0255-5 |
| Gustafsson 1988; Lundh & Söderlind 1988; Lundh 1988; Söderlind 1988 | Gustafsson, Lundh & Söderlind 1988 | K. Gustafsson, M. Lundh, G. Söderlind, "A PI stepsize control for the numerical solution of ordinary differential equations", BIT 28:270-287, 1988 | 10.1007/BF01934091 |
| Gustafsson 1994 | Gustafsson 1994 | K. Gustafsson, "Control-theoretic techniques for stepsize selection in implicit Runge-Kutta methods", ACM TOMS 20(4):496-517, 1994 | 10.1145/198429.198437 |
| Söderlind 2002; Soderlind 2002 | Söderlind 2002 | G. Söderlind, "Automatic control and adaptive time-stepping", Numerical Algorithms 31:281-310, 2002 | 10.1023/A:1021160023092 |
| Söderlind 2003; Soderlind 2003 | Söderlind 2003 | G. Söderlind, "Digital filters in adaptive time-stepping", ACM TOMS 29(1):1-26, 2003 | 10.1145/641876.641877 |
| Söderlind & Wang 2006; Soderlind & Wang 2006; Söderlind 2006 | Söderlind & Wang 2006 | G. Söderlind, L. Wang, "Adaptive time-stepping and computational stability", J. Comput. Appl. Math. 185:225-243, 2006 | 10.1016/j.cam.2005.03.008 |
| Valli 2002; Carey & Coutinho 2002; Carey 2002; Coutinho 2002 | Valli, Carey & Coutinho 2002 | A. M. P. Valli, G. F. Carey, A. L. G. A. Coutinho, "Control strategies for timestep selection in simulation of coupled viscous flow and heat transfer", Commun. Numer. Meth. Engng 18(2):131-139, 2002 | 10.1002/cnm.475 |
| Brenan 1996; Campbell & Petzold 1996; Campbell 1996; Petzold 1996 | Brenan, Campbell & Petzold 1996 | K. E. Brenan, S. L. Campbell, L. R. Petzold, *Numerical Solution of Initial-Value Problems in Differential-Algebraic Equations*. SIAM, 1996 | 10.1137/1.9781611971224 |
| Hairer & Wanner 1996; Hairer 1996; Wanner 1996 | Hairer & Wanner 1996, II.4 | E. Hairer, G. Wanner, *Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems*. Springer, 1996 | DOI not confirmed locally |
| Burman & Zunino 2006; Burman 2006; Zunino 2006 | Burman & Zunino 2006 | E. Burman, P. Zunino, weighted interior penalty for advection-diffusion with discontinuous coefficients, SIAM J. Numer. Anal. 44:1612-1638, 2006 | DOI not confirmed locally |
| Ern 2009; Stephansen & Zunino 2009; Stephansen 2009; Zunino 2009 | Ern, Stephansen & Zunino 2009 | A. Ern, A. F. Stephansen, P. Zunino, symmetric weighted interior penalty with harmonic flux weighting, IMA J. Numer. Anal. 29:235-256, 2009 | DOI not confirmed locally |
| Hansbo & Hansbo 2004; Hansbo 2004 | Hansbo & Hansbo 2004 | A. Hansbo, P. Hansbo, element-doubling / XFEM-equivalent formulation, Comput. Methods Appl. Mech. Engrg. 193(33-35):3523-3540, 2004 | 10.1016/j.cma.2003.12.041 |
| Karypis & Kumar 1998; Karypis 1998; Kumar 1998 | Karypis & Kumar 1998 | G. Karypis, V. Kumar, multilevel nested dissection, the algorithm behind METIS, 1998 | DOI not confirmed locally |
| George 1973 | George 1973 | A. George, "Nested Dissection of a Regular Finite Element Mesh", SIAM J. Numer. Anal., 1973 | DOI not confirmed locally |
| Hirani & Krishnamoorthy 2011; Krishnamoorthy 2011; Hirani 2011 | Dey et al. 2011 | see **Dey 2011** — the same work, cited in `src/homology/doc/` by its second and third authors | 10.1137/100800245 |
| Freedman 2010; Freedman 2011 | Chen & Freedman 2010 | see **Chen & Freedman 2010** — the same work cited by its second author | 10.1007/s00454-010-9322-8 |
| Geuzaine & Legros 1999; Geuzaine 1999; Legros 1999 | Dular et al. 1999 | see **Dular 1999** — cited in `src/homology/doc/` by its second and third authors | 10.1109/20.767308 |
| Kaczynski 2004; Kaczynski & Mischaikow & Mrozek 2004 | Kaczynski et al. 2004 | T. Kaczynski, K. Mischaikow, M. Mrozek, *Computational Homology*, Applied Mathematical Sciences. Springer, 2004 — the chain complex reduction (CCR) and Smith normal form behind `fn_Smith.cpp` | DOI not confirmed locally |
| O'Neal 1965; O'Neil 1965; O'Neal 1964; O'Neil 1964 | O'Neal et al. 1965 | H. R. O'Neal et al., Phys. Rev. 137, A748, 1965 — the materials bibliography keys this as `O'Neil 1965`, and the code cites it as 1964; the author is O'Neal | 10.1103/PhysRev.137.A748 |
| Launey 1959; Launay 1959; De Launay 1959 | De Launay 1959 | J. De Launay et al., J. Phys. Chem. Solids, 1959 — spelled `Launey` in `cl_Material_Copper.cpp` | 10.1016/0022-3697(59)90038-1 |

| Amestoy 2000; Amestoy & Duff 2000; Duff & L'Excellent 2000; Duff 2000; L'Excellent 2000 | Amestoy, Duff & L'Excellent 2000 | P. R. Amestoy, I. S. Duff, J.-Y. L'Excellent, "Multifrontal parallel distributed symmetric and unsymmetric solvers", Comput. Methods Appl. Mech. Engrg. 184:501-520, 2000 — the MUMPS paper, reference [24] of Messe et al. 2023 | 10.1016/S0045-7825(99)00242-X |
| Ghysels 2016; Ghysels & Li 2016 | Ghysels et al. 2016 | P. Ghysels, X. S. Li, F.-H. Rouet, S. Williams, A. Napov, "An Efficient Multicore Implementation of a Novel HSS-Structured Multifrontal Solver Using Randomized Sampling", SIAM J. Sci. Comput. 38(5):S358-S384, 2016 — the STRUMPACK paper, reference [26] of Messe et al. 2023 | 10.1137/15M1010117 |
| Grigorieff 1983 | Grigorieff 1983 | R. D. Grigorieff, "Stability of multistep-methods on variable grids", Numer. Math. 42:359-377, 1983 — the variable-grid zero-stability bound the BDF growth clamp rests on | 10.1007/BF01389580 |

| Nédélec 1980; Nedelec 1980 | Nédélec 1980 | J. C. Nédélec, "Mixed finite elements in R³", Numer. Math. 35:315-341, 1980 — the first-kind edge-element family | 10.1007/BF01396415 |
| Nédélec 1986; Nedelec 1986 | Nédélec 1986 | J. C. Nédélec, "A new family of mixed finite elements in R³", Numer. Math. 50:57-81, 1986 — the second-kind family | 10.1007/BF01389668 |
| Zaglmayr 2006 | Zaglmayr 2006 | S. Zaglmayr, *High Order Finite Element Methods for Electromagnetic Field Computation*, PhD thesis, Johannes Kepler University Linz, 2006 | thesis |
| Schoof & Yarberry 1994; Schoof 1994; Yarberry 1994 | Schoof & Yarberry 1994 | L. A. Schoof, V. R. Yarberry, "EXODUS II: A Finite Element Data Model", SAND92-2137, Sandia National Laboratories, 1994 — the node and edge numbering BELFEM follows | SAND92-2137 |
| Bíró & Preis 1989; Biro & Preis 1989; Bíró 1989; Preis 1989 | Bíró & Preis 1989 | O. Bíró, K. Preis, "On the use of the magnetic vector potential in the finite-element analysis of three-dimensional eddy currents", IEEE Trans. Magn. 25(4):3145-3159, 1989 — the classical Coulomb-gauge penalty for the A-formulation | 10.1109/20.34388 |
| Juntunen & Stenberg 2009; Juntunen & Stenberg 2008; Juntunen 2009; Stenberg 2009 | Juntunen & Stenberg 2009 | M. Juntunen, R. Stenberg, "Nitsche's method for general boundary conditions", Math. Comp. 78(267):1353-1374 — printed in the 2009 volume; Crossref records the 2008 online date | 10.1090/S0025-5718-08-02183-2 |
| Shahbazi 2005 | Shahbazi 2005 | K. Shahbazi, "An explicit expression for the penalty parameter of the interior penalty method", J. Comput. Phys. 205:401-407, 2005 | 10.1016/j.jcp.2004.11.017 |
| Hayasaka & Ito 2019; Hayasaka 2019; Ito 2019 | Hayasaka & Ito 2019 | R. Hayasaka, S. Ito, H. Hashizume, "Evaluation of Interface Resistance in a REBCO Tape at Different Temperatures by Contact-Probing Current Transfer Length Method", IEEE Trans. Appl. Supercond. 29:1-5, 2019 | 10.1109/TASC.2019.2895603 |

| Leachman 2009; Leachman & Jacobsen 2009; Penoncello & Lemmon 2009; Penoncello 2009; Lemmon 2009 | Leachman et al. 2009 | J. W. Leachman, R. T Jacobsen, S. G. Penoncello, E. W. Lemmon, fundamental Helmholtz equation of state for parahydrogen, normal hydrogen and orthohydrogen, J. Phys. Chem. Ref. Data, 2009 — validity to 1000 K and 2000 MPa | 10.1063/1.3160306 |
| Muzny 2013; Muzny & Huber 2013; Huber 2013; Kazakov 2013 | Muzny et al. 2013 | C. D. Muzny, M. L. Huber, A. F. Kazakov, viscosity correlation for normal hydrogen, J. Chem. Eng. Data 58:969-979, 2013. **BELFEM applies the 2022 erratum**, J. Chem. Eng. Data 67:2855, 10.1021/acs.jced.2c00523 — without it the correlation is wrong | 10.1021/je301273j |
| Assael 2011; Assael & Huber 2011 | Assael et al. 2011 | M. J. Assael, S. K. Assael, M. L. Huber, R. A. Perkins, Y. Takata, thermal conductivity of hydrogen, J. Phys. Chem. Ref. Data 40:033101, 2011 — separate coefficient tables for normal and parahydrogen | 10.1063/1.3606499 |
| Lemmon & Jacobsen 2004; Lemmon 2004; Jacobsen 2004 | Lemmon & Jacobsen 2004 | E. W. Lemmon, R. T Jacobsen, viscosity and thermal conductivity of nitrogen, argon, oxygen and air, Int. J. Thermophys. 25:21-69, 2004 — one functional form, four fluids | 10.1023/B:IJOT.0000022327.04529.f3 |
| Friend 1989; Friend & Ely 1989; Ely 1989; Ingham 1989 | Friend et al. 1989 | D. G. Friend, J. F. Ely, H. Ingham, "Thermophysical Properties of Methane", J. Phys. Chem. Ref. Data 18(2):583-638, 1989 | 10.1063/1.555828 |

| Soave 1972 | Soave 1972 | G. Soave, "Equilibrium Constants from a Modified Redlich-Kwong Equation of State", Chem. Eng. Sci. 27(6):1197-1203, 1972 | 10.1016/0009-2509(72)80096-4 |
| Peng & Robinson 1976; Peng 1976 | Peng & Robinson 1976 | D.-Y. Peng, D. B. Robinson, "A New Two-Constant Equation of State", Ind. Eng. Chem. Fundam. 15(1):59-64, 1976 | 10.1021/i160057a011 |
| Setzmann & Wagner 1991; Setzmann 1991; Wagner 1991 | Setzmann & Wagner 1991 | U. Setzmann, W. Wagner, "A New Equation of State and Tables of Thermodynamic Properties for Methane", J. Phys. Chem. Ref. Data 20(6):1061-1155, 1991 | 10.1063/1.555898 |
| Schmidt & Wagner 1985; Schmidt 1985; Wagner 1985 | Schmidt & Wagner 1985 | R. Schmidt, W. Wagner, "A New Form of the Equation of State for Pure Substances", Fluid Phase Equilib. 19(3):175-200, 1985 | 10.1016/0378-3812(85)87016-3 |
| Span 2000 | Span et al. 2000 | R. Span, E. W. Lemmon, R. T Jacobsen, W. Wagner, A. Yokozeki, "A Reference Equation of State for the Thermodynamic Properties of Nitrogen", J. Phys. Chem. Ref. Data 29(6):1361-1433, 2000 | 10.1063/1.1349047 |
| Wilke 1950 | Wilke 1950 | C. R. Wilke, "A Viscosity Equation for Gas Mixtures", J. Chem. Phys. 18:517-519, 1950 | 10.1063/1.1747673 |
| Gordon & McBride 1994; Gordon 1994; McBride 1994 | Gordon & McBride 1994 | S. Gordon, B. J. McBride, "Computer Program for Calculation of Complex Chemical Equilibrium Compositions and Applications", NASA Reference Publication 1311, 1994 — the CEA program; mixture rules Eqs. (5.3)-(5.7) | NASA RP-1311 |
| Prausnitz & Poling 1987; Reid 1987; Poling 1987 | Reid, Prausnitz & Poling 1987 | R. C. Reid, J. M. Prausnitz, B. E. Poling, *The Properties of Gases and Liquids*, 4th ed. McGraw-Hill, 1987 | book |
| Prausnitz & O'Connell 2001; Poling 2001; O'Connell 2001 | Poling, Prausnitz & O'Connell 2001 | B. E. Poling, J. M. Prausnitz, J. P. O'Connell, *The Properties of Gases and Liquids*, 5th ed. McGraw-Hill, 2001 | book |
| Stewart & Lightfoot 2007; Bird 2007; Lightfoot 2007 | Bird, Stewart & Lightfoot 2007 | R. B. Bird, W. E. Stewart, E. N. Lightfoot, *Transport Phenomena*, 2nd ed. Wiley, 2007 | book |
| Smith & Missen 1982; Missen 1982 | Smith & Missen 1982 | W. R. Smith, R. W. Missen, *Chemical Reaction Equilibrium Analysis*. Wiley-Interscience, 1982 | book |
| Anderson 2003 | Anderson 2003 | J. D. Anderson, *Modern Compressible Flow*, 3rd ed. McGraw-Hill, 2003 | book |
| Shapiro 1953 | Shapiro 1953 | A. H. Shapiro, *The Dynamics and Thermodynamics of Compressible Fluid Flow*. Ronald Press, 1953 | book |
| Cuthill & McKee 1969; Cuthill 1969; McKee 1969 | Cuthill & McKee 1969 | E. Cuthill, J. McKee, "Reducing the bandwidth of sparse symmetric matrices", Proc. 1969 24th ACM National Conference, 157-172 | 10.1145/800195.805928 |
| Liu & Sherman 1976; Sherman 1976 | Liu & Sherman 1976 | W.-H. Liu, A. H. Sherman, "Comparative Analysis of the Cuthill-McKee and the Reverse Cuthill-McKee Ordering Algorithms for Sparse Matrices", SIAM J. Numer. Anal. 13:198-213, 1976 | 10.1137/0713020 |
| Stockmeyer 1976; Gibbs 1976; Poole 1976; Gibbs & Poole 1976 | Gibbs, Poole & Stockmeyer 1976 | N. E. Gibbs, W. G. Poole Jr., P. K. Stockmeyer, "An Algorithm for Reducing the Bandwidth and Profile of a Sparse Matrix", SIAM J. Numer. Anal. 13:236-250, 1976 — the RCM start-vertex heuristic | 10.1137/0713023 |
| Dunavant 1985 | Dunavant 1985 | D. A. Dunavant, "High degree efficient symmetrical Gaussian quadrature rules for the triangle", Int. J. Numer. Methods Eng. 21:1129-1148, 1985 | 10.1002/nme.1620210612 |
| Stroud 1971 | Stroud 1971 | A. H. Stroud, *Approximate Calculation of Multiple Integrals*. Prentice-Hall, 1971 | book |
| Gear 1971 | Gear 1971 | C. W. Gear, *Numerical Initial Value Problems in Ordinary Differential Equations*. Prentice-Hall, 1971 — the BDF methods | book |
| Jin 2014 | Jin 2014, Ch. 6 | J.-M. Jin, *The Finite Element Method in Electromagnetics*, 3rd ed. Wiley, 2014 | book |
| Zienkiewicz & Zhu 1992; Zienkiewicz 1992; Zhu 1992 | Zienkiewicz & Zhu 1992 | O. C. Zienkiewicz, J. Z. Zhu, "The superconvergent patch recovery and a posteriori error estimates", Int. J. Numer. Methods Eng. 33(7): Part 1 1331-1364 (10.1002/nme.1620330702), Part 2 1365-1382 — the SPR algorithm the Maxwell postprocessor implements | 10.1002/nme.1620330703 |
| Zienkiewicz & Zhu 1987; Zienkiewicz 1987; Zhu 1987 | Zienkiewicz & Zhu 1987 | O. C. Zienkiewicz, J. Z. Zhu, "A simple error estimator and adaptive procedure for practical engineering analysis", Int. J. Numer. Methods Eng. 24(2):337-357, 1987 | 10.1002/nme.1620240206 |

<!-- END CITATION INDEX -->

**The Agélas ambiguity.** The retired alias table lists `F1` as Agélas, Di Pietro & Masson 2008,
"A Symmetric and Coercive Finite Volume Scheme for Multiphase Porous Media Flow" (FVCA V,
DOI:10.1515/JNUM.2008.006). The local extraction, `agelas2008.txt`, is a different work — the
two-author C. R. Acad. Sci. note above. Both are real; BELFEM routes to the coercivity result in
the local one. Cite whichever work you actually read, and keep the DOIs apart.

---

## Overview

BELFEM development is informed by four complementary literature categories:

1. **FEM Fundamentals** — Authoritative textbooks and lecture notes on finite element theory, implementation, and numerical methods
2. **BELFEM FEM Specialization** — Research papers specific to HTS modeling, h-φ formulations, and thin-shell approximations
3. **BELFEM FVM Methods** — Research papers on finite volume methods, specifically MPFA for diffusion problems
4. **Computational Topology and Optimal Cuts** — Algebraic-topology foundations and optimization theory behind cut generation in multiply-connected domains (the `homology/` module)

---

## Quick Reference: Which Literature Do I Need?

### For BELFEM-Specific Questions

Start with **research papers** if your question relates to:
- h-φ formulation and interface coupling
- Thin-shell models for HTS tapes
- Transport current constraints and cuts (thin/thick)
- Cohomology and multiply-connected domains
- HTS material models (E-J power law, Kim model)
- Validation benchmarks for superconductor applications

### For BELFEM FVM Questions

Start with **FVM research papers** if your question relates to:
- Multipoint Flux Approximation (MPFA) methods
- Anisotropic diffusion on general grids
- Transmissibility matrix computation
- Physical vs reference space evaluation
- Enhanced BDDF₁ spaces for 3D hexahedra
- Coercivity and convergence conditions

### For Cuts, Cohomology, and Topology Questions

Start with the **topology references** if your question relates to:
- Why cuts exist and what they are topologically
- How a cohomology basis is computed on an FE mesh (BELFEM's thick cuts)
- Whether a *minimal* cut is computable, and at what cost
- Reduction/coreduction preprocessing before Smith normal form
- The graph and LP algorithms underneath any of the above

### For General FEM Theory

Start with **textbooks** if your question relates to:
- Weak formulations and Galerkin methods
- Element technology (isoparametric, mixed, reduced integration)
- Locking phenomena (volumetric, shear, membrane)
- Inf-sup stability and mixed formulations
- Solver strategies (direct, iterative, eigenvalue)
- Nonlinear solvers (Newton-Raphson, arc-length)
- Time integration and dynamic analysis

---

## FEM Textbooks

### Core References

**Bathe, K.-J.** *Finite Element Procedures*, 2nd ed., published by the author, Watertown MA, copyright 2014; the copy held here is the 2016 4th printing. ISBN 978-0-9790049-5-7. Prentice Hall published only the first edition.
- **Use for:** Robustness, reliability, production code implementation
- **Key chapters:** Ch. 4 (formulation), Ch. 5 (elements), Ch. 8 (solvers), Ch. 9 (time integration)
- **Personality:** Disciplined engineer — emphasizes proven, robust methods

**Hughes, T.J.R.** *The Finite Element Method: Linear Static and Dynamic Finite Element Analysis*. Dover, 2000 (reprint of 1987 original).
- **Use for:** Clear linear FEM implementations, algorithm design
- **Key chapters:** Ch. 1 (weak forms), Ch. 4 (locking), Ch. 5 (isoparametric elements)
- **Personality:** Algorithm designer — clean, systematic approach

**Zienkiewicz, O.C., Taylor, R.L., and Zhu, J.Z.** *The Finite Element Method: Its Basis and Fundamentals*, 7th ed. Butterworth-Heinemann, 2013.
- **Use for:** Comprehensive reference, adaptivity
- **Personality:** Encyclopedia — covers everything

**Zienkiewicz, O.C., and Taylor, R.L.** *The Finite Element Method for Solid and Structural Mechanics*, 7th ed. Butterworth-Heinemann, 2014.
- **Use for:** Nonlinear implementation, exact shell elements
- **Key chapters:** Ch. 5 (nonlinear theory), Ch. 8 (contact)

### Advanced References

**Belytschko, T., Liu, W.K., Moran, B., and Elkhodary, K.** *Nonlinear Finite Elements for Continua and Structures*, 2nd ed. Wiley, 2014.
- **Use for:** Nonlinear mechanics, XFEM, stability theory
- **Not for:** First-time learning (dense presentation)

**Brenner, S.C., and Scott, L.R.** *The Mathematical Theory of Finite Element Methods*, 3rd ed. Springer, 2008.
- **Use for:** Rigorous proofs, convergence theory
- **Key chapters:** Ch. 4-5 (convergence), Ch. 8 (mixed methods)

**Bronshtein, I.N., Semendyayev, K.A., Musiol, G., and Mühlig, H.** *Handbook of Mathematics*, 6th ed. Springer, 2015.
- **Use for:** Quick formula verification and lookup

### Lecture Notes

**Evans, J.A.** *Mathematical Foundations of the Finite Element Method*, Lecture Notes, University of Colorado Boulder, 2017.
- **Use for:** The functional-analysis layer under everything else — Sobolev spaces, distributions, Lax-Milgram, Babuška-Brezzi, interpolation error estimates
- **Cite as:** "Evans 2017, Ch. N" or "Evans 2017, §N.M"
- **Pairs with:** Brenner & Scott (proofs in book form), Boffi et al. (mixed extensions)

**Felippa, C.A.** *Introduction to Finite Element Methods (IFEM)*, Course Notes, University of Colorado Boulder.
- **Use for:** Step-by-step learning, clear derivations; the penalty-method and Lagrange-multiplier tutorials, and static condensation
- **Cite as:** "Felippa IFEM, Ch. N, §X.Y"
- **Available:** Online as course materials

### Specialized Advanced References

**Arnold, D.N.** *Finite Element Exterior Calculus*, CBMS-NSF Regional Conference Series in Applied Mathematics 93. SIAM, 2018. ISBN 978-1-61197-553-6. Arnold is the sole author; Falk and Winther co-authored the Acta Numerica articles, not this book.
- **Use for:** FEEC theory, topological foundations of FEM spaces
- **Key topics:** de Rham complex, structure preservation, why FEM spaces must exist
- **Personality:** Topologist — explains the mathematical structure underlying FEM

**Boffi, D., Brezzi, F., and Fortin, M.** *Mixed Finite Element Methods and Applications*. Springer, 2013.
- **Use for:** Mixed FEM theory, saddle-point problems, inf-sup stability
- **Key chapters:** Ch. 3-5 (inf-sup conditions), Ch. 11 (Maxwell eigenvalues)
- **Personality:** Mixed FEM theorist — comprehensive treatment of stability

**Monk, P.** *Finite Element Methods for Maxwell's Equations*. Oxford University Press, 2003.
- **Use for:** Maxwell equations, edge elements (Nédélec), electromagnetic scattering
- **Key chapters:** Ch. 5-6 (edge elements), Ch. 7 (convergence), Ch. 9-12 (scattering)
- **Personality:** Maxwell specialist — definitive source for electromagnetics FEM

### Element Distortion on Quadrilateral/Hexahedral Meshes

Why QUAD/HEX elements in the magnetic solve must be perfectly rectangular: see
`src/fem/interpolation/doc/nedelec.md` §6.6 and the maxwell usage guide §1.7. On
non-affine (bi-/trilinearly mapped) quads and hexes, the mapped H(curl)/H(div) spaces lose
completeness, and the lowest-order elements can lose convergence entirely. Monk 2003
develops the hexahedral edge-element theory only for parallelepipeds (§6.1) and cites this
line of work as the reason (§8.2–8.3). Boffi et al. 2013 summarize it in §2.2.4/§2.5.5.

**Arnold, D.N., Boffi, D., Falk, R.S., and Gastaldi, L.** "Finite element approximation on
quadrilateral meshes". *Communications in Numerical Methods in Engineering*, 17(11):805-812,
2001. DOI:10.1002/cnm.450
- Necessary and sufficient conditions for approximation order under bilinear maps

**Arnold, D.N., Boffi, D., and Falk, R.S.** "Approximation by quadrilateral finite
elements". *Mathematics of Computation*, 71(239):909-922, 2002.
DOI:10.1090/S0025-5718-02-01439-4
- Scalar theory showing that distorted quads already break serendipity elements

**Arnold, D.N., Boffi, D., and Falk, R.S.** "Quadrilateral H(div) finite elements".
*SIAM Journal on Numerical Analysis*, 42(6):2429-2451, 2005. DOI:10.1137/S0036142903431924
- The vector-valued case: RT[0] divergence does not converge on general quad meshes

**Falk, R.S., Gatto, P., and Monk, P.** "Hexahedral H(div) and H(curl) finite elements".
*ESAIM: Mathematical Modelling and Numerical Analysis*, 45(1):115-143, 2011.
DOI:10.1051/m2an/2010034
- Extends the negative results to 3D hexahedral H(curl), the case that bites h-φ directly

---

## Software Craft

The four books behind `doc/coding_philosophy.md` and `doc/commenting_guidelines.md`. Two are about
scientific code, two are general craft books. All four are in `literature/books/coding/` when that
directory is present; its `index.md` routes questions to chapter, section, Topic or Tip.

**Oliveira, S. and Stewart, D.E.** *Writing Scientific Software: A Guide to Good Style*. Cambridge
University Press, 2006. ISBN 978-0-521-85896-0. DOI: 10.1017/CBO9780511617973
- Cited for: comments and documentation (§7.9–7.10), interface design and scratch space (§7.6),
  memory management (Ch. 14–15), array layout and cache (§8.3, §12.5), testing and debugging (Ch. 9)
- Cite as: "Oliveira & Stewart 2006, §7.9"

**Rouson, D., Xia, J., and Xu, X.** *Scientific Software Design: The Object-Oriented Way*. Cambridge
University Press, 2011. ISBN 978-0-521-88813-4. DOI: 10.1017/CBO9780511977381
- Cited for: style rules (§1.7, Rules 1.1–1.3), comment intent as complexity grows (§7.2.1),
  design patterns for solvers (Ch. 5–9), C++/Fortran interoperability (Ch. 11)
- Cite as: "Rouson et al. 2011, §1.7"

**Martin, R.C.** *Clean Code: A Handbook of Agile Software Craftsmanship*, 2nd ed. Pearson, 2025.
ISBN 978-0-13-539857-9
- Cited for: naming (Ch. 4), comments (Ch. 5, pp. 89–116), the small-function position BELFEM
  declines inside numerical kernels (Ch. 7, Ch. 10), the debate with Ousterhout (Appendix)
- **Edition note:** BELFEM documents cite this second edition. The first edition (Prentice Hall,
  2008) numbers the comments chapter as Ch. 4; sections are unnumbered in both, so cite the chapter
  and the section title, with the printed page when precision matters
- Cite as: "Martin 2025, Ch. 5 'Journal Comments', p. 102"

**Thomas, D. and Hunt, A.** *The Pragmatic Programmer: Your Journey to Mastery*, 20th anniversary
ed. Addison-Wesley, 2020. ISBN 978-0-13-595705-9
- Cited for: documentation built in, non-API comments discuss *why* (Topic 7, Tip 13); DRY as a rule
  about knowledge, duplication in documentation (Topic 9); design by contract (Topic 23); crash
  early (Topic 24); assertions check the impossible and stay on (Topic 25); resource ownership
  (Topic 26); regression tests and automation (Topic 51)
- Cite by Topic (1–53) and Tip (1–100), never by page: the local extraction's page markers are
  reflow pages
- Cite as: "Thomas & Hunt 2020, Topic 25" or "Thomas & Hunt 2020, Tip 39"

---

## BELFEM Research Papers

### Tier 1: Core BELFEM

**Messe et al. 2023** — "BELFEM: a special purpose FE code for magnetodynamic modeling of HTS tapes" (primary reference)
- BELFEM architecture, h-φ formulation, thin-shell implementation
- Static condensation preferred over Lagrange multipliers
- Solver strategies: STRUMPACK vs MUMPS
- Nonlinear iteration: hybrid Picard/Quasi-Newton/Newton
- Checkerboarding prevention: tight tolerances (ε < 10⁻¹¹) — the paper's recommendation, tightened per deck; BELFEM's default `tolerance` is `1e-6`

**Messe 2022** — "A Special Purpose Finite-Element Framework for HTS Applications"
- Project motivation and scope

### Tier 2: h-φ Thin-Shell Focus

**Alves et al. 2022b** — "Thin-shell approach for modeling superconducting tapes in H-φ formulation"
- Full derivation with interface conditions
- Cohomology basis theory (thick cuts)
- N=1 vs N>1 analysis
- **Key finding:** N>1 needed for stacked tapes to capture top/bottom losses

**Alves et al. 2024** — "2-D Thin-Shell Model Based on H-φ-Formulation in COMSOL"
- Thin cut implementation
- Alternative approach to Alves et al. 2022b

**Schnaubelt et al. 2023** — "Electromagnetic Simulation of No-Insulation Coils Using H-φ TSA"
- Practical interpretation for coils
- Transport current handling

### Tier 3: Related Methods

**Arsenault et al. 2023** — "Magnetodynamic H-φ Formulation"
- **Recommended variant:** Magnetodynamic (H-φ/D) over quasi-static
- Interface coupling methodology
- Includes ∂(μh)/∂t term for correct transients

**Alves et al. 2022a** — "3D FE Thin-Shell Model for HTS Tapes"
- Adaptable to h-φ
- Relative homology background

**Arsenault et al. 2021** — "Implementation of H-φ Formulation in COMSOL Multiphysics"
- Performance comparison with pure H-formulation
- Element order recommendations (quadratic for 2D)

### Tier 4: Supporting Theory

**Dular et al. 2021** — "Stability of Mixed FE Formulations for HTS"
- Inf-sup conditions for h-φ
- Hierarchical basis functions
- HTS-specific stability analysis

**Riva et al. 2023** — "H-φ in Sparselizard with DDM"
- Domain decomposition methods
- **Note:** Different software (NOT BELFEM)

### Tier 5: HTS Material Models (E-J laws)

**Rhyner 1993** — "Magnetic properties and AC-losses of superconductors with power-law
current-voltage characteristics"
- Physica C 212, pp. 292–300
- The E-J power law underlying all three BELFEM resistivity laws
- Cited in `src/physics/materials/powerlaws.hpp`

**Plummer & Evetts 1987** — "Dependence of the shape of the resistive transition on composite
inhomogeneity in multifilamentary wires"
- IEEE Trans. Magn. 23 (2), pp. 1179–1182
- Phenomenology of the n-value characterizing the transition sharpness

**Duron et al. 2004** — "Modelling the E–J relation of high-Tc superconductors in an arbitrary
current range"
- J. Duron, F. Grilli, B. Dutoit, S. Stavrev — Physica C 401 (1-4), pp. 231–235
- The parallel combination ρ = ρ_PL·ρ_n/(ρ_PL+ρ_n) used by the `powerlaw` and `riva` laws
- **DOI:** 10.1016/j.physc.2003.09.044

**Riva 2021** — "Quench behavior of high-temperature superconductor tapes for power applications:
a strategy toward resilience"
- N. Riva, EPFL doctoral thesis no. 8754 (2021)
- Measured overcritical-regime resistivity of REBCO tapes 77–90 K (Ch. 3–5); the ρ_ηβ
  overcritical model and its continuity analysis (§5.1); impact of the law choice on simulated
  quench speed (Ch. 6)
- Namesake of the `resistivity type : riva` law (which implements the Duron parallel form the
  thesis builds on, not the fitted ρ_ηβ model itself)
- Local copy (ephemeral): `./tmp/EPFL_TH8754.pdf`

See `src/physics/materials/doc/resistivity_laws.md` for how these map onto the three laws.

### Tier 6: Elastic Constants of the Pure Metals

This section covers the `Metal::create_mech` closure and the datasets used to fit its constants
(`src/physics/materials/doc/materials_usage_guide.md` §9.7, `material_property_sources.md`).

**Garai & Laugier 2007** — "The temperature dependence of the isothermal bulk modulus at 1 bar pressure"
- J. Appl. Phys. 101, 023514
- B_T(T) = B_T(0)·exp(−∫α_V δ_T dT), Eq. 8–10, derived for the isothermal bulk modulus; BELFEM uses this
  form with constant δ for K and, as its own extension, for G
- **DOI:** 10.1063/1.2424535

**Ledbetter 1981** — "Elastic constants of polycrystalline copper at low temperatures"
- phys. stat. sol. (a) 66, pp. 477–484
- Polycrystalline Cu, 5–295 K: longitudinal and shear modulus, E, B, ν (Table 2); the copper anchor
- **DOI:** 10.1002/pssa.2210660209

**Ledbetter & Naimon 1974** — "Elastic properties of metals and alloys. II. Copper"
- J. Phys. Chem. Ref. Data 3, pp. 897–935
- Review; room-temperature polycrystal averages (Table 10)
- **DOI:** 10.1063/1.3253150

**Ledbetter & Reed 1973** — "Elastic properties of metals and alloys, I. Iron, nickel, and iron-nickel alloys"
- J. Phys. Chem. Ref. Data 2, pp. 531–618
- Compilation; §15 on the ΔE effect (saturated vs demagnetized moduli), the basis for Nickel's data choice
- **DOI:** 10.1063/1.3253127

**Alers, Neighbours & Sato 1960** — "Temperature dependent magnetic contributions to the high field elastic constants of nickel and an Fe-Ni alloy"
- J. Phys. Chem. Solids 13, pp. 40–55
- Ni single-crystal constants at 10 kOe, 0–760 K (Table 2); the data used for nickel
- **DOI:** 10.1016/0022-3697(60)90125-6

**Kamm & Alers 1964** — "Low-temperature elastic moduli of aluminum"
- J. Appl. Phys. 35, pp. 327–330 — Al single-crystal data 0–300 K (Table I)
- **DOI:** 10.1063/1.1713309

**Neighbours & Alers 1958** — "Elastic constants of silver and gold"
- Phys. Rev. 111, pp. 707–712 — Ag single-crystal data 0–300 K (Table I)
- **DOI:** 10.1103/PhysRev.111.707

**Waldorf & Alers 1962** — "Low-temperature elastic moduli of lead"
- J. Appl. Phys. 33, pp. 3266–3269 — Pb single-crystal data 0–300 K (Table I); interim source, see Lead's header
- **DOI:** 10.1063/1.1931149

**Rayne & Chandrasekhar 1960** — "Elastic constants of β tin from 4.2°K to 300°K"
- Phys. Rev. 120, pp. 1658–1663 — six tetragonal constants at 300, 77, 4.2 K (Table IV)
- **DOI:** 10.1103/PhysRev.120.1658

**Rayne & Chandrasekhar 1961** — "Elastic constants of iron from 4.2 to 300°K"
- Phys. Rev. 122, pp. 1714–1716 — smoothed zero-field constants (Table I)
- **DOI:** 10.1103/PhysRev.122.1714

**Palmer & Lee 1971** — "The elastic constants of chromium"
- Phil. Mag. 24, pp. 311–318 — zero-field constants 4.2–345 K through the spin-flip and Néel transitions
- **DOI:** 10.1080/14786437108227390

**Kim & Ledbetter 1998** — "Low-temperature elastic coefficients of polycrystalline indium"
- Mater. Sci. Eng. A 252, pp. 139–143 — Varshni fits of C_l, G, B, E, ν, 5–300 K (Table 1); the indium dataset
- **DOI:** 10.1016/S0921-5093(98)00490-0

---

## BELFEM FVM Papers

### Primary References

**Aavatsmark 2002** — "An Introduction to Multipoint Flux Approximations for Quadrilateral Grids"
- **DOI:** 10.1023/A:1021291114475
- **Use for:** MPFA O-method fundamentals, transmissibility matrices (start here)
- **Key topics:** Interaction volumes, T = CA⁻¹B - D, K-orthogonality, monotonicity

**Ingram et al. 2010** — "A Multipoint Flux Mixed Finite Element Method on Hexahedra"
- **DOI:** 10.1137/090766176
- **Use for:** Enhanced BDDF₁ for 3D hexahedra (critical for implementation)
- **Key topics:** 4 DOF/face, curl enrichment, trapezoidal quadrature, superconvergence

### Convergence and Stability

**Agélas et al. 2008** — "A Symmetric and Coercive Finite Volume Scheme"
- **DOI:** 10.1515/JNUM.2008.006
- **Use for:** Coercivity conditions, when your mesh will converge
- **Key topics:** Local coercivity coer(D,Λ), general meshes, L∞ coefficients

**Klausen & Winther 2006** — "Robust Convergence of Multi Point Flux Approximation on Rough Grids"
- **DOI:** 10.1007/s00211-006-0023-4
- **Use for:** Physical vs reference space, robust convergence on rough grids
- **Key topics:** Jacobian at cell center, stability condition (27), broken RT

### Theory and Connection

**Wheeler et al. 2011** — "A Family of Multipoint Flux Mixed Finite Element Methods"
- **DOI:** 10.1016/j.procs.2011.04.097
- **Use for:** MFMFE overview, variational framework
- **Key topics:** BDM₁/BDDF₁ with quadrature, FEM-FVM connection, h²-perturbations

---

## Computational Topology and Optimal Cuts

The references behind the `homology/` module: why cuts exist, how BELFEM computes them, and what is
and is not tractable if you want them minimal.

### Theory and Method

**Gross, P.W., and Kotiuga, P.R.** *Electromagnetic Theory and Computation: A Topological Approach*. MSRI Publications, vol. 48, Cambridge University Press, 2004. ISBN 0-521-80160-5.
- **Use for:** Why cuts must exist at all — cohomology for electromagnetics, branch cuts for the magnetic scalar potential, duality theorems
- **Cite as:** "Gross & Kotiuga, §X" or "Gross & Kotiuga, Ch. N"
- **Answers:** existence and construction in principle; not minimality

**Pellikka, M., Suuriniemi, S., Kettunen, L., and Geuzaine, C.** "Homology and Cohomology Computation in Finite Element Modeling", *SIAM Journal on Scientific Computing*, 2013.
- **DOI:** 10.1137/130906556
- **Use for:** **BELFEM's actual method** — chain complex → reduction → Smith normal form → thick cuts (the Gmsh solver)
- **Key topics:** cohomology basis on an FE mesh, thick cuts vs thin cuts

**Giard, G., et al.** "Generalized Pellikka algorithm for cohomology computation", **in preparation** (drafted 2025, shelved, resumed 2026).
- **DOI:** none — unpublished draft, not in `./literature/`
- **Use for:** the algorithm BELFEM actually runs by default — `reduce_complexPellikkaGeneralized` / `coreduce_complexPellikkaGeneralized`, selected at `cl_MaxwellFactory.hpp:48`
- **Key topics:** generalized Pellikka reduction and coreduction; the cohomology-side 0-cochain omit loop
- **Note:** this is the reason the cohomology core is closed to AI edits (`doc/ai_collaboration_protocol.md` §7.1). A modified algorithm with no published implementation gives a model no correct prior — its nearest neighbor is textbook Pellikka, so deliberate departures read as defects. Ask Gregory for the draft; do not reconstruct the method from the code.
- **Author's name is Giard.** The module documentation carried "Giarda" in nine places until 2026-08-31; if it reappears, it is a typo.

**Mrozek, M., and Batko, B.** "Coreduction Homology Algorithm", *Discrete & Computational Geometry*, 2009.
- **DOI:** 10.1007/s00454-008-9073-y
- **Use for:** Why the reduction step is fast — linear-time coreduction collapses the complex before Smith normal form

### Optimality and Complexity

**Dey, T.K., Hirani, A.N., and Krishnamoorthy, B.** "Optimal Homologous Cycles, Total Unimodularity, and Linear Programming", *SIAM Journal on Computing*, 2011.
- **DOI:** 10.1137/100800245
- **Use for:** When a *minimal* cut is computable in polynomial time — over **Z**, with a totally unimodular boundary matrix, via LP

**Chen, C., and Freedman, D.** "Hardness Results for Homology Localization", *Discrete & Computational Geometry*, 2011 (also Proc. ACM-SIAM SODA 2010).
- **DOI:** 10.1007/s00454-010-9322-8
- **Use for:** The negative result — over **Z₂** the problem is NP-hard to approximate within any constant factor. Read before proposing a "just minimize it" cut algorithm

**Dunfield, N.M., and Hirani, A.N.** "The Least Spanning Area of a Knot and the Optimal Bounding Chain Problem", *Proc. 27th Annual Symposium on Computational Geometry (SoCG '11)*, 2011.
- **DOI:** 10.1145/1998196.1998218
- **Use for:** Minimal spanning surfaces — NP-complete in general, polynomial when H₂ = 0

**Costantini, M.** "A Novel Phase Unwrapping Method Based on Network Programming", *IEEE Transactions on Geoscience and Remote Sensing*, 1998.
- **DOI:** 10.1109/36.673674
- **Use for:** The concrete grid analogue — branch cuts as min-cost network flow, i.e. a TU-structured LP solved at scale

**Haken, W.** "Theorie der Normalflächen: Ein Isotopiekriterium für den Kreisknoten", *Acta Mathematica*, 1961 (German).
- **DOI:** 10.1007/BF02559591
- **Use for:** Foundational — normal surface theory, decidability of unknot recognition. Historical context for why surface-finding is hard

### Algorithms Reference

**Cormen, T.H., Leiserson, C.E., Rivest, R.L., and Stein, C.** *Introduction to Algorithms*, 3rd ed. MIT Press, 2009. ISBN 978-0-262-03384-8.
- **Use for:** The graph and optimization algorithms underneath the above — BFS/DFS, spanning trees, shortest paths, min-cost flow, LP duality, total unimodularity, NP-completeness
- **Cite as:** "CLRS, §X.Y" or "Cormen et al., Ch. N"

---

## Key Implementation Guidelines from Literature

### FVM Implementations (from FVM papers)

**MPFA O-Method (Aavatsmark 2002 - Aavatsmark 2002)**
- Use interaction volumes at vertices with continuity points at subedge midpoints
- Compute transmissibility T = CA⁻¹B - D for local flux-pressure relationship
- Check K-orthogonality for TPFA feasibility

**Physical Space Evaluation (Klausen & Winther 2006 - Klausen 2006)**
- Evaluate Jacobian Dₖ = D(x̂ₖ) at cell CENTER, not at integration points
- Accept non-symmetric system matrix for robustness on rough grids
- Check stability condition det(Λ + Λᵀ) ≥ γ₀ > 0

**Enhanced BDDF₁ for 3D (Ingram et al. 2010 - Ingram 2010)**
- Standard BDDF₁ has only 3 DOF/face → insufficient for MPFA elimination
- Must use enhanced BDDF₁ with 6 curl terms → 4 DOF/face
- Requires h²-parallelepipeds for full O(h) convergence
- Superconvergence O(h²) at cell centers on regular grids

**Convergence Check (Agélas et al. 2008 - Agélas 2008)**
- Compute coercivity condition coer(D,Λ) ≥ θ > 0 before running
- Condition is computable for any mesh + tensor combination
- Predicts convergence on general polygonal/polyhedral meshes

### FEM Implementations (from FEM papers)

**Static Condensation (Messe et al. 2023)**
- **Preferred method** for hanging DOFs
- Avoids zero diagonal blocks in system matrix
- Better conditioning for direct solvers (MUMPS/STRUMPACK)

### Magnetodynamic Coupling (Arsenault et al. 2023)
- Use H-φ/D formulation (not quasi-static)
- Must include ∂(μh)/∂t term
- Better convergence than pure quasi-static

### Thin-Shell Element Order (Alves et al. 2022b, Arsenault et al. 2021)
- **N=1:** Linear through-thickness distribution (faster, similar to T-A)
- **N>1:** Needed when tangential field penetration matters
- Use N>1 for closely-packed tapes

### Convergence Tolerance (Messe et al. 2023)
- Drive ε < 10⁻¹¹ to prevent checkerboarding
- Tighter than typical FEM tolerances
- Critical for nonlinear HTS problems

### Solver Selection (Messe et al. 2023)
- **STRUMPACK:** Typically 2× faster
- **MUMPS:** More robust for difficult problems
- Try STRUMPACK first, fall back to MUMPS if needed

---

## Question Type → Literature Mapping

| Question | Primary Source | Secondary Sources |
|----------|----------------|-------------------|
| How does BELFEM implement X? | Messe et al. 2023, Messe 2022 | Specific papers |
| Theory behind h-φ formulation? | Arsenault et al. 2023, Arsenault et al. 2021 | Bathe Ch. 4, Hughes Ch. 3 |
| Why isn't this converging? | Messe et al. 2023 Section 2.7, Bathe Ch. 8 | Dular et al. 2021 (stability) |
| How do thin-shell models work? | Messe et al. 2023, Alves et al. 2022b | Alves et al. 2022a, Alves et al. 2024, Schnaubelt et al. 2023 |
| What are cuts for? | Alves et al. 2022b, Alves et al. 2024, Schnaubelt et al. 2023 | Alves et al. 2022a (homology) |
| How does MPFA work? | Aavatsmark 2002 | Wheeler et al. 2011 |
| Will my FVM mesh converge? | Agélas et al. 2008, Klausen & Winther 2006 | Aavatsmark 2002 (fundamentals) |
| How to implement 3D MPFA? | Ingram et al. 2010 | Aavatsmark 2002, Klausen & Winther 2006 |
| Physical vs reference space? | Klausen & Winther 2006 | Ingram et al. 2010 (3D details) |
| How to implement element X? | Hughes Ch. 5, Bathe Ch. 5 | Zienkiewicz Vol. 1 |
| What's the weak form for Y? | BELFEM papers | Bathe Ch. 3, Hughes Ch. 1 |
| Why is my element locking? | Hughes Ch. 4, Bathe §4.3-4.5 | Zienkiewicz Ch. 10 |
| How to prove convergence? | Brenner Ch. 4-5 | Bathe Ch. 4 |
| Which solver should I use? | Messe et al. 2023, Bathe Ch. 8 | Hughes Ch. 11 |
| How much Ic defect can a coil tolerate? | Wozniak et al. 2025 | Badel et al. 2019 |
| Will quench detection see a defect in time? | Wozniak et al. 2025 §IV-C | Badel et al. 2019 §5 |
| How do I control the time step in a nonlinear transient? | Badel et al. 2019 §4.3.2 (integral controller) | Bathe Ch. 9 |
| How do I warm-start Newton-Raphson between steps? | Badel et al. 2019 §4.3.1 | Messe et al. 2023 §3 |
| How is current shared between the tape and its stabilizer? | Badel et al. 2019 §3.1 | Wozniak et al. 2025 §III |

---

## Common Pitfalls and Solutions

| Pitfall | Symptom | Solution | Reference |
|---------|---------|----------|-----------|
| Loose convergence tolerance | Checkerboarding | Drive ε < 10⁻¹¹ | Messe et al. 2023 |
| N=1 for stacked tapes | Missing losses | Use N > 1 | Alves et al. 2022b |
| Identical FE orders at interface | Oscillations | Hierarchical enrichment | Dular et al. 2021 |
| Missing cuts | Wrong current flow | Add cuts (thin/thick) | Alves et al. 2022b, Alves et al. 2024 |
| Lagrange multipliers | Zero diagonal | Use static condensation | Messe et al. 2023 |
| QS formulation + linear elements | Poor convergence | Use magnetodynamic | Arsenault et al. 2023 |
| Neglecting ∂(μh)/∂t | Incorrect transient | Include time derivative | Arsenault et al. 2023, Alves et al. 2022b |
| TPFA on non-K-orthogonal grids | O(1) error | Use MPFA | Aavatsmark 2002 |
| Reference space on rough grids | Divergence | Use physical space Jacobian | Klausen & Winther 2006 |
| Standard BDDF₁ for hexahedra | Only 3 DOF/face | Use enhanced BDDF₁ | Ingram et al. 2010 |
| Violating coercivity | No convergence | Check coer(D,Λ) ≥ θ | Agélas et al. 2008 |

---

## Citation Guidelines

### When Referencing Papers in Code/Documentation

**DO:**
- Cite by **author and year**, with the exact section or equation: "Messe et al. 2023, §2.7", "Bathe, §8.6"
- Use full academic citations with author names, journal, volume, pages, year in reference lists
- Include DOIs as clickable links: `DOI: [10.xxxx/xxxxx](https://doi.org/10.xxxx/xxxxx)`
- Reference the published work directly
- Example: "Following the approach of Arsenault et al. 2023 (IEEE Trans. Appl. Supercond.)..."

**DO NOT:**
- Write a new `paperN` or `FN` alias — they were retired on 2026-08-11 and are decoded above for older records only
- Reference local file paths (e.g., "alves2022a.txt")
- Include extracted text content
- Assume literature files will be available in public repository

**Rationale:** Proper citations ensure academic attribution and provide public access via DOIs.

---

## How to Access These References

### Research Papers
Most BELFEM-related papers are published in:
- IEEE Transactions on Applied Superconductivity
- Superconductor Science and Technology
- COMPEL - The International Journal for Computation and Mathematics in Electrical and Electronic Engineering
- SIAM Journal on Scientific Computing

Search for papers using DOIs (when provided) or author names and keywords like:
- "h-phi formulation superconductor"
- "thin-shell HTS tape"
- "cohomology cuts transport current"

### Textbooks
Standard references available from:
- University libraries
- Publishers: Springer, Wiley, Dover, Prentice Hall
- Some older editions available as Dover reprints (Hughes, Zienkiewicz)
- Felippa's IFEM notes available online as course materials

---

## Workflow Examples

### Implementing a New Feature
1. Understand physics → BELFEM papers for HTS-specific, books for general FEM
2. Find weak form → Research papers
3. Check element technology → Hughes Ch. 5, Bathe Ch. 5
4. Verify stability → Dular et al. 2021 (h-φ specific) or Bathe/Brenner (general)
5. Implement → Cross-reference papers (equations) with books (algorithms)
6. Validate → Papers contain benchmark cases

### Debugging Convergence Issues
1. Identify symptom → Oscillations (Dular et al. 2021, Hughes Ch. 4) or Divergence (Bathe Ch. 8, Messe et al. 2023)
2. Check formulation → Interface coupling (Messe et al. 2023, Arsenault et al. 2023, Alves et al. 2022b)
3. Verify implementation → Compare with cited equations
4. Adjust solver strategy → Messe et al. 2023 (nonlinear), Bathe Ch. 8 (linear)

---

## Special Topics

### h-φ Formulation
Consult in this order:
1. Messe et al. 2023 — BELFEM's implementation choices
2. Arsenault et al. 2023 — Magnetodynamic coupling (recommended)
3. Arsenault et al. 2021 — Performance comparison
4. Alves et al. 2022b — Interface conditions and cuts
5. Bathe Ch. 4-5 — Mixed formulation theory
6. Dular et al. 2021 — Stability analysis

### Thin-Shell Implementation
1. Messe et al. 2023 — BELFEM's macro-element approach
2. Alves et al. 2022b — Full derivation
3. Alves et al. 2022a — 3D extension
4. Alves et al. 2024 — Alternative COMSOL approach
5. Hughes Ch. 5 — Isoparametric element foundations

### Cuts and Cohomology
1. Alves et al. 2022b — Cohomology basis theory (thick cuts)
2. Alves et al. 2024 — Thin cut implementation
3. Schnaubelt et al. 2023 — Practical interpretation for coils
4. Alves et al. 2022a — Relative homology background

---

## Tags for Searching

### BELFEM FEM
`#h-phi-formulation` `#thin-shell` `#cuts` `#cohomology` `#transport-current` `#hts-material` `#belfem-architecture` `#magnetodynamic` `#checkerboarding`

### BELFEM FVM
`#mpfa` `#transmissibility` `#coercivity` `#physical-space` `#enhanced-bddf1` `#k-orthogonal` `#mfmfe` `#anisotropic-diffusion`

### Computational Topology
`#cohomology-computation` `#thick-cuts` `#coreduction` `#smith-normal-form` `#total-unimodularity` `#np-hardness` `#min-cost-flow` `#normal-surfaces`

### General FEM
`#weak-form` `#isoparametric` `#locking` `#mixed-methods` `#inf-sup` `#newton-raphson` `#time-integration` `#eigenvalue` `#convergence`

---

## See Also

- [Documentation Guidelines](documentation_guidelines.md) - How to cite literature in BELFEM docs
- Project README: `../README.md`
- CLAUDE.md: `../CLAUDE.md` (includes literature directory information)

**Note:** If you have access to the full `./literature/` repository, see `./literature/README.md` for detailed navigation guides, routing tables, and search patterns for the complete text extractions.
