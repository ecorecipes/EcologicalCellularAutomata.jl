# Flavescence dorée at Faleyras, Adrakey et al. (2023)

`vines.csv` holds the 6,056 vines of three plots at Faleyras (Gironde, France) surveyed for flavescence dorée in
2018 and 2019, from the deposit of Adrakey, Gibson, Eveillard, Malembic-Maher and Fabre (2023), "Bayesian
inference for spatio-temporal stochastic transmission of plant disease in the presence of roguing: a case study to
characterise the dispersal of Flavescence dorée", PLOS Computational Biology 19(9), e1011399
(doi 10.1371/journal.pcbi.1011399); data deposited at Recherche Data Gouv, doi 10.57745/YXOEHX, under the Licence
Ouverte / Open Licence 2.0 (Etalab).

Columns: `x`, `y` (m, east and north of the westernmost and southernmost vine; the source gives RGF93 Lambert-93
coordinates, from which 445553.052 m east and 6411006.086 m north were subtracted), `plot` (CS, M1, M2), `cultivar`
(`cabernet_sauvignon`, `merlot`), `status` and `titre` (the phytoplasma titre class 1-4 of a symptomatic vine, 0
otherwise). `status`, from the source's `val` and `detection` columns: `absent_2018` (no vine at the 2018 survey),
`symptomatic_2018`, `absent_2019` (present in 2018 without symptoms, gone by the 2019 survey; in plot CS, which was
uprooted after 2018, every vine not symptomatic in 2018), `symptomatic_2019` and `healthy` (never seen with
symptoms). The data are otherwise unchanged.

This file is distributed under the Licence Ouverte / Open Licence 2.0, not under this package's licence.
