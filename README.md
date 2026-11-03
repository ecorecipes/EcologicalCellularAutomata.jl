# EcologicalCellularAutomata.jl

[![Build Status](https://github.com/ecorecipes/EcologicalCellularAutomata.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/ecorecipes/EcologicalCellularAutomata.jl/actions/workflows/CI.yml)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://ecorecipes.github.io/EcologicalCellularAutomata.jl/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A model zoo of published ecological cellular automata, from plant diseases and agricultural and forest pests to coastal dunes, each runnable on any compatible `CellularAutomataCore.jl` backend.

Part of the ecorecipes cellular-automata ecosystem: `CellularAutomataCore.jl` → `CategoricalCellularAutomata.jl` (applied category theory) and `CellularAutomataCore.jl` → `EcologicalCellularAutomata.jl` (model zoo).

## Features

- A registry of models (`list_models`, `model_info`, `build`, `problem`), each with metadata recording its
  source, what the source validated it against, its fidelity (`reproduction`, `adaptation` or `inspired`) and
  every deviation from the source.
- Twenty-five models from the plant-health and bio-geomorphology literature, checked against their primary papers:
  - Rhizoctonia damping-off with primary and secondary infection;
  - a multicomponent Domany-Kinzel automaton for viral fronts;
  - fungus dispersal among permissive and resistant hosts;
  - predator and prey mites in Huffaker's orange universe;
  - citrus variegated chlorosis with Levy-flight vectors;
  - coffee leaf rust in fragmented landscapes, with neutral landscape generators;
  - two competing lepidopterans in a Bt cotton refuge with threshold spraying;
  - ant nests on a coffee farm and an entomopathogenic fungus;
  - pine wilt disease on a forest raster, with felling;
  - grapevine leafroll in a vineyard, with roguing, testing and net present values;
  - cassava mealybugs and the release of green lacewings;
  - releases of the egg parasitoid *Telenomus podisi* against the brown stink bug in soybean;
  - a spruce forest growing, harvested and struck by drought and bark beetle outbreaks;
  - potato tuber moths in an Andean valley, carried between villages by farmers;
  - mountain pine beetles on a grid of hexagons, with population dynamics within cells and pooled dispersal between them;
  - Western X disease in a cherry orchard, with removal, replanting, testing and net present values;
  - sudden oak death spreading through forests by short- and long-range spores, with seasons and weather;
  - the emerald ash borer spreading through ash, by its foraging and the arrangement of its host;
  - slowing emerald ash borer outliers by removing ash, girdling trap trees or insecticide;
  - the dispersal of released biological control agents over a map, gated by the environment (iPyDisp);
  - the tomato leafminer *Tuta absoluta* spreading across Asia from Turkey, gated by monthly temperature and humidity;
  - the western corn rootworm in a region of rotated and continuous maize, and laws making rotation compulsory;
  - the evolution of pesticide resistance across a region of sprayed fields and refuges, with economic thresholds;
  - a Bt crop with its target pest, a secondary pest and their natural enemy on three stacked lattices;
  - a coastal foredune growing from wind-blown sand trapped by plants, with the sea resetting the beach.
- Every model is a `CellularAutomataCore` `CAModel`, so it runs with any scheme and backend that supports it.

## Installation

The ecosystem packages are not registered, and this package finds its siblings through `[sources]` paths, so
clone them side by side:

```sh
mkdir -p ca && cd ca
git clone https://github.com/ecorecipes/CellularAutomataCore.jl
git clone https://github.com/ecorecipes/EcologicalCellularAutomata.jl
cd EcologicalCellularAutomata.jl
julia --project -e 'using Pkg; Pkg.instantiate()'
```

Requires Julia ≥ 1.12.

## Vignettes

| Vignette | Topic |
|---|---|
| [01 Using the zoo](vignettes/01_using_the_zoo/01_using_the_zoo.md) | listing models, reading metadata, building and running models on several backends |
| [02 Rhizoctonia](vignettes/02_rhizoctonia_contact/02_rhizoctonia_contact.md) | primary and secondary infection; mean and variance against microcosm data |
| [03 Viral fronts](vignettes/03_domany_kinzel_viral/03_domany_kinzel_viral.md) | extinction thresholds and the saturation of replicative ability |
| [04 Fungus dispersal](vignettes/04_fungus_dispersal/04_fungus_dispersal.md) | resistant hosts, their arrangement, propagule load and seasonality |
| [05 Huffaker's mites](vignettes/05_huffaker_mites/05_huffaker_mites.md) | predator-prey oscillations, patches and posts |
| [06 Citrus variegated chlorosis](vignettes/06_cvc_citrus/06_cvc_citrus.md) | Levy flights, stress and seasons; calibrating to the field |
| [07 Coffee rust](vignettes/07_coffee_rust/07_coffee_rust.md) | neutral landscapes, spores and canopy, clustering |
| [08 Bt cotton refuge](vignettes/08_bt_cotton_refuge/08_bt_cotton_refuge.md) | competition, threshold spraying and resistance scenarios |
| [09 Ants and a fungus](vignettes/09_coffee_ant_fungus/09_coffee_ant_fungus.md) | nest clusters, update order, the fungus as density-dependent killer |
| [10 Pine wilt](vignettes/10_pine_wilt_casim/10_pine_wilt_casim.md) | ring-weighted neighbourhoods, felling rates, a stand-in forest |
| [11 Grapevine leafroll](vignettes/11_grapevine_leafroll/11_grapevine_leafroll.md) | roguing and testing strategies ranked by net present value |
| [12 Cassava mealybug](vignettes/12_cassava_mealybug/12_cassava_mealybug.md) | within-plant life cycles, surveys and lacewing releases |
| [13 Telenomus releases](vignettes/13_telenomus_release/13_telenomus_release.md) | an individual-based lattice of stink bugs and wasps; release strategies against the economic threshold |
| [14 Bark beetles in a spruce forest](vignettes/14_bark_beetle_forest/14_bark_beetle_forest.md) | a forest that grows from a few trees; contamination, strip harvest, drought, and the timing of a bark beetle outbreak |
| [15 Potato moths and farmers](vignettes/15_potato_moth_farmers/15_potato_moth_farmers.md) | a moth CA coupled to farmer agents; how travel, numbers and pest control knowledge speed an invasion |
| [16 Mountain pine beetles across scales](vignettes/16_mpb_multiscale/16_mpb_multiscale.md) | a multiscale CA on hexagons; parameter sweeps, regimes, seasons and barriers |
| [17 Western X disease in a cherry orchard](vignettes/17_western_x_orchard/17_western_x_orchard.md) | a five-state orchard CA; calibrating unreported rates, and the value of removal, replanting and testing |
| [18 Sudden oak death](vignettes/18_sudden_oak_death/18_sudden_oak_death.md) | spore dispersal by two Cauchy kernels as FFT kernel fields; weather scenarios and curbing long-distance spread |
| [19 Emerald ash borer](vignettes/19_emerald_ash_borer/19_emerald_ash_borer.md) | a coupled map lattice of phloem and beetles; foraging behaviour and the arrangement of ash |
| [20 Coastal dunes](vignettes/20_dune_formation/20_dune_formation.md) | a sand-slab automaton with vegetation and the sea; how far slabs hop decides how a foredune grows |
| [21 Managing emerald ash borer outliers](vignettes/21_eab_management/21_eab_management.md) | ash removal, girdled trap trees and insecticide around an outlier; the radius, timing and intensity of treatment |
| [22 Dispersal of biological control agents](vignettes/22_ipydisp_dispersal/22_ipydisp_dispersal.md) | iPyDisp's automaton on a graph of map cells; travel distance, environmental constraints and arrival times |
| [23 The tomato leafminer across Asia](vignettes/23_tuta_absoluta_asia/23_tuta_absoluta_asia.md) | a three-state invasion on monthly climate; exposed and invaded cells, the reach per month, and records of the pest |
| [24 Crop rotation against the western corn rootworm](vignettes/24_rootworm_rotation/24_rootworm_rotation.md) | random crop change that keeps the shares of maize, density-dependent dispersal on a torus, and mandatory rotation |
| [25 Pesticide resistance across a region](vignettes/25_resistance_regional/25_resistance_regional.md) | logistic growth, selection and economic thresholds in each field, migration between fields, and patches of resistance |
| [26 A Bt crop, two pests and their enemy](vignettes/26_bt_crop_three_layer/26_bt_crop_three_layer.md) | three stacked layers, moves and births in contiguous pairs as a Margolus block rule, Penna individuals, and the release of a secondary pest |
| [27 Dune landforms](vignettes/27_dune_landforms/27_dune_landforms.md) | the dune model on landforms the source never ran, built from the same profile pieces: an embryo dune, a seawall, a breached foredune, and one's own |
| [28 Wind over dunes](vignettes/28_wind_over_dunes/28_wind_over_dunes.md) | the dune model with a wind steered by the topography (linear theory of flow over low hills): speed-up on dune faces, funnelling into a breach, narrower dunes, and why a breach still fills |

Each vignette is also published as a tutorial in the
[documentation](https://ecorecipes.github.io/EcologicalCellularAutomata.jl/); the sources are the `.qmd` files in
[`vignettes/`](vignettes/).

## References

The full bibliography for the ecosystem is `docs/src/references.bib` (also in
`vignettes/references.bib`); the rendered list is the References page of the
[documentation](https://ecorecipes.github.io/EcologicalCellularAutomata.jl/).
