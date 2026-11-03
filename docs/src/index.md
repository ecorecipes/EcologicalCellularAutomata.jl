# EcologicalCellularAutomata.jl

A model zoo of published ecological cellular automata, from plant diseases and agricultural and forest pests to coastal dunes, each runnable on any compatible `CellularAutomataCore.jl` backend.

Part of the ecorecipes cellular-automata ecosystem: `CellularAutomataCore.jl` → `CategoricalCellularAutomata.jl` (applied category theory) and `CellularAutomataCore.jl` → `EcologicalCellularAutomata.jl` (model zoo).

Each model reimplements a published cellular automaton for a plant disease, a
vector of plant disease, an agricultural or forest pest, or a landform shaped by
vegetation (bio-geomorphology), as a backend-agnostic
`CellularAutomataCore` model, together with metadata recording the citation, the
mechanism, the validation reported by the source, and how faithfully the
reimplementation follows it. The literature these models come from is
reviewed in [Tonnang2017](@cite) among others.

## The models

| id | Source | Organism | Fidelity |
|---|---|---|---|
| `rhizoctonia_contact` | [Kleczkowski1997](@cite) | *Rhizoctonia solani* on radish | reproduction |
| `domany_kinzel_viral` | [Capitan2011](@cite) | a virus spreading through host tissue | reproduction |
| `fungus_dispersal` | [Ferreira2013](@cite) | a fungus among permissive and resistant hosts | reproduction |
| `huffaker_mites` | [Zsoldos2024](@cite) | predatory and prey mites on oranges | adaptation |
| `cvc_citrus` | [Martins2000](@cite) | citrus variegated chlorosis | adaptation |
| `coffee_rust` | [Beasley2022](@cite) | coffee leaf rust | adaptation |
| `bt_cotton_refuge` | [Malaquias2021](@cite) | two lepidopterans in a Bt cotton refuge | adaptation |
| `coffee_ant_fungus` | [Jackson2009](@cite) | ant nests and an entomopathogenic fungus | adaptation |
| `pine_wilt_casim` | [Huang2012](@cite) | pine wilt disease | adaptation |
| `grapevine_leafroll` | [Atallah2015](@cite) | grapevine leafroll with roguing and testing | adaptation |
| `cassava_mealybug` | [Promrak2017](@cite) | cassava mealybugs and green lacewings | adaptation |
| `telenomus_release` | [Weber2022](@cite) | releases of an egg parasitoid against a stink bug | reproduction |
| `bark_beetle_forest` | [Konigswieser2020](@cite) | a spruce forest and bark beetles | adaptation |
| `potato_moth_farmers` | [Rebaudo2011](@cite) | potato tuber moths spread by farmers | reproduction |
| `mpb_multiscale` | [Liu2025](@cite) | mountain pine beetles, within and between cells | reproduction |
| `western_x_orchard` | [Ghorbani2026](@cite) | Western X disease in a cherry orchard | adaptation |
| `sudden_oak_death` | [Meentemeyer2011](@cite) | sudden oak death in California | adaptation |
| `emerald_ash_borer` | [Mercader2011](@cite) | the emerald ash borer in ash stands | reproduction |
| `eab_management` | [Mercader2011b](@cite) | treating emerald ash borer outliers | reproduction |
| `ipydisp_dispersal` | [Agboka2024](@cite) | released biological control agents over a map | reproduction |
| `tuta_absoluta_asia` | [Guimapi2020](@cite) | the tomato leafminer across Asia | adaptation |
| `rootworm_rotation` | [Szalai2014](@cite) | the western corn rootworm and crop rotation | reproduction |
| `resistance_regional` | [Peck1997](@cite) | pesticide resistance across a region of fields | adaptation |
| `bt_crop_three_layer` | [Yang2009](@cite) | a Bt crop, a target and a secondary pest and their enemy | adaptation |
| `dune_formation` | [Teixeira2023](@cite) | coastal foredunes, sand and dune plants | adaptation |

A *reproduction* implements the source's model with its parameters; an
*adaptation* fills gaps or departs from it in ways its metadata lists. Each
model's vignette compares it with the source's results and states its
deviations.

```julia
using EcologicalCellularAutomata, CellularAutomataCore
list_models()
model_info(:fungus_dispersal)                 # citation, fidelity, deviations, ...
sol = solve(problem(:fungus_dispersal; seed=1))
model = build(:fungus_dispersal)              # the CAModel, for any scheme and backend
```

Some sources update individuals one at a time or need decisions on the whole
lattice (spraying when a pest passes a control level, spores walking between
patches). Those parts are callbacks stored with the problem, so `solve` runs
them, drawing their random numbers from the problem's seed.

## Errors

Every exception the package defines subtypes the abstract root
[`EcologicalCellularAutomataError`](@ref).

See the API Reference for docstrings and the [References](references.md) page
for full bibliographic entries.
