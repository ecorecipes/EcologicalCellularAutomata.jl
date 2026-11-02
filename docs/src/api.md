# Registry and exceptions

```@docs
EcologicalCellularAutomata
```

## The registry

Every model has a builder, which returns a backend-agnostic `CAModel`, a
problem constructor with the setup of its source, and metadata in
`models/<id>/metadata.toml`, read and checked when the package is compiled.

```@docs
ModelInfo
list_models
model_info
build
problem
```

Each model's own functions are listed by category: [plant diseases](api_plant_diseases.md), [agricultural pests](api_agricultural_pests.md), [forest pests](api_forest_pests.md) and [bio-geomorphology](api_bio_geomorphology.md).

## Exceptions

Every exception the package defines subtypes the abstract root
[`EcologicalCellularAutomataError`](@ref). Invalid arguments and keywords raise
Base's `ArgumentError`, which is outside it.

```@docs
EcologicalCellularAutomataError
UnknownModelError
MetadataError
```
