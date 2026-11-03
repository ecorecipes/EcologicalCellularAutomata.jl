# Bio-geomorphology

Builders, problem constructors and helpers of the models whose `category` is "bio-geomorphology": landforms
shaped by vegetation.

## Coastal dunes

```@docs
build_dune_formation
problem_dune_formation
dune_step!
dune_callback
dune_runup
dune_profile
dune_landform
dune_water_levels
dune_toe
dune_measures
dune_wind
```

## Landform profiles

Cross-shore profiles built by joining pieces from the sea to the land; two pieces join where the first ends
at the elevation where the second starts.

```@docs
ProfileSegment
Base.:*(::ProfileSegment, ::ProfileSegment)
shoreface_segment
slope_segment
flat_segment
ridge_segment
profile_heights
profile_kinds
reference_profile
```
