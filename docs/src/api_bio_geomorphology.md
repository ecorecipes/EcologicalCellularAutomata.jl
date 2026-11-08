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
scheveningen_water_levels
dune_toe
dune_measures
dune_wind
```

## Barrier islands

```@docs
build_barrier_island
problem_barrier_island
barrier_step!
barrier_callback
barrier_event!
barrier_crest
barrier_profile
barrier_storms
barrier_storm_series
barrier_florence_site
barrier_shoreface
barrier_diffusivity
barrier_measures
brier_skill_score
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
