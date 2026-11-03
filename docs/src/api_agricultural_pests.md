# Agricultural pests

Builders, problem constructors and helpers of the models whose `category` is "agricultural pest".

## Huffaker's mites

```@docs
build_huffaker_mites
problem_huffaker_mites
huffaker_step!
huffaker_step_callback
HUFFAKER_POSTS
huffaker_counts
count_peaks
huffaker_outcome
```

## Two lepidopterans in a Bt cotton refuge

```@docs
build_bt_cotton_refuge
problem_bt_cotton_refuge
bt_cotton_survival
bt_cotton_spray_callback
bt_cotton_dominance
```

## Ant nests and a fungus on a coffee farm

```@docs
build_coffee_ant_fungus
problem_coffee_ant_fungus
ant_iterations
nest_density
mean_variance_ratio
```

## Cassava mealybugs and green lacewings

```@docs
build_cassava_mealybug
problem_cassava_mealybug
cassava_sweep!
cassava_release!
cassava_management_callback
cassava_yield
```

## Telenomus releases against the brown stink bug

```@docs
build_telenomus_release
problem_telenomus_release
StinkBugEggMass
telenomus_dispersal_steps
telenomus_step!
telenomus_callback
telenomus_release_cells
telenomus_et
telenomus_counts
```

## Potato tuber moths and farmers

```@docs
build_potato_moth_farmers
problem_potato_moth_farmers
potato_moth_survival
potato_moth_farmers_step!
potato_moth_farmers_callback
potato_moth_landscape
potato_moth_counts
potato_moth_invasion_time
```

## A Bt crop, two pests and their enemy on three layers

```@docs
build_bt_crop_three_layer
problem_bt_crop_three_layer
bt_crop_summary
```
