# Forest pests

Builders, problem constructors and helpers of the models whose `category` is "forest pest".

## Pine wilt disease

```@docs
build_pine_wilt_casim
problem_pine_wilt_casim
pine_wilt_weights
pine_landscape
pine_wilt_infected
```

## A spruce forest and bark beetles

```@docs
build_bark_beetle_forest
problem_bark_beetle_forest
bark_beetle_step!
bark_beetle_callback
bark_beetle_counts
bark_beetle_growth_rate
```

## Mountain pine beetles across scales

```@docs
build_mpb_multiscale
problem_mpb_multiscale
mpb_disperse!
mpb_callback
mpb_pattern
mpb_complexity
mpb_radius
```

## Emerald ash borer

```@docs
build_emerald_ash_borer
problem_emerald_ash_borer
eab_disperse!
eab_callback
eab_environment
eab_spread
```

## Managing emerald ash borer outliers

```@docs
build_eab_management
problem_eab_management
eab_manage!
eab_management_callback
eab_attraction
eab_beyond
```
