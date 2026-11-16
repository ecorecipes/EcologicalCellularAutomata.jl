# Plant diseases

Builders, problem constructors and helpers of the models whose `category` is "plant disease".

## Rhizoctonia solani: primary and secondary infection

```@docs
build_rhizoctonia_contact
problem_rhizoctonia_contact
rhizoctonia_inoculum
rhizoctonia_pathozone
rhizoctonia_days
rhizoctonia_infected
```

## Viral fronts: a multicomponent Domany–Kinzel automaton

```@docs
build_domany_kinzel_viral
problem_domany_kinzel_viral
viral_transition
viral_density
viral_mean_ability
viral_diversity
```

## Fungus dispersal among permissive and resistant hosts

```@docs
build_fungus_dispersal
problem_fungus_dispersal
fungus_kernel
fungus_hosts
fungus_infected
```

## Citrus variegated chlorosis

```@docs
build_cvc_citrus
problem_cvc_citrus
cvc_flight_kernel
cvc_incidence
cvc_field_curve
```

## Coffee leaf rust in fragmented landscapes

```@docs
build_coffee_rust
problem_coffee_rust
coffee_rust_spores!
coffee_rust_spores_archived!
coffee_rust_spore_callback
coffee_rust_landscape
random_cluster_landscape
random_element_landscape
classify_landscape
coffee_rust_prevalence
coffee_rust_rate
```
