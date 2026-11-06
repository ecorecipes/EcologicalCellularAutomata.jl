"""
    EcologicalCellularAutomata

A model zoo of published ecological cellular automata, from plant diseases and agricultural and forest pests to coastal dunes, each runnable on any compatible `CellularAutomataCore.jl` backend.

Each model reimplements a published cellular automaton for a plant disease, a
vector of plant disease, an agricultural or forest pest, or a landform shaped by
vegetation (bio-geomorphology), as a backend-agnostic
`CellularAutomataCore` model, together with metadata recording the citation, the
mechanism, the validation reported by the source, and how faithfully the
reimplementation follows it. The literature these models come from is
reviewed in [Tonnang2017](@cite) among others.

Part of the ecorecipes cellular-automata ecosystem.
"""
module EcologicalCellularAutomata

using CellularAutomataCore
using CellularAutomataCore: CellularAutomataCore as CAC
using FFTW: FFTW
using Random: Random
using TOML: TOML

# The registry
export ModelInfo, list_models, model_info, build, problem
# Models: builders, problem constructors and model-specific helpers
export build_fungus_dispersal, problem_fungus_dispersal, fungus_kernel, fungus_hosts,
       fungus_infected
export build_rhizoctonia_contact, problem_rhizoctonia_contact, rhizoctonia_inoculum,
       rhizoctonia_pathozone, rhizoctonia_days, rhizoctonia_infected
export build_domany_kinzel_viral, problem_domany_kinzel_viral, viral_density,
       viral_mean_ability, viral_diversity, viral_transition
export build_bt_cotton_refuge, problem_bt_cotton_refuge, bt_cotton_survival,
       bt_cotton_spray_callback, bt_cotton_dominance
export build_coffee_rust, problem_coffee_rust, coffee_rust_spores!,
       coffee_rust_spores_archived!,
       coffee_rust_spore_callback,
       classify_landscape, random_cluster_landscape, random_element_landscape,
       coffee_rust_landscape, coffee_rust_prevalence, coffee_rust_rate
export build_huffaker_mites, problem_huffaker_mites, huffaker_step!, huffaker_step_callback,
       HUFFAKER_POSTS, huffaker_counts, count_peaks, huffaker_outcome
export build_cvc_citrus, problem_cvc_citrus, cvc_flight_kernel, cvc_incidence,
       cvc_field_curve
export build_coffee_ant_fungus, problem_coffee_ant_fungus, ant_iterations, nest_density,
       mean_variance_ratio
export build_pine_wilt_casim, problem_pine_wilt_casim, pine_wilt_weights, pine_landscape,
       pine_wilt_infected
export build_grapevine_leafroll, problem_grapevine_leafroll, grapevine_manage!,
       grapevine_strategy_callback, grapevine_revenue, grapevine_npv, grapevine_half_life
export build_cassava_mealybug, problem_cassava_mealybug, cassava_sweep!, cassava_release!,
       cassava_management_callback, cassava_yield
export build_telenomus_release, problem_telenomus_release, StinkBugEggMass,
       telenomus_dispersal_steps, telenomus_step!, telenomus_callback,
       telenomus_release_cells, telenomus_et, telenomus_counts
export build_bark_beetle_forest, problem_bark_beetle_forest, bark_beetle_step!,
       bark_beetle_callback, bark_beetle_counts, bark_beetle_growth_rate
export build_potato_moth_farmers, problem_potato_moth_farmers, potato_moth_survival,
       potato_moth_farmers_step!, potato_moth_farmers_callback, potato_moth_landscape,
       potato_moth_counts, potato_moth_invasion_time
export build_mpb_multiscale, problem_mpb_multiscale, mpb_disperse!, mpb_callback,
       mpb_pattern, mpb_complexity, mpb_radius
export build_western_x_orchard, problem_western_x_orchard, western_x_manage!,
       western_x_callback, western_x_economics
export build_sudden_oak_death, problem_sudden_oak_death, sod_kernel, sod_temperature_index,
       sod_weather, sod_infected_area
export build_emerald_ash_borer, problem_emerald_ash_borer, eab_disperse!, eab_callback,
       eab_environment, eab_spread
export build_eab_management, problem_eab_management, eab_manage!, eab_management_callback,
       eab_attraction, eab_beyond
export build_ipydisp_dispersal, problem_ipydisp_dispersal, ipydisp_distance, ipydisp_grid,
       ipydisp_neighbours, ipydisp_suitable, ipydisp_example, ipydisp_counts
export build_dune_formation, problem_dune_formation, dune_step!, dune_callback, dune_runup,
       dune_profile, dune_landform, dune_water_levels, scheveningen_water_levels, dune_toe,
       dune_measures, dune_wind
export ProfileSegment, shoreface_segment, slope_segment, flat_segment, ridge_segment,
       profile_heights, profile_kinds, reference_profile
# Exceptions: the root and the concrete types under it
export build_tuta_absoluta_asia, problem_tuta_absoluta_asia, tuta_climate, tuta_cell,
       tuta_arrival
export build_rootworm_rotation, problem_rootworm_rotation, rootworm_probabilities,
       rootworm_rings, rootworm_step!, rootworm_callback, rootworm_risk, rootworm_output
export build_resistance_regional, problem_resistance_regional, resistance_summary,
       resistance_time
export build_bt_crop_three_layer, problem_bt_crop_three_layer, bt_crop_summary
export EcologicalCellularAutomataError, UnknownModelError, MetadataError

include("errors.jl")
include("registry.jl")
include("landforms.jl")
include("models/rhizoctonia_contact.jl")
include("models/huffaker_mites.jl")
include("models/cvc_citrus.jl")
include("models/fungus_dispersal.jl")
include("models/domany_kinzel_viral.jl")
include("models/coffee_rust.jl")
include("models/bt_cotton_refuge.jl")
include("models/coffee_ant_fungus.jl")
include("models/pine_wilt_casim.jl")
include("models/grapevine_leafroll.jl")
include("models/cassava_mealybug.jl")
include("models/telenomus_release.jl")
include("models/bark_beetle_forest.jl")
include("models/potato_moth_farmers.jl")
include("models/mpb_multiscale.jl")
include("models/western_x_orchard.jl")
include("models/sudden_oak_death.jl")
include("models/emerald_ash_borer.jl")
include("models/eab_management.jl")
include("models/ipydisp_dispersal.jl")
include("models/dune_formation.jl")
include("models/tuta_absoluta_asia.jl")
include("models/rootworm_rotation.jl")
include("models/resistance_regional.jl")
include("models/bt_crop_three_layer.jl")

end # module
