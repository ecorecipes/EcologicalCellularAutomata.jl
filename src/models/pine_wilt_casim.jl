# Huang, Yu, Zhou & Bei (2012): PWD-CASim, pine wilt disease (Bursaphelenchus xylophilus, carried by
# Monochamus beetles) spreading through a pine forest on a 10 m raster. A healthy pine is infected
# with probability sum_n w_n C_n / (8n), where C_n is the number of infected pines in the square ring
# at distance n; infections may pass through a latent period; symptomatic pines are felled with
# probability alpha a year. The source gives no rule for natural death; here a symptomatic pine dies
# after four years.

const PWD_NONFOREST = Int8(0)
const PWD_OTHER = Int8(1)
const PWD_HEALTHY = Int8(2)
const PWD_LATENT = Int8(3)
const PWD_SYMPTOMATIC = Int8(4)
const PWD_DEAD = Int8(5)

"""
    pine_wilt_weights(k) -> Vector{Float64}

The ring weights of the source's Moore-`k` neighbourhoods, `k` in 1:5, as read
from its simulation tool: `[1]`, `[1, 0.8]`, `[1, 0.8, 0.6]`,
`[1, 0.8, 0.6, 0.4]` and `[1, 0.8, 0.5, 0.3, 0.1]`.
"""
function pine_wilt_weights(k::Integer)
    k in 1:5 ||
        throw(ArgumentError("the source gives ring weights for Moore 1 to 5, got $k"))
    return ([1.0], [1.0, 0.8], [1.0, 0.8, 0.6], [1.0, 0.8, 0.6, 0.4],
            [1.0, 0.8, 0.5, 0.3, 0.1])[k]
end

"""
    build_pine_wilt_casim(; weights=pine_wilt_weights(1), incubation=false) -> CAModel

The pine wilt disease automaton of [Huang2012](@cite), in yearly steps. Each
cell is a 10 m pixel: non-forest (0), other forest (1), healthy pine (2),
latent (3) or symptomatic (4) infected pine, or dead pine (5); the layer
`years` counts years in the latent or symptomatic state.

- A healthy pine is infected with probability `min(1, Σₙ wₙ Cₙ / (8n))`, where
  `Cₙ` is the number of latent or symptomatic pines in the square ring at
  distance `n` (8n cells) and `weights` gives `wₙ` for rings `1:length(weights)`.
  With `incubation = true` a new infection is latent with probability
  `p_latent` and symptomatic otherwise; without, it is symptomatic.
- A latent pine stays latent with probability `p_stay`, for at most
  `max_latent` years, and otherwise shows symptoms.
- A symptomatic pine is felled with probability `α` a year and becomes other
  forest; otherwise it dies after `max_symptomatic` years.

Other cells do not change. Nothing lies beyond a landscape: edges reflect.
"""
function build_pine_wilt_casim(; weights::AbstractVector{<:Real}=pine_wilt_weights(1),
                               incubation::Bool=false)
    k = length(weights)
    k >= 1 || throw(ArgumentError("give at least one ring weight"))
    all(>=(0), weights) || throw(ArgumentError("ring weights must not be negative"))
    offsets = [CartesianIndex(i, j) for j in (-k):k for i in (-k):k if (i, j) != (0, 0)]
    ring = [max(abs(o[1]), abs(o[2])) for o in offsets]
    scale = [weights[r] / (8r) for r in ring]
    rule = LocalRule(Offsets(offsets)) do cell, nbrs, ctx
        p = ctx.p
        s = cell.cover
        if s == PWD_HEALTHY
            prob = 0.0
            for m in eachindex(scale)
                c = nbrs[m].cover
                (c == PWD_LATENT || c == PWD_SYMPTOMATIC) && (prob += scale[m])
            end
            bernoulli(ctx, min(1.0, prob), 1) || return cell
            latent = incubation && bernoulli(ctx, p.p_latent, 2)
            return (cover=latent ? PWD_LATENT : PWD_SYMPTOMATIC, years=Int8(0))
        elseif s == PWD_LATENT
            if cell.years + 1 <= p.max_latent && bernoulli(ctx, p.p_stay, 3)
                return (cover=PWD_LATENT, years=cell.years + Int8(1))
            end
            return (cover=PWD_SYMPTOMATIC, years=Int8(0))
        elseif s == PWD_SYMPTOMATIC
            bernoulli(ctx, p.α, 4) && return (cover=PWD_OTHER, years=Int8(0))
            years = cell.years + Int8(1)
            return years >= p.max_symptomatic ? (cover=PWD_DEAD, years=Int8(0)) :
                   (cover=PWD_SYMPTOMATIC, years)
        end
        return cell
    end
    return CAModel(rule; name=:pine_wilt_casim,
                   states=[StateVariable(:cover, PWD_HEALTHY; domain=Int8(0):Int8(5),
                                         description="0 non-forest, 1 other forest, 2 healthy, 3 latent, 4 symptomatic, 5 dead pine"),
                           StateVariable(:years, Int8(0);
                                         description="years latent or symptomatic")],
                   parameters=[Parameter(:α, 0.6;
                                         description="felling rate of symptomatic pines per year"),
                               Parameter(:p_latent, 0.8;
                                         description="share of new infections that are latent (incubation)"),
                               Parameter(:p_stay, 0.1;
                                         description="probability a latent pine stays latent a year"),
                               Parameter(:max_latent, 2;
                                         description="longest latent period (years)"),
                               Parameter(:max_symptomatic, 4;
                                         description="years from symptoms to natural death")],
                   description="Huang et al. (2012) PWD-CASim, Moore $k, incubation $incubation")
end

"""
    pine_landscape(rng; n=300, pine=0.6, other=0.25, clustering=0.3) -> Matrix{Int8}

A substitute for the source's forest map: an `n × n` neutral landscape
([`random_cluster_landscape`](@ref)) classified into pine (2), other forest
(1) and non-forest (0) in the given shares.
"""
function pine_landscape(rng; n::Integer=300, pine::Real=0.6, other::Real=0.25,
                        clustering::Real=0.3)
    0 < pine && 0 <= other && pine + other <= 1 ||
        throw(ArgumentError("pine and other forest shares must be positive and sum to at most 1"))
    classes = classify_landscape(random_cluster_landscape(rng, n, clustering),
                                 [pine, other, 1 - pine - other])
    return Int8[c == 0 ? PWD_HEALTHY : c == 1 ? PWD_OTHER : PWD_NONFOREST for c in classes]
end

"""
    problem_pine_wilt_casim(; cover=nothing, n=300, moore=1, weights=pine_wilt_weights(moore),
                            incubation=false, symptomatic=1400, latent=200, tspan=10,
                            seed=0, pmap=(;), kwargs...) -> CAProblem

The setup of [Huang2012](@cite): 1400 symptomatic and 200 latent pines placed
at random among the pines of a landscape, for 10 years. `cover` is a matrix of
land cover (0 non-forest, 1 other forest, 2 pine), for example from a raster;
by default a [`pine_landscape`](@ref) of size `n`, drawn from `seed`, stands in
for the source's unpublished map.
"""
function problem_pine_wilt_casim(; cover=nothing, n::Integer=300, moore::Integer=1,
                                 weights=pine_wilt_weights(moore), incubation::Bool=false,
                                 symptomatic::Integer=1400, latent::Integer=200, tspan=10,
                                 seed::Integer=0, pmap=(;), kwargs...)
    rng = Random.Xoshiro(seed)
    land = cover === nothing ? pine_landscape(rng; n) : Int8.(cover)
    all(in((PWD_NONFOREST, PWD_OTHER, PWD_HEALTHY)), land) ||
        throw(ArgumentError("cover must hold only 0 (non-forest), 1 (other forest) and 2 (pine)"))
    pines = findall(==(PWD_HEALTHY), land)
    symptomatic + latent <= length(pines) ||
        throw(ArgumentError("$(symptomatic + latent) infected pines do not fit among $(length(pines)) pines"))
    picks = pines[Random.randperm(rng, length(pines))[1:(symptomatic + latent)]]
    land[picks[1:symptomatic]] .= PWD_SYMPTOMATIC
    land[picks[(symptomatic + 1):end]] .= PWD_LATENT
    space = Grid(size(land)...; boundary=Reflecting())
    return CAProblem(build_pine_wilt_casim(; weights, incubation), (; cover=land), tspan,
                     pmap; space, seed, kwargs...)
end

"""
    pine_wilt_infected(u) -> Int

The number of infected pines, latent or symptomatic.
"""
pine_wilt_infected(u) = count(c -> c == PWD_LATENT || c == PWD_SYMPTOMATIC, u.cover)

_register!(:pine_wilt_casim; build=build_pine_wilt_casim, problem=problem_pine_wilt_casim)
