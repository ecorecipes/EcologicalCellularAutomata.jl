# Szalai, Kiss, Kövér & Toepfer (2014): the western corn rootworm, Diabrotica virgifera virgifera, in a region
# of maize and other crops, and the effect of crop rotation and of laws making it compulsory. Each cell is a
# field, in another crop, first-year maize or continuous maize; crops change at random each year with
# probabilities that keep their shares. The pest grows only in continuous maize; adults emigrate from it, in
# proportion to their density, to the maize fields within a dispersal range. The crop change is the model's
# rule; growth and dispersal, which need the year's common random factor and the whole lattice, are a
# callback stored with the problem. The model follows the source's published R script.

"""
    rootworm_probabilities(a, b, Q, policy) -> NamedTuple

The yearly probabilities with which a field stays in or turns to maize in the
model of [Szalai2014](@cite): `pE` that a field in another crop is sown with
maize, `pF` that first-year maize stays maize, and `pC` that continuous maize
stays maize, so that a share `a` of the fields is maize and a share `b` of the
maize is continuous, with `Q = pF / pC`. With mandatory rotation after
`policy + 2` years of maize (`policy` 1, 2 or 3; 0 for none), `pC` solves the
polynomial for the steady state, as in the source's script. Returns `nothing`
if a probability falls outside `[0, 1]`.
"""
function rootworm_probabilities(a::Real, b::Real, Q::Real, policy::Integer)
    policy in 0:3 || throw(ArgumentError("policy must be 0, 1, 2 or 3, got $policy"))
    A = a * b / (Q * a * (1 - b))
    pE = a * (1 - b) / (1 - a)
    pC = if policy == 0
        1 / (Q * (1 / b - 1) + 1)
    elseif policy == 1
        (sqrt(4A + 1) - 1) / 2
    elseif policy == 2
        z = cbrt(sqrt(27A^2 + 14A + 3) / (2 * 3^1.5) + (27A + 7) / 54)
        z - 2 / (9z) - 1 / 3
    else
        z = cbrt(sqrt(256A^3 + 203A^2 + 88A + 16) / (2 * 3^1.5) + (45A + 20) / 54)
        w = sqrt((36z^2 - 15z - 48A - 8) / z)
        sqrt(15 / (2w) - z + (12A + 2) / (9z) - 5 / 6) / 2 - w / 12 - 1 / 4
    end
    pF = Q * pC
    all(x -> 0 <= x <= 1, (pE, pF, pC)) || return nothing
    return (; pE, pF, pC)
end

"""
    build_rootworm_rotation() -> CAModel

The western corn rootworm model of [Szalai2014](@cite), in yearly steps on a
torus of fields. Each cell holds its `crop`, the years it has been in maize
(0 another crop, 1 first-year maize, 2 or more continuous maize), and the
rootworm density `pop` (adults per plant). Each year, in every cell at once,
continuous maize stays maize with probability `pC` (and its count rises),
first-year maize with probability `pF`, and another crop turns to maize with
probability `pE` ([`rootworm_probabilities`](@ref)); otherwise the field turns
to another crop. With mandatory rotation (`policy` 1 to 3) a field that would
be in maize for more than `policy + 2` years turns to another crop. Growth and
dispersal are [`rootworm_step!`](@ref), a callback stored with the problem.
"""
function build_rootworm_rotation()
    rule = LocalRule(Offsets((1, 0))) do c, nbrs, ctx
        p = ctx.p
        u = uniform(ctx, 1)
        crop = c.crop > 1 ? (u < p.pC ? c.crop + Int8(1) : Int8(0)) :
               c.crop == 1 ? (u < p.pF ? Int8(2) : Int8(0)) :
               (u < p.pE ? Int8(1) : Int8(0))
        p.policy > 0 && crop > p.policy + 2 && (crop = Int8(0))
        return (; crop, pop=c.pop)
    end
    states = [StateVariable(:crop, Int8(0);
                            description="years in maize: 0 another crop, 1 first-year maize, 2+ continuous maize"),
              StateVariable(:pop, 0.0; description="rootworm adults per plant")]
    parameters = [Parameter(:a, 0.4; description="share of fields in maize"),
                  Parameter(:b, 0.5;
                            description="share of maize fields in continuous maize"),
                  Parameter(:Q, 1.0;
                            description="farmers' decisions to rotate first-year or continuous maize (pF / pC)"),
                  Parameter(:policy, 0;
                            description="mandatory rotation after policy + 2 years of maize (0 for none)"),
                  Parameter(:pE, 0.0;
                            description="probability another crop turns to maize"),
                  Parameter(:pF, 0.0;
                            description="probability first-year maize stays maize"),
                  Parameter(:pC, 0.0;
                            description="probability continuous maize stays maize"),
                  Parameter(:fieldsize, 15.0; description="area of a field (ha)"),
                  Parameter(:range, 2000.0; description="dispersal range of adults (m)"),
                  Parameter(:threshold, 1.0;
                            description="economic threshold (adults per plant)"),
                  Parameter(:initial, 0.8;
                            description="mean initial density in continuous maize (share of the threshold)"),
                  Parameter(:initial_sd, 0.3;
                            description="standard deviation of the initial density (share of its mean)"),
                  Parameter(:grA, 4.0; description="growth rate at low density"),
                  Parameter(:grB, 1.0; description="growth rate at high density"),
                  Parameter(:grA_density, 10.0;
                            description="density below which growth is grA (adults per plant)"),
                  Parameter(:grB_density, 40.0;
                            description="density above which growth is grB (adults per plant)"),
                  Parameter(:growth_sd, 0.3;
                            description="standard deviation of the growth rate (share of its mean)"),
                  Parameter(:year_effect, 0.3;
                            description="half-width of the year's uniform factor on growth (share of the mean)"),
                  Parameter(:emigration_min, 0.2;
                            description="share emigrating at low density"),
                  Parameter(:emigration_max, 0.6;
                            description="share emigrating at high density"),
                  Parameter(:emigration_min_density, 0.6;
                            description="density below which emigration is emigration_min (adults per plant)"),
                  Parameter(:emigration_max_density, 10.0;
                            description="density above which emigration is emigration_max (adults per plant)"),
                  Parameter(:emigration_sd, 0.3;
                            description="standard deviation of the emigrating share (share of its mean)")]
    return CAModel(rule; name=:rootworm_rotation, states, parameters,
                   description="Szalai et al. (2014) western corn rootworm and crop rotation")
end

"""
    rootworm_rings(p) -> Int

The dispersal range of [Szalai2014](@cite) in fields: `range` over the side of
a field, `sqrt(fieldsize × 10⁴)` m, rounded, and at least 1.
"""
rootworm_rings(p) = max(1, round(Int, p.range / sqrt(p.fieldsize * 10_000)))

# Sum over the cells at Chebyshev distance 1..k on a torus of w(d) * X, with w(d) = k - d + 1 (no self).
function _rootworm_ring_sum(X::AbstractMatrix, k::Integer)
    n1, n2 = size(X)
    out = zeros(n1, n2)
    for j in 1:n2, i in 1:n1
        s = 0.0
        for dj in (-k):k, di in (-k):k
            d = max(abs(di), abs(dj))
            d == 0 && continue
            s += (k - d + 1) * X[mod1(i + di, n1), mod1(j + dj, n2)]
        end
        out[i, j] = s
    end
    return out
end

"""
    rootworm_step!(u, p, rng; grow=true) -> u

The rootworm's year in the model of [Szalai2014](@cite), in place, after the
crops have changed. If `grow`, populations outside continuous maize are lost,
and those in continuous maize grow by a rate falling linearly from `grA` at
`grA_density` to `grB` at `grB_density`, times the year's factor (uniform
within `year_effect` of 1, the same for every field) and a normal factor of
mean 1 and standard deviation `growth_sd` for each field. Then a share of the
adults of each continuous maize field with maize within the range, rising
linearly from `emigration_min` at `emigration_min_density` to
`emigration_max` at `emigration_max_density` and times a normal factor of
mean 1 and standard deviation `emigration_sd` (kept within the population),
leaves, and lands in the maize fields within [`rootworm_rings`](@ref) fields
(Chebyshev distance on the torus), in proportion to `k - d + 1` at distance
`d`. As in the source's script, densities below zero are set to zero.
"""
function rootworm_step!(u, p, rng; grow::Bool=true)
    crop, pop = u.crop, u.pop
    if grow
        year = p.year_effect == 0 ? 1.0 : 1 - p.year_effect + 2p.year_effect * rand(rng)
        slope = (p.grB - p.grA) / (p.grB_density - p.grA_density)
        for I in eachindex(pop)
            crop[I] > 1 || (pop[I]=0.0; continue)
            x = pop[I]
            x == 0 && continue
            rate = min(p.grA, max(p.grB, slope * (x - p.grA_density) + p.grA))
            pop[I] = max(x * rate * year * (1 + p.growth_sd * randn(rng)), 0.0)
        end
    end
    k = rootworm_rings(p)
    maize = Float64.(crop .> 0)
    weight = _rootworm_ring_sum(maize, k)
    slope = (p.emigration_max - p.emigration_min) /
            (p.emigration_max_density - p.emigration_min_density)
    out = zeros(size(pop))
    share = zeros(size(pop))
    for I in CartesianIndices(pop)
        (crop[I] > 1 && weight[I] > 0) || continue
        e = min(p.emigration_max,
                max(p.emigration_min,
                    slope * (pop[I] - p.emigration_min_density) +
                    p.emigration_min)) * pop[I]
        e = clamp(e * (1 + p.emigration_sd * randn(rng)), 0.0, pop[I])
        out[I] = e
        share[I] = e / weight[I]
    end
    arrive = _rootworm_ring_sum(share, k)
    for I in eachindex(pop)
        pop[I] += (crop[I] > 0 ? arrive[I] : 0.0) - out[I]
    end
    return u
end

"""
    rootworm_callback() -> DiscreteCallback

A callback that runs [`rootworm_step!`](@ref) after each year's crop change,
with a random number generator seeded from the problem's seed and the year.
"""
function rootworm_callback()
    affect! = function (integ)
        rng = Random.Xoshiro(hash((integ.prob.seed, integ.t, :rootworm_rotation)))
        rootworm_step!(integ.u, integ.p, rng)
        return nothing
    end
    return DiscreteCallback((u, t, integ) -> true, affect!; save_positions=(false, false))
end

"""
    problem_rootworm_rotation(; n=100, a=0.4, rotated=0.5, Q=1, policy=0, years=10, seed=0,
                              pmap=(;), kwargs...) -> CAProblem

The setup of [Szalai2014](@cite): a torus of `n × n` fields with a share `a`
in maize, of which a share `rotated` is first-year (rotated) maize, crops drawn
at random from their long-run shares, and initial rootworm densities in
continuous maize drawn from a normal distribution of mean `initial ×
threshold` and standard deviation `initial_sd` times that (negative draws set
to zero), all from `seed`. The state at time 0 is the first year after
dispersal; the run lasts `years` years in all. The other parameters default to
the source's most likely values (its Table 1). Throws an `ArgumentError` for a
combination the crop dynamics cannot sustain.
"""
function problem_rootworm_rotation(; n::Integer=100, a::Real=0.4, rotated::Real=0.5,
                                   Q::Real=1, policy::Integer=0, years::Integer=10,
                                   seed::Integer=0, pmap=(;), kwargs...)
    0 < a < 1 || throw(ArgumentError("a must lie in (0, 1), got $a"))
    0 < rotated < 1 || throw(ArgumentError("rotated must lie in (0, 1), got $rotated"))
    years >= 1 || throw(ArgumentError("years must be at least 1, got $years"))
    b = 1 - rotated
    probs = rootworm_probabilities(a, b, Q, policy)
    probs === nothing &&
        throw(ArgumentError("a = $a, rotated = $rotated, Q = $Q and policy $policy give crop probabilities outside [0, 1]"))
    model = build_rootworm_rotation()
    given = merge((; a=Float64(a), b=Float64(b), Q=Float64(Q), policy), probs,
                  (; (Symbol(k) => v for (k, v) in pairs(pmap))...))
    p = merge(defaults(model), given)
    rng = Random.Xoshiro(hash((seed, 0, :rootworm_rotation)))
    pE, pF, pC = probs.pE, probs.pF, probs.pC
    shares = [1 - a, a * (1 - b), a * (1 - b) * pF, a * (1 - b) * pF * pC,
              a * (1 - b) * pF * pC^2, a * (1 - b) * pF * pC^3][1:(policy + 2)]
    push!(shares, 1 - sum(shares))
    cum = cumsum(shares)
    crop = [Int8(min(searchsortedfirst(cum, rand(rng)), length(cum)) - 1)
            for _ in 1:n, _ in 1:n]
    m = p.initial * p.threshold
    pop = [c > 1 ? max(m + p.initial_sd * m * randn(rng), 0.0) : 0.0 for c in crop]
    u0 = (; crop, pop)
    rootworm_step!(u0, p, rng; grow=false)
    return CAProblem(model, u0, years - 1, given; space=Grid(n, n; boundary=Periodic()),
                     seed,
                     callback=rootworm_callback(), kwargs...)
end

"""
    rootworm_risk(u; threshold=1) -> NamedTuple

The shares of maize fields (`all`), of first-year maize (`first_year`) and of
continuous maize (`continuous`) whose rootworm density exceeds `threshold`
(adults per plant) in the state `u`, the source's yearly measure (`NaN` where
there are no such fields).
"""
function rootworm_risk(u; threshold::Real=1)
    over = u.pop .> threshold
    share(mask) = (k=count(mask); k == 0 ? NaN : count(over .& mask) / k)
    return (; all=share(u.crop .> 0), first_year=share(u.crop .== 1),
            continuous=share(u.crop .> 1))
end

"""
    rootworm_output(sol; threshold=1) -> Float64

The source's output: the mean over the years after the first (the saved
states at times 1 onwards) of the share of maize fields above `threshold`
([`rootworm_risk`](@ref)).
"""
function rootworm_output(sol; threshold::Real=1)
    length(sol.u) > 1 || throw(ArgumentError("the solution needs more than the first year"))
    return sum(rootworm_risk(u; threshold).all for u in sol.u[2:end]) / (length(sol.u) - 1)
end

_register!(:rootworm_rotation; build=build_rootworm_rotation,
           problem=problem_rootworm_rotation)
