# Rebaudo, Crespo-Pérez, Silvain & Dangles (2011): the Andean potato tuber moth Tecia solanivora spreading
# through a valley of Ecuador, by flight between neighbouring cells and in sacks of seed potatoes carried
# by farmers between villages, as in the authors' published NetLogo model. One step is a moth generation.
# The moths' generation is the rule; the farmers, who act one after another, are a callback stored with
# the problem. A farmer is a bit: the cells' `farmers` layer marks who stands there, and `carriers` those
# of them carrying moths, so the farmers live in the state and the model runs on any backend.

const PTM_NONE = Int8(0)
const PTM_BIOLOGICAL = Int8(1)
const PTM_RED = Int8(2)
const PTM_BLUE = Int8(3)
const PTM_MAX_FARMERS = 64

# Sharpe–DeMichele survival of a stage at `T` kelvin, with coefficients (a, b, c, d, e, f).
function _ptm_stage_survival(T, (a, b, c, d, e, f), R)
    return a * T / 298.16 * exp(b * (1 / 298.16 - 1 / T)) /
           (1 + exp(c / R * (1 / d - 1 / T)) + exp(e / R * (1 / f - 1 / T)))
end

"""
    potato_moth_survival(T, p) -> Float64

The share of immature *Tecia solanivora* surviving a generation at `T` °C
because of temperature: the product of the Sharpe–DeMichele survival of eggs,
larvae and pupae, with the coefficients of [Rebaudo2011](@cite) (`egg`,
`larva`, `pupa` in `p`) and the gas constant `p.R`.
"""
function potato_moth_survival(T, p)
    K = T + 273.15
    return _ptm_stage_survival(K, p.egg, p.R) * _ptm_stage_survival(K, p.larva, p.R) *
           _ptm_stage_survival(K, p.pupa, p.R)
end

# Immatures after a generation's mortality, as the source's code subtracts them.
function _ptm_survivors(immature, T, p)
    crude = 1 - exp(-(p.mort_dispersal + p.mort_predation) * p.mort_time)
    a = immature - immature * crude
    return a - a * (1 - potato_moth_survival(T, p))
end

# Gravid females from the survivors: half are females, all mated, then farmers kill `control`.
function _ptm_gravid(survivors, p)
    g = survivors / 2
    return g - g * p.control
end

"""
    build_potato_moth_farmers() -> CAModel

The potato tuber moth model of [Rebaudo2011](@cite), as in the authors'
published NetLogo code, in steps of one moth generation on a lattice of
500 m cells. Each cell holds immature moths (`immature`, a real number), the
origin of its moths (`source`: 0 none, 1 spread by flight from the first
infestation, 2 and 3 brought by a farmer of the red or the blue team), the
temperature the coming generation meets (`temperature_now`), and fixed layers
of mean temperature (`temperature`, °C) and village population
(`population`). The farmers stand on cells as bits of `farmers`, and those
carrying moths as bits of `carriers` (see [`potato_moth_farmers_step!`](@ref)).
Each generation, in every cell at once:

1. immatures die of crude mortality `1 - exp(-(mort_dispersal +
   mort_predation) mort_time)` and then of temperature
   ([`potato_moth_survival`](@ref) at `temperature_now`); half the survivors
   are gravid females, of which farmers kill `control`;
2. an empty cell (no surviving immatures) takes the `source` of a neighbour
   with gravid females, chosen at random;
3. a share `diffusion` of the gravid females leaves each cell, a quarter to
   each of its four neighbours; shares that would leave the lattice stay;
4. gravid females, capped at `capacity`, give immatures by the logistic map
   `capacity G / (G + (capacity - G) exp(-growth))`;
5. the temperature the next generation meets is drawn from a normal
   distribution about `temperature` with standard deviation `temperature_sd`.

With `neighbourhood = false`, moths do not leave their cells (2 and 3 are
skipped).
"""
function build_potato_moth_farmers()
    nbhd = Offsets([CartesianIndex(-1, 0), CartesianIndex(1, 0), CartesianIndex(0, -1),
                    CartesianIndex(0, 1)])
    rule = LocalRule(nbhd) do c, nbrs, ctx
        p = ctx.p
        i, j = Tuple(ctx.index)
        n1, n2 = p.size
        inside = (i > 1, i < n1, j > 1, j < n2)
        survivors = _ptm_survivors(c.immature, c.temperature_now, p)
        G = _ptm_gravid(survivors, p)
        source = c.source
        if p.neighbourhood
            gs = ntuple(k -> inside[k] ?
                             _ptm_gravid(_ptm_survivors(nbrs[k].immature,
                                                        nbrs[k].temperature_now, p), p) :
                             0.0, Val(4))
            if survivors == 0
                m = count(>(0), gs)
                if m > 0
                    pick = min(floor(Int, uniform(ctx, 1) * m) + 1, m)
                    for k in 1:4
                        gs[k] > 0 && (pick -= 1)
                        if pick == 0
                            source = nbrs[k].source
                            break
                        end
                    end
                end
            end
            kept = (4 - count(inside)) * G
            G = (1 - p.diffusion) * G + p.diffusion / 4 * (sum(gs) + kept)
        end
        K = p.capacity
        G = min(G, K)
        immature = K * G / (G + (K - G) * exp(-p.growth))
        z = sqrt(-2 * log(1 - uniform(ctx, 2))) * cospi(2 * uniform(ctx, 3))
        return (immature, source, temperature_now=c.temperature + p.temperature_sd * z,
                temperature=c.temperature, population=c.population, farmers=c.farmers,
                carriers=c.carriers)
    end
    states = [StateVariable(:immature, 0.0; description="immature moths in the cell"),
              StateVariable(:source, PTM_NONE; domain=Int8(0):Int8(3),
                            description="origin of the cell's moths: 0 none, 1 flight, 2 red team, 3 blue team"),
              StateVariable(:temperature_now, 0.0;
                            description="temperature the coming generation meets (°C)"),
              StateVariable(:temperature, 0.0; description="mean temperature (°C)"),
              StateVariable(:population, Int16(0);
                            description="population of the village in the cell, 0 if none"),
              StateVariable(:farmers, UInt64(0);
                            description="the farmers standing in the cell, as bits"),
              StateVariable(:carriers, UInt64(0);
                            description="the farmers in the cell carrying moths, as bits")]
    parameters = [Parameter(:mort_dispersal, 0.06;
                            description="crude mortality rate from dispersal"),
                  Parameter(:mort_predation, 0.145;
                            description="crude mortality rate from predation"),
                  Parameter(:mort_time, 2.0;
                            description="time over which crude mortality acts"),
                  Parameter(:egg, (0.822, -758.5, -212100.0, 281.9, 405200.0, 303.8);
                            description="Sharpe–DeMichele coefficients of egg survival"),
                  Parameter(:larva, (0.758, -180.2, -475700.0, 282.7, 1298000.0, 301.5);
                            description="Sharpe–DeMichele coefficients of larval survival"),
                  Parameter(:pupa, (0.9, -73.72, -1263000.0, 286.5, 1095000.0, 306.3);
                            description="Sharpe–DeMichele coefficients of pupal survival"),
                  Parameter(:R, 1.987; description="gas constant (cal/(mol K))"),
                  Parameter(:temperature_sd, 0.1;
                            description="standard deviation of a generation's temperature (°C)"),
                  Parameter(:control, 0.1;
                            description="share of gravid females farmers kill each generation"),
                  Parameter(:neighbourhood, true;
                            description="whether moths fly to neighbouring cells"),
                  Parameter(:diffusion, 0.3;
                            description="share of gravid females leaving a cell each generation"),
                  Parameter(:capacity, 1000.0; description="most gravid females in a cell"),
                  Parameter(:growth, 5.0;
                            description="logistic growth rate per generation"),
                  Parameter(:size, (80, 80); description="size of the lattice"),
                  Parameter(:farmers, 10; description="number of farmers"),
                  Parameter(:moves, 6; description="journeys of each farmer a generation"),
                  Parameter(:p_nearest, 0.8;
                            description="probability a farmer's journey follows the source's choice rule rather than a random village"),
                  Parameter(:population_scale, 3000.0;
                            description="population dividing village populations in the choice rule"),
                  Parameter(:inoculum, 90.0;
                            description="immatures a carrying farmer leaves in an uninfested cell"),
                  Parameter(:knowledge, ();
                            description="each farmer's pest control knowledge (0-100)"),
                  Parameter(:infection, :threshold;
                            description=":threshold (knowledge below threshold) or :probability (1 - knowledge/100)"),
                  Parameter(:threshold, 100.0;
                            description="knowledge below which a farmer buys infested potatoes")]
    return CAModel(rule; name=:potato_moth_farmers, states, parameters,
                   description="Rebaudo et al. (2011) potato tuber moth spread by flight and by farmers")
end

# The team of farmer `f` of `n`: the first half red, the second blue.
_ptm_team(f, n) = f <= n ÷ 2 ? PTM_RED : PTM_BLUE

function _ptm_catches(p, f, rng)
    κ = p.knowledge[f]
    return p.infection === :threshold ? κ < p.threshold : rand(rng) < 1 - κ / 100
end

"""
    potato_moth_farmers_step!(u, p, rng) -> u

The farmers' part of a generation, applied to the state `u` in place after
the moths' step: `moves` times, every farmer, in a random order, first

- if it carries no moths and its cell has immatures, makes every farmer there
  whose knowledge is below `threshold` (`infection = :threshold`, the
  source's code) or, with `infection = :probability`, each with probability
  `1 - knowledge/100`, a carrier;
- if it carries moths and its cell has none, leaves `inoculum` immatures there,
  marked with its team;

and then every farmer, in a random order, travels to another village (a cell
with a population). With probability `p_nearest` it follows the source's
rule: it goes to the village `j` minimising `w n pop_j d_j / population_scale`,
where `w` is its cell's population over `population_scale`, `n` its
immatures over `capacity` and `d_j` the distance, ties broken at random (so a
farmer standing where there are no moths, or outside a village, picks a
village at random). Otherwise it picks one at random. Carriers stay carriers.
"""
function potato_moth_farmers_step!(u, p, rng)
    n = p.farmers
    n == 0 && return u
    pos = fill(CartesianIndex(0, 0), n)
    carrying = falses(n)
    for I in CartesianIndices(u.farmers)
        bits = u.farmers[I]
        bits == 0 && continue
        for f in 1:n
            if (bits >> (f - 1)) & 1 == 1
                pos[f] = I
                carrying[f] = (u.carriers[I] >> (f - 1)) & 1 == 1
            end
        end
    end
    missing = findfirst(==(CartesianIndex(0, 0)), pos)
    missing === nothing ||
        throw(ArgumentError("farmer $missing of $n stands on no cell of the state"))
    villages = [I for I in CartesianIndices(u.population) if u.population[I] > 0]
    for _ in 1:p.moves
        for f in Random.shuffle(rng, 1:n)
            I = pos[f]
            if !carrying[f]
                if u.immature[I] > 0
                    for g in 1:n
                        pos[g] == I && !carrying[g] && _ptm_catches(p, g, rng) &&
                            (carrying[g] = true)
                    end
                end
            elseif u.immature[I] == 0
                u.immature[I] = p.inoculum
                u.source[I] = _ptm_team(f, n)
            end
        end
        for f in Random.shuffle(rng, 1:n)
            I = pos[f]
            options = filter(!=(I), villages)
            if rand(rng, 0:99) < 100 * p.p_nearest
                w = u.population[I] / p.population_scale * u.immature[I] / p.capacity
                score(J) = w * (u.population[J] / p.population_scale *
                                hypot((Tuple(J) .- Tuple(I))...))
                best = minimum(score, options)
                pos[f] = rand(rng, filter(J -> score(J) == best, options))
            else
                pos[f] = rand(rng, options)
            end
        end
    end
    fill!(u.farmers, 0)
    fill!(u.carriers, 0)
    for f in 1:n
        bit = UInt64(1) << (f - 1)
        u.farmers[pos[f]] |= bit
        carrying[f] && (u.carriers[pos[f]] |= bit)
    end
    return u
end

"""
    potato_moth_farmers_callback() -> DiscreteCallback

A callback that runs [`potato_moth_farmers_step!`](@ref) after each
generation's moth step, with a random number generator seeded from the
problem's seed and the generation.
"""
function potato_moth_farmers_callback()
    affect! = function (integ)
        rng = Random.Xoshiro(hash((integ.prob.seed, Int(integ.t), :potato_moth)))
        potato_moth_farmers_step!(integ.u, integ.p, rng)
        return nothing
    end
    return DiscreteCallback((u, t, integ) -> true, affect!; save_positions=(false, false))
end

"""
    potato_moth_landscape() -> NamedTuple

The landscape of the source's published model of the Simiatug valley, Ecuador:
80 × 80 cells of 500 m, indexed `[x + 1, y + 1]` from the source's patch
coordinates, with `elevation` (m), mean `temperature` (°C, whole degrees),
`village` (an identifier, 0 where there is none) and `population` (500 for a
village, 3000 for Simiatug; the source's code calls these populations
invented). The data come with the source's code (Rebaudo 2013, Academic Free
License 3.0).
"""
function potato_moth_landscape()
    path = joinpath(MODELS_DIR, "potato_moth_farmers", "data", "simiatug.csv")
    elevation, temperature = zeros(Int, 80, 80), zeros(Int, 80, 80)
    village, population = zeros(Int, 80, 80), zeros(Int, 80, 80)
    for line in Iterators.drop(eachline(path), 1)
        x, y, e, t, v, pop = parse.(Int, split(line, ','))
        elevation[x + 1, y + 1], temperature[x + 1, y + 1] = e, t
        village[x + 1, y + 1], population[x + 1, y + 1] = v, pop
    end
    return (; elevation, temperature, village, population)
end

"""
    problem_potato_moth_farmers(; farmers=10, moves=6, knowledge=:uniform,
                                infection=:threshold, generations=100, seed=0,
                                pmap=(;), kwargs...) -> CAProblem

The setup of the published model of [Rebaudo2011](@cite): the Simiatug valley
([`potato_moth_landscape`](@ref)), with 90 immature moths in the cell of
Simiatug village, and `farmers` farmers (an even number, half in each team;
0 for none) placed at random, none carrying moths, for `generations`
generations. Each farmer's pest control knowledge is drawn from `seed`:

- `knowledge = :uniform`, as in the source's code: a whole number from 0 to 99;
- `knowledge = (mean, sd)`: a normal draw, kept within 0 and 100, as in the
  source's experiments;
- or a vector giving each farmer's knowledge.

`infection` is `:threshold` (the source's code: a farmer buys infested
potatoes if its knowledge is below `threshold`, 100 by default, so every
farmer does) or `:probability` (with probability `1 - knowledge/100`). The
temperature the first generation meets is drawn from `seed`. Nothing lies
beyond the valley.
"""
function problem_potato_moth_farmers(; farmers::Integer=10, moves::Integer=6,
                                     knowledge=:uniform, infection::Symbol=:threshold,
                                     generations::Integer=100, seed::Integer=0, pmap=(;),
                                     kwargs...)
    (0 <= farmers <= PTM_MAX_FARMERS && iseven(farmers)) ||
        throw(ArgumentError("farmers must be an even number from 0 to $PTM_MAX_FARMERS, got $farmers"))
    moves >= 0 || throw(ArgumentError("moves must not be negative, got $moves"))
    infection in (:threshold, :probability) ||
        throw(ArgumentError("infection must be :threshold or :probability, got :$infection"))
    rng = Random.Xoshiro(hash((seed, 0, :potato_moth)))
    κ = if knowledge === :uniform
        Tuple(Float64.(rand(rng, 0:99, farmers)))
    elseif knowledge isa Tuple{Real,Real}
        μ, σ = knowledge
        Tuple(clamp.(μ .+ σ .* randn(rng, farmers), 0.0, 100.0))
    elseif knowledge isa AbstractVector{<:Real} && length(knowledge) == farmers
        Tuple(Float64.(knowledge))
    else
        throw(ArgumentError("knowledge must be :uniform, (mean, sd) or a vector of $farmers values"))
    end
    land = potato_moth_landscape()
    dims = size(land.temperature)
    model = build_potato_moth_farmers()
    given = merge((; farmers, moves, knowledge=κ, infection, size=dims),
                  (; (Symbol(k) => v for (k, v) in pairs(pmap))...))
    p = merge(defaults(model), given)
    temperature = Float64.(land.temperature)
    immature = zeros(dims)
    source = zeros(Int8, dims)
    simiatug = CartesianIndex(50, 41)
    immature[simiatug], source[simiatug] = p.inoculum, PTM_BIOLOGICAL
    temperature_now = temperature .+ p.temperature_sd .* randn(rng, dims)
    bits = zeros(UInt64, dims)
    for f in 1:farmers
        bits[rand(rng, CartesianIndices(dims))] |= UInt64(1) << (f - 1)
    end
    u0 = (; immature, source, temperature_now, temperature,
          population=Int16.(land.population), farmers=bits)
    return CAProblem(model, u0, generations, given; space=Grid(dims...), seed,
                     callback=potato_moth_farmers_callback(), kwargs...)
end

"""
    potato_moth_counts(u; threshold=0) -> NamedTuple

The cells with more than `threshold` immature moths (`infested`; the source
does not say what makes a cell infested, and its code counts any moths), those whose moths came by flight
from the first infestation (`biological`) or were brought by the `red` or the
`blue` team (the source's colours, which a cell keeps once given), and the
farmers carrying moths (`carriers`).
"""
function potato_moth_counts(u; threshold::Real=0)
    return (infested=count(>(threshold), u.immature),
            biological=count(==(PTM_BIOLOGICAL), u.source),
            red=count(==(PTM_RED), u.source), blue=count(==(PTM_BLUE), u.source),
            carriers=sum(count_ones, u.carriers))
end

"""
    potato_moth_invasion_time(prob; fraction=0.05, threshold=0) -> Union{Int,Nothing}

The source's measure of invasion speed: the first generation at which at
least `fraction` of the cells have more than `threshold` immature moths,
running `prob` until then, or `nothing` if that does not happen within its
time span.
"""
function potato_moth_invasion_time(prob; fraction::Real=0.05, threshold::Real=0)
    0 < fraction <= 1 || throw(ArgumentError("fraction must lie in (0, 1], got $fraction"))
    integ = init(prob; save_everystep=false, save_start=false)
    target = fraction * length(integ.u.immature)
    while potato_moth_counts(integ.u; threshold).infested < target
        integ.t >= prob.tspan[2] && return nothing
        step!(integ)
    end
    return Int(integ.t)
end

_register!(:potato_moth_farmers; build=build_potato_moth_farmers,
           problem=problem_potato_moth_farmers)
