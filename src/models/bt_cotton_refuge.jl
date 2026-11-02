# Malaquias et al. (2021): competition between Spodoptera frugiperda and Helicoverpa armigera in the
# non-Bt refuge of Bt cotton. A cell is a plant, empty or holding one of the two species; a step is a
# generation. Occupied cells lose their larvae to natural mortality, to insecticide (in sprayed
# generations) and to competition, whose survival depends on local density and on meeting the other
# species; empty cells are colonised from their Moore neighbours in proportion to fecundity. The
# spray decision is global: a step callback stored with the problem sets `spray` whenever either
# species occupies a control level of all cells.

const BT_EMPTY = Int8(0)
const BT_SF = Int8(1)
const BT_HA = Int8(2)

"""
    bt_cotton_survival(source=:figure1) -> NamedTuple

The larval survival probabilities of the refuge model of [Malaquias2021](@cite),
as parameters: `sL_*` at low and `sH_*` at high local density; `*_Sf` and `*_Ha`
for one species alone, `*_Sf_Ha` for *S. frugiperda* meeting *H. armigera* and
`*_Ha_Sf` the reverse. `:figure1` gives the values of the source's Figure 1 and
rule text (the default); `:table1` those of its Table 1, which contradict them.
"""
function bt_cotton_survival(source::Symbol=:figure1)
    source === :figure1 &&
        return (sL_Sf=0.75, sL_Ha=0.80, sL_Sf_Ha=0.60, sL_Ha_Sf=0.49, sH_Sf=0.48,
                sH_Ha=0.22, sH_Sf_Ha=0.18, sH_Ha_Sf=0.47)
    source === :table1 &&
        return (sL_Sf=0.63, sL_Ha=0.46, sL_Sf_Ha=0.58, sL_Ha_Sf=0.90, sH_Sf=0.31,
                sH_Ha=0.38, sH_Sf_Ha=0.11, sH_Ha_Sf=0.16)
    return throw(ArgumentError("source must be :figure1 or :table1, got :$source"))
end

"""
    build_bt_cotton_refuge() -> CAModel

The refuge competition model of [Malaquias2021](@cite). Each cell is a non-Bt
cotton plant, empty (0) or holding *S. frugiperda* (1) or *H. armigera* (2); a
step is a generation, and cells update synchronously from their Moore
neighbours.

An occupied cell empties with probability `mortality` (natural enemies and
dispersal); then, in a sprayed generation (`spray`), with probability `kill`,
unless its species is resistant (`resistant_Sf`, `resistant_Ha`), in which case
the larva escapes with probability `escape`; then it survives competition with
probability `s`. The survival table ([`bt_cotton_survival`](@ref)) is the
high-density one if at least `threshold * 8` neighbours are occupied, and `s`
is the value for meeting the other species with probability `n_other / 8`, the
share of neighbours holding it, and for the species alone otherwise.

An empty cell with `n₁` and `n₂` neighbours of the two species is colonised by
*S. frugiperda* with probability `f₁ = emergence n₁ R₁ / (n₁ R₁ + n₂ R₂)`, and
failing that by *H. armigera* with probability `f₂`, defined alike, where `R₁`
and `R₂` are the fecundities `fecundity_Sf` and `fecundity_Ha`.

`spray` is set by the callback [`bt_cotton_spray_callback`](@ref), which the
problem constructor stores with the problem.
"""
function build_bt_cotton_refuge()
    rule = LocalRule(Moore(1)) do cell, nbrs, ctx
        p = ctx.p
        n1 = count(==(BT_SF), nbrs)
        n2 = count(==(BT_HA), nbrs)
        if cell == BT_EMPTY
            w1, w2 = n1 * p.fecundity_Sf, n2 * p.fecundity_Ha
            w1 + w2 == 0 && return BT_EMPTY
            bernoulli(ctx, p.emergence * w1 / (w1 + w2), 5) && return BT_SF
            bernoulli(ctx, p.emergence * w2 / (w1 + w2), 6) && return BT_HA
            return BT_EMPTY
        end
        sf = cell == BT_SF
        bernoulli(ctx, p.mortality, 1) && return BT_EMPTY
        if p.spray
            resistant = sf ? p.resistant_Sf : p.resistant_Ha
            pkill = resistant ? p.kill * (1 - p.escape) : p.kill
            bernoulli(ctx, pkill, 2) && return BT_EMPTY
        end
        high = n1 + n2 >= p.threshold * 8
        meets = bernoulli(ctx, (sf ? n2 : n1) / 8, 3)
        s = if sf
            high ? (meets ? p.sH_Sf_Ha : p.sH_Sf) : (meets ? p.sL_Sf_Ha : p.sL_Sf)
        else
            high ? (meets ? p.sH_Ha_Sf : p.sH_Ha) : (meets ? p.sL_Ha_Sf : p.sL_Ha)
        end
        return bernoulli(ctx, s, 4) ? cell : BT_EMPTY
    end
    parameters = [Parameter(:mortality, 0.25;
                            description="natural mortality and dispersal"),
                  Parameter(:kill, 0.99; description="insecticide mortality"),
                  Parameter(:spray, false;
                            description="this generation is sprayed (set by the callback)"),
                  Parameter(:resistant_Sf, false;
                            description="S. frugiperda resists the insecticide"),
                  Parameter(:resistant_Ha, false;
                            description="H. armigera resists the insecticide"),
                  Parameter(:escape, 1.0;
                            description="probability a resistant larva escapes the spray"),
                  Parameter(:threshold, 0.5;
                            description="occupied share of neighbours for high density"),
                  Parameter(:fecundity_Sf, 400.0;
                            description="reproductive capacity of S. frugiperda"),
                  Parameter(:fecundity_Ha, 725.06;
                            description="reproductive capacity of H. armigera"),
                  Parameter(:emergence, 0.86;
                            description="adult emergence, scaling colonisation")]
    for (k, v) in pairs(bt_cotton_survival())
        push!(parameters, Parameter(k, v; description="larval survival"))
    end
    return CAModel(rule; name=:bt_cotton_refuge,
                   states=[StateVariable(:plant, BT_EMPTY; domain=Int8(0):Int8(2),
                                         description="0 empty, 1 S. frugiperda, 2 H. armigera")],
                   parameters,
                   description="Malaquias et al. (2021) S. frugiperda and H. armigera in a Bt cotton refuge")
end

"""
    bt_cotton_spray_callback(; control_level=0.2, insecticide=true) -> DiscreteCallback

A callback that, after every generation, sprays the next one (sets the
parameter `spray`) if either species occupies at least `control_level` of all
cells, as in [Malaquias2021](@cite), and stops spraying otherwise. With
`insecticide = false` it never sprays.
"""
function bt_cotton_spray_callback(; control_level::Real=0.2, insecticide::Bool=true)
    affect! = function (integ)
        u = integ.u
        spray = insecticide &&
                max(count(==(BT_SF), u), count(==(BT_HA), u)) >= control_level * length(u)
        integ.p = merge(integ.p, (; spray))
        return nothing
    end
    return DiscreteCallback((u, t, integ) -> true, affect!; save_positions=(false, false))
end

"""
    problem_bt_cotton_refuge(; n=1000, scenario=:no_insecticide, initial=500,
                             control_level=0.2, tspan=80, seed=0, pmap=(;), kwargs...)
        -> CAProblem

The setup of [Malaquias2021](@cite): an `n × n` torus with `initial`
individuals of each species on distinct random cells (drawn from `seed`), for
80 generations, in one of the source's five scenarios: `:no_insecticide`,
`:insecticide` (neither species resistant), `:Sf_resistant`, `:Ha_resistant`
or `:both_resistant`. The spray callback is stored with the problem; `pmap`
overrides any parameter, including the resistance of a scenario.
"""
function problem_bt_cotton_refuge(; n::Integer=1000, scenario::Symbol=:no_insecticide,
                                  initial::Integer=500, control_level::Real=0.2, tspan=80,
                                  seed::Integer=0, pmap=(;), kwargs...)
    scenarios = (:no_insecticide, :insecticide, :Sf_resistant, :Ha_resistant,
                 :both_resistant)
    scenario in scenarios ||
        throw(ArgumentError("scenario must be one of $scenarios, got :$scenario"))
    2initial <= n^2 ||
        throw(ArgumentError("$(2initial) individuals do not fit on $n × $n cells"))
    resistance = (resistant_Sf=scenario in (:Sf_resistant, :both_resistant),
                  resistant_Ha=scenario in (:Ha_resistant, :both_resistant))
    plant = fill(BT_EMPTY, n, n)
    cells = Random.randperm(Random.Xoshiro(seed), n^2)
    plant[cells[1:initial]] .= BT_SF
    plant[cells[(initial + 1):(2initial)]] .= BT_HA
    callback = bt_cotton_spray_callback(; control_level,
                                        insecticide=scenario !== :no_insecticide)
    return CAProblem(build_bt_cotton_refuge(), (; plant), tspan,
                     merge(resistance, (; (Symbol(k) => v for (k, v) in pairs(pmap))...));
                     space=Torus(n, n), seed, callback, kwargs...)
end

"""
    bt_cotton_dominance(u) -> NamedTuple

The relative dominance of each species, its share of the occupied cells:
`(Sf=..., Ha=...)`, both `NaN` when no cell is occupied.
"""
function bt_cotton_dominance(u)
    a, b = count(==(BT_SF), u), count(==(BT_HA), u)
    return a + b == 0 ? (Sf=NaN, Ha=NaN) : (Sf=a / (a + b), Ha=b / (a + b))
end

_register!(:bt_cotton_refuge; build=build_bt_cotton_refuge,
           problem=problem_bt_cotton_refuge)
