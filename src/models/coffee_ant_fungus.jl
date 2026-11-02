# Jackson, Vandermeer & Perfecto (2009): nests of the ant Azteca instabilis on the shade trees of a
# coffee farm, whose scale insects can be infected by the fungus Lecanicillium lecanii, which kills
# the nest. A cell is a shade tree: empty, a nest whose scale insects are free of the fungus, or a
# nest whose scale insects are infected. Built on the ant automaton of Vandermeer, Perfecto & Philpott
# (2008), whose density-dependent (phorid) mortality is available as `mortality = :phorid`.

const ANT_EMPTY = Int8(0)
const ANT_NEST = Int8(1)
const ANT_INFECTED = Int8(2)

"""
    build_coffee_ant_fungus(; mortality=:fungus, phases=:two) -> CAModel

The ant and fungus automaton of [Jackson2009](@cite). Each cell is a shade
tree: empty (0), an *Azteca* nest whose scale insects are free of the fungus
(1), or a nest whose scale insects are infected (2). A step is six months.
With `N` occupied and `F` infected Moore neighbours:

- an empty tree with `N > 0` gains a satellite nest, uninfected, with
  probability `s0 + s1 N`;
- an uninfected nest is infected with probability `min(1, t0 + t1 F)` if
  `F > 0`, and with the background probability `a` otherwise;
- an infected nest dies with probability `v` (virulence); uninfected nests do
  not die.

With `mortality = :phorid` the model is the fungus-free automaton of
[Vandermeer2008](@cite) that it extends: a nest dies with probability
`m0 + m1 N` (attack by phorid flies), and nests are never infected.

The source's iteration draws colonisation for every empty tree and mortality
for every nest. With `phases = :two` (the default) an iteration is two steps:
the even step colonises, and the odd step infects and kills, counting the
nests just founded among the neighbours; [`ant_iterations`](@ref) converts
steps to iterations. With `phases = :one` both happen in one step, from the
same state. Within a step the source updates "cell by cell": solve with
`Synchronous()` (the default) or `RandomSequential()` for an in-place sweep.
"""
function build_coffee_ant_fungus(; mortality::Symbol=:fungus, phases::Symbol=:two)
    mortality in (:fungus, :phorid) ||
        throw(ArgumentError("mortality must be :fungus or :phorid, got :$mortality"))
    phases in (:one, :two) ||
        throw(ArgumentError("phases must be :one or :two, got :$phases"))
    fungus = mortality === :fungus
    two = phases === :two
    rule = LocalRule(Moore(1)) do cell, nbrs, ctx
        p = ctx.p
        N = count(!=(ANT_EMPTY), nbrs)
        colonise = !two || iseven(ctx.t)
        perish = !two || isodd(ctx.t)
        if cell == ANT_EMPTY
            return colonise && N > 0 && bernoulli(ctx, p.s0 + p.s1 * N, 1) ? ANT_NEST :
                   ANT_EMPTY
        elseif !perish
            return cell
        elseif !fungus
            return bernoulli(ctx, min(1.0, p.m0 + p.m1 * N), 3) ? ANT_EMPTY : cell
        elseif cell == ANT_NEST
            F = count(==(ANT_INFECTED), nbrs)
            pt = F > 0 ? min(1.0, p.t0 + p.t1 * F) : p.a
            return bernoulli(ctx, pt, 2) ? ANT_INFECTED : ANT_NEST
        else
            return bernoulli(ctx, p.v, 3) ? ANT_EMPTY : ANT_INFECTED
        end
    end
    return CAModel(rule; name=:coffee_ant_fungus,
                   states=[StateVariable(:tree, ANT_EMPTY; domain=Int8(0):Int8(2),
                                         description="0 empty, 1 nest, 2 nest with infected scale insects")],
                   parameters=[Parameter(:s0, 0.0035;
                                         description="satellite expansion intercept"),
                               Parameter(:s1, 0.035;
                                         description="satellite expansion per occupied neighbour"),
                               Parameter(:v, 0.35;
                                         description="death probability of an infected nest"),
                               Parameter(:t0, 0.0; description="transmission intercept"),
                               Parameter(:t1, 0.63;
                                         description="transmission per infected neighbour"),
                               Parameter(:a, 0.007;
                                         description="background infection probability"),
                               Parameter(:m0, 0.116;
                                         description="phorid mortality intercept (:phorid)"),
                               Parameter(:m1, 0.036;
                                         description="phorid mortality per occupied neighbour (:phorid)")],
                   description="Jackson et al. (2009) Azteca nests and Lecanicillium, mortality $mortality, $phases phase(s)")
end

"""
    problem_coffee_ant_fungus(; rows=90, cols=120, density=0.08, mortality=:fungus,
                              phases=:two, iterations=1000, seed=0, pmap=(;), kwargs...)
        -> CAProblem

The setup of [Jackson2009](@cite): a 90 × 120 torus of shade trees, each a
nest with probability `density` (drawn from `seed`, all uninfected), for 1000
six-month iterations (2000 steps with two phases).
"""
function problem_coffee_ant_fungus(; rows::Integer=90, cols::Integer=120,
                                   density::Real=0.08, mortality::Symbol=:fungus,
                                   phases::Symbol=:two, iterations::Integer=1000,
                                   seed::Integer=0, pmap=(;), kwargs...)
    0 <= density <= 1 || throw(ArgumentError("density must lie in [0, 1], got $density"))
    rng = Random.Xoshiro(seed)
    tree = Int8[rand(rng) < density ? ANT_NEST : ANT_EMPTY for _ in 1:rows, _ in 1:cols]
    steps = phases === :two ? 2iterations : iterations
    return CAProblem(build_coffee_ant_fungus(; mortality, phases), (; tree), steps, pmap;
                     space=Torus(rows, cols), seed, kwargs...)
end

"""
    ant_iterations(t; phases=:two) -> iterations

The six-month iterations of the source that `t` steps of the ant model make:
`t / 2` with two phases per iteration, `t` with one. Iterations are complete at
even steps.
"""
ant_iterations(t; phases::Symbol=:two) = phases === :two ? t ./ 2 : t

"""
    nest_density(u) -> Float64

The share of shade trees holding a nest, infected or not.
"""
nest_density(u) = count(!=(ANT_EMPTY), u) / length(u)

"""
    mean_variance_ratio(u; quadrat=3) -> Float64

The mean over the variance of the number of nests in non-overlapping
`quadrat × quadrat` squares tiling the lattice (rows and columns beyond the
last whole square are left out): about 1 for nests placed at random and below
1 for clumped nests. `NaN` if every square holds the same number.
"""
function mean_variance_ratio(u; quadrat::Integer=3)
    quadrat >= 1 || throw(ArgumentError("quadrat must be at least 1, got $quadrat"))
    nr, nc = size(u) .÷ quadrat
    counts = [count(!=(ANT_EMPTY),
                    view(u, ((i - 1) * quadrat + 1):(i * quadrat),
                         ((j - 1) * quadrat + 1):(j * quadrat))) for i in 1:nr, j in 1:nc]
    m = sum(counts) / length(counts)
    v = sum(abs2, counts .- m) / (length(counts) - 1)
    return v == 0 ? NaN : m / v
end

_register!(:coffee_ant_fungus; build=build_coffee_ant_fungus,
           problem=problem_coffee_ant_fungus)
