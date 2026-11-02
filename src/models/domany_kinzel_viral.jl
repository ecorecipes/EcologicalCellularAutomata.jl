# Capitán, Cuesta, Manrubia & Aguirre (2011): a multicomponent Domany-Kinzel automaton for the front
# of a viral infection spreading through host tissue. The front is a ring of L cells; time is the
# second spatial dimension. Each cell holds the replicative ability r in 0:R of the virion in it, and
# is infected at the next generation from the offspring of its two parents, cells i and i + 1.

"""
    build_domany_kinzel_viral(; R=6) -> CAModel

The multicomponent Domany–Kinzel automaton of [Capitan2011](@cite), for the
front of a viral infection in host tissue: a one-dimensional ring of cells,
updated synchronously, whose successive generations are successive rows of
tissue. A cell holds the replicative ability `r ∈ 0:R` of the virion in it, the
number of infective offspring it produces (0 for an empty cell). Cell `i` of the
next generation receives the `n = r₁ + r₂` offspring of cells `i` and `i + 1`:

- with probability `π^n` every virion fails to infect it (host resistance), and
  it is empty;
- otherwise one virion infects it, an offspring of the first parent with
  probability `r₁ / n`, and the offspring's ability is that of its parent minus
  one with probability `p`, plus one with probability `q` (only below `R`), and
  unchanged otherwise.

`R` is fixed when the model is built; `π`, `p` and `q` are parameters. With
`p = 0` and every cell at `R` this is directed bond percolation, and with only
classes 0 and 1 the Domany–Kinzel automaton.
"""
function build_domany_kinzel_viral(; R::Integer=6)
    1 <= R <= typemax(Int8) - 1 || throw(ArgumentError("R must lie in 1:126, got $R"))
    Rmax = Int8(R)
    rule = LocalRule(Offsets([CartesianIndex(1)])) do r1, nbrs, ctx
        p = ctx.p
        r2 = nbrs[1]
        n = Int(r1) + Int(r2)
        n == 0 && return Int8(0)
        bernoulli(ctx, p.π^n, 1) && return Int8(0)                  # resisted
        r = uniform(ctx, 2) * n < r1 ? r1 : r2                       # the parent
        u = uniform(ctx, 3)                                          # its offspring's mutation
        u < p.p && return r - Int8(1)
        return r < Rmax && u < p.p + p.q ? r + Int8(1) : r
    end
    return CAModel(rule; name=:domany_kinzel_viral,
                   states=[StateVariable(:r, Rmax; domain=Int8(0):Rmax,
                                         description="replicative ability (0: empty)")],
                   parameters=[Parameter(:π, 0.1;
                                         description="host resistance: probability that one virion fails to infect"),
                               Parameter(:p, 0.2;
                                         description="probability of a deleterious mutation"),
                               Parameter(:q, 0.01;
                                         description="probability of a beneficial mutation")],
                   description="Capitán et al. (2011) multicomponent Domany-Kinzel automaton, R = $R")
end

"""
    problem_domany_kinzel_viral(; L=10_000, R=6, tspan=1000, seed=0, pmap=(;), kwargs...)
        -> CAProblem

The setup of [Capitan2011](@cite): a ring of `L` cells, every one infected by a
virion of the greatest ability `R`. Solver keywords (`saveat`, `callback`, ...)
are stored with the problem.
"""
function problem_domany_kinzel_viral(; L::Integer=10_000, R::Integer=6, tspan=1000,
                                     seed::Integer=0, pmap=(;), kwargs...)
    L >= 2 || throw(ArgumentError("the ring needs at least 2 cells, got $L"))
    return CAProblem(build_domany_kinzel_viral(; R), (;), tspan, pmap; space=Torus(L), seed,
                     kwargs...)
end

"""
    viral_density(u) -> Float64

The density of active cells, those holding a virion able to replicate
(`r ≥ 1`): the order parameter ρ of [Capitan2011](@cite), zero at extinction.
"""
viral_density(u) = count(>(0), u) / length(u)

"""
    viral_mean_ability(u) -> Float64

The mean replicative ability ⟨r⟩ of the active cells (`NaN` when there are none).
"""
function viral_mean_ability(u)
    n = count(>(0), u)
    return n == 0 ? NaN : sum(Int, u) / n
end

"""
    viral_diversity(u) -> Float64

The phenotypic diversity `D`, the standard deviation of the replicative
ability among active cells (`NaN` when there are none).
"""
function viral_diversity(u)
    m = viral_mean_ability(u)
    isnan(m) && return NaN
    active = [Int(r) for r in u if r > 0]
    return sqrt(sum((r - m)^2 for r in active) / length(active))
end

"""
    viral_transition(r, r1, r2; π, p, q, R) -> Float64

The probability that a cell of the next generation holds ability `r` when its
parents hold `r1` and `r2`: equations 2 and 3 of [Capitan2011](@cite), from
which the automaton draws.
"""
function viral_transition(r::Integer, r1::Integer, r2::Integer; π::Real, p::Real, q::Real,
                          R::Integer)
    n = r1 + r2
    n == 0 && return Float64(r == 0)
    mutate(ri) = r == ri - 1 ? p :
                 r == ri ? (ri < R ? 1 - p - q : 1 - p) :
                 (r == ri + 1 && ri < R ? q : 0.0)
    return π^n * (r == 0) + (1 - π^n) * (r1 * mutate(r1) + r2 * mutate(r2)) / n
end

_register!(:domany_kinzel_viral; build=build_domany_kinzel_viral,
           problem=problem_domany_kinzel_viral)
