# Peck & Ellner (1997): the evolution of pesticide resistance across a region of fields. Each field holds a
# pest population with a single resistance locus; each generation it grows logistically, pays a cost of
# resistance, is sprayed if it is a treated field above an economic threshold, and exchanges a fixed share
# of its pests with its neighbours. After a random start the model is deterministic. One synchronous rule
# does it all: each cell works out its own and its neighbours' growth and selection, then the migrants.

"""
    build_resistance_regional(; neighbourhood=:von_neumann) -> CAModel

The regional model of [Peck1997](@cite), one generation a step. Each field
holds `N` pests, a frequency `q` of the resistance allele, whether it is
`sprayed` (a treated field rather than a refuge) or `rotated`, its carrying
capacity `K` and its number of `neighbours` in the lattice. In every field at
once:

1. the population grows logistically, `N + ⌊r (1 - N/K) N⌋` (at least 0);
2. genotypes in Hardy-Weinberg proportions pay the cost of resistance
   (fitnesses `1`, `1 - s h` and `1 - s` for susceptible homozygotes,
   heterozygotes and resistant homozygotes);
3. a sprayed field whose population exceeds `delta × Ks` (`Ks = k N0`) is
   sprayed, killing a share `rho` of susceptible homozygotes and
   `rho (1 - h)` of heterozygotes;
4. migrants leave and arrive at once: a share `M` of each field's population
   (rounded down), split equally among its neighbours in the lattice, carrying
   their field's allele frequency (`migration = :total`), or `⌊M N⌋` to each
   neighbour, as the source's equations read literally
   (`migration = :per_neighbour`);
5. a rotated field is emptied every `psi_period` generations.

Neighbours are the von Neumann or Moore neighbours (`neighbourhood`); the
lattice has hard edges.
"""
function build_resistance_regional(; neighbourhood::Symbol=:von_neumann)
    neighbourhood in (:von_neumann, :moore) ||
        throw(ArgumentError("neighbourhood must be :von_neumann or :moore, got :$neighbourhood"))
    nb = neighbourhood === :von_neumann ? VonNeumann(1) : Moore(1)
    rule = LocalRule(nb) do c, nbrs, ctx
        p = ctx.p
        N, q = _resistance_local(c, p)
        out = _resistance_migrants(N, c.neighbours, p) * c.neighbours
        alleles = (N - out) * q
        total = N - out
        for n in nbrs
            n.neighbours > 0 || continue          # beyond the edge
            Nn, qn = _resistance_local(n, p)
            m = _resistance_migrants(Nn, n.neighbours, p)
            total += m
            alleles += m * qn
        end
        q = total > 0 ? alleles / total : q
        p.psi_period > 0 && c.rotated && iszero(mod(ctx.t + 1, p.psi_period)) &&
            (total = 0.0)
        return (; c..., N=total, q)
    end
    states = [StateVariable(:N, 0.0; description="pests in the field"),
              StateVariable(:q, 0.0; description="frequency of the resistance allele"),
              StateVariable(:sprayed, false;
                            description="a treated field (true) or a refuge (false)"),
              StateVariable(:rotated, false; description="a field in crop rotation"),
              StateVariable(:K, 0.0; description="carrying capacity of the field"),
              StateVariable(:neighbours, Int8(0);
                            description="neighbours of the field in the lattice")]
    parameters = [Parameter(:r, 1.5; description="intrinsic growth rate of the logistic"),
                  Parameter(:Ks, 2e7;
                            description="carrying capacity of treated fields, k N0"),
                  Parameter(:delta, 0.45;
                            description="economic threshold, the share of Ks above which fields are sprayed"),
                  Parameter(:rho, 0.75; description="share of susceptibles a spray kills"),
                  Parameter(:s, 5e-4; description="cost of resistance (proportion killed)"),
                  Parameter(:h, 0.05; description="dominance of resistance"),
                  Parameter(:M, 0.05; description="share of a field's pests that migrate"),
                  Parameter(:migration, :total;
                            description=":total (M leaves, split among neighbours) or :per_neighbour (M to each)"),
                  Parameter(:psi_period, 2;
                            description="rotated fields are emptied every psi_period generations (0 for never)")]
    return CAModel(rule; name=:resistance_regional, states, parameters,
                   description="Peck & Ellner (1997) evolution of pesticide resistance across a region")
end

# Growth, cost of resistance and spraying in one field: returns its population and allele frequency.
function _resistance_local(c, p)
    N = c.N
    N > 0 || return (0.0, c.q)
    N = max(N + floor(p.r * (1 - N / c.K) * N), 0.0)
    q = c.q
    P = 1 - q
    GG, Gg, gg = P^2 * N, 2P * q * N, q^2 * N
    Gg *= 1 - p.s * p.h
    gg *= 1 - p.s
    if c.sprayed && GG + Gg + gg > p.delta * p.Ks
        GG *= 1 - p.rho
        Gg *= 1 - p.rho * (1 - p.h)
    end
    N = GG + Gg + gg
    return N > 0 ? (N, (gg + Gg / 2) / N) : (0.0, q)
end

# Migrants a field sends to each of its B neighbours.
function _resistance_migrants(N, B, p)
    B > 0 || return 0.0
    return p.migration === :per_neighbour ? floor(p.M * N) : floor(p.M * N / B)
end

"""
    problem_resistance_regional(; n=50, N0=1e7, k=2, kappa=0.5, q0=0.005, lambda=0.55, psi=0,
                                neighbourhood=:von_neumann, generations=250, seed=0, pmap=(;),
                                kwargs...) -> CAProblem

The setup of [Peck1997](@cite): an `n × n` region of fields, a share `lambda`
of them treated (drawn at random) and the rest refuges, and a share `psi` in
rotation; carrying capacities `Ks = k N0` in treated fields and `kappa Ks` in
refuges; each field starting with `⌊N0 + σ⌋` pests, `σ` uniform on
`[-10⁶, 10⁶]`, and a frequency `q0` of the resistance allele; all drawn from
`seed`, for `generations` generations. The other parameters (`r`, `delta`,
`rho`, `s`, `h`, `M`, `migration`, `psi_period`) default to the values of the
source's figures and are set through `pmap`.
"""
function problem_resistance_regional(; n::Integer=50, N0::Real=1e7, k::Real=2,
                                     kappa::Real=0.5,
                                     q0::Real=0.005, lambda::Real=0.55, psi::Real=0,
                                     neighbourhood::Symbol=:von_neumann,
                                     generations::Integer=250, seed::Integer=0, pmap=(;),
                                     kwargs...)
    n >= 1 || throw(ArgumentError("n must be positive, got $n"))
    0 <= lambda <= 1 || throw(ArgumentError("lambda must lie in [0, 1], got $lambda"))
    0 <= psi <= 1 || throw(ArgumentError("psi must lie in [0, 1], got $psi"))
    0 <= q0 <= 1 || throw(ArgumentError("q0 must lie in [0, 1], got $q0"))
    model = build_resistance_regional(; neighbourhood)
    given = merge((; Ks=Float64(k * N0)), (; (Symbol(a) => v for (a, v) in pairs(pmap))...))
    p = merge(defaults(model), given)
    p.migration in (:total, :per_neighbour) ||
        throw(ArgumentError("migration must be :total or :per_neighbour, got $(repr(p.migration))"))
    rng = Random.Xoshiro(hash((seed, 0, :resistance_regional)))
    cells = n * n
    sprayed = falses(n, n)
    sprayed[Random.randperm(rng, cells)[1:round(Int, lambda * cells)]] .= true
    rotated = falses(n, n)
    rotated[Random.randperm(rng, cells)[1:round(Int, psi * cells)]] .= true
    K = [s ? p.Ks : kappa * p.Ks for s in sprayed]
    N = [max(floor(N0 + 1e6 * (2rand(rng) - 1)), 0.0) for _ in 1:n, _ in 1:n]
    moore = neighbourhood === :moore
    neighbours = [Int8(count(((a, b),) -> (a, b) != (0, 0) && 1 <= i + a <= n &&
                                          1 <= j + b <= n &&
                                          (moore || abs(a) + abs(b) == 1),
                             Iterators.product(-1:1, -1:1)))
                  for i in 1:n, j in 1:n]
    u0 = (; N, q=fill(Float64(q0), n, n), sprayed, rotated, K, neighbours)
    return CAProblem(model, u0, generations, given; space=Grid(n, n), seed, kwargs...)
end

"""
    resistance_summary(u) -> NamedTuple

The source's regional measures of the state `u`: the `mean` and `variance`
over fields of the frequency of the resistance allele (unweighted), and its
Moran's I (`morans_i`) with neighbours sharing an edge (`NaN` when every
field has the same frequency).
"""
function resistance_summary(u)
    q = u.q
    n1, n2 = size(q)
    m = sum(q) / length(q)
    v = sum(abs2, q .- m) / length(q)
    num = 0.0
    W = 0
    for j in 1:n2, i in 1:n1, (a, b) in ((1, 0), (-1, 0), (0, 1), (0, -1))
        (1 <= i + a <= n1 && 1 <= j + b <= n2) || continue
        num += (q[i, j] - m) * (q[i + a, j + b] - m)
        W += 1
    end
    I = !allequal(q) && W > 0 ? (length(q) / W) * num / sum(abs2, q .- m) : NaN
    return (; mean=m, variance=v, morans_i=I)
end

"""
    resistance_time(sol) -> Union{Int,Nothing}

The source's time to resistance: the first saved time at which the regional
mean frequency of the resistance allele ([`resistance_summary`](@ref))
exceeds 0.5, or `nothing` if it never does.
"""
function resistance_time(sol)
    for (t, u) in zip(sol.t, sol.u)
        resistance_summary(u).mean > 0.5 && return Int(t)
    end
    return nothing
end

_register!(:resistance_regional; build=build_resistance_regional,
           problem=problem_resistance_regional)
