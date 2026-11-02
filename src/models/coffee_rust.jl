# Beasley, Aristizábal, Bueno & White (2022): coffee leaf rust in fragmented landscapes. Coffee farms
# are cells of a landscape whose other cells carry a forest canopy cover. Within coffee, rust spreads
# by a probabilistic CA on the Moore neighbourhood; between patches, each infected coffee cell on the
# edge of a patch releases a spore into the matrix, which walks until the canopy it crosses uses up its
# movement credit and then infects a susceptible coffee neighbour of its end point. The CA is the
# model's rule; the spores are a step callback, stored with the problem, drawing from a sequential RNG
# seeded from the problem's seed and the time, so results do not depend on the backend.

const RUST_MATRIX = Int8(0)
const RUST_SUSCEPTIBLE = Int8(1)
const RUST_INFECTED = Int8(2)

const MOORE_OFFSETS = ((-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1))

"""
    build_coffee_rust(; ca_rule=:archived) -> CAModel

The within-patch automaton of the coffee rust model of [Beasley2022](@cite).
Each cell is matrix (0), with a canopy cover in the layer `canopy`, or a coffee
farm, susceptible (1) or infected (2). A susceptible farm with `N ≥ 1` infected
Moore neighbours is infected in a step with probability `scale (N + 1) / 10`,
the mean of the authors' code's `scale × Beta(N + 1, 9 - N)` draw
(`ca_rule = :archived`), or `N / 8`, the mean of the paper's `Beta(N, 8 - N)`
(`ca_rule = :paper`). Infected farms stay infected. Nothing lies beyond the
edges.

Between-patch spread by spores is a step callback,
[`coffee_rust_spore_callback`](@ref), which [`problem_coffee_rust`](@ref)
stores with the problem; it reads the parameters `credit` and
`spore_infection`.
"""
function build_coffee_rust(; ca_rule::Symbol=:archived)
    ca_rule in (:archived, :paper) ||
        throw(ArgumentError("ca_rule must be :archived or :paper, got :$ca_rule"))
    archived = ca_rule === :archived
    rule = LocalRule(Moore(1)) do cell, nbrs, ctx
        cell.coffee == RUST_SUSCEPTIBLE || return cell
        N = count(c -> c.coffee == RUST_INFECTED, nbrs)
        N == 0 && return cell
        prob = archived ? ctx.p.scale * (N + 1) / 10 : N / 8
        return bernoulli(ctx, prob, 1) ? (coffee=RUST_INFECTED, canopy=cell.canopy) : cell
    end
    return CAModel(rule; name=:coffee_rust,
                   states=[StateVariable(:coffee, RUST_MATRIX; domain=Int8(0):Int8(2),
                                         description="0 matrix, 1 susceptible coffee, 2 infected coffee"),
                           StateVariable(:canopy, 0.0;
                                         description="canopy cover of a matrix cell")],
                   parameters=[Parameter(:scale, 0.1;
                                         description="scale of the archived CA probability"),
                               Parameter(:credit, 10.0;
                                         description="movement credit of a spore"),
                               Parameter(:spore_infection, 0.5;
                                         description="probability a landing spore infects")],
                   description="Beasley et al. (2022) coffee rust within patches, CA rule $ca_rule")
end

# The admissible moves of a spore at (i, j): in range, matrix, with canopy > 0.
function _spore_moves!(moves, coffee, canopy, i, j)
    empty!(moves)
    n, m = size(coffee)
    for (di, dj) in MOORE_OFFSETS
        a, b = i + di, j + dj
        (1 <= a <= n && 1 <= b <= m) || continue
        coffee[a, b] == RUST_MATRIX && canopy[a, b] > 0 && push!(moves, (a, b))
    end
    return moves
end

"""
    coffee_rust_spores!(coffee, canopy, rng; credit=10.0, spore_infection=0.5) -> Int

One round of between-patch spread in the coffee rust model of
[Beasley2022](@cite), in place: every infected farm with a matrix neighbour
releases one spore into one of them, chosen at random; each spore in turn
moves to a random admissible Moore neighbour (in range, matrix, canopy above
zero), paying the canopy cover of each cell it enters, until its `credit` is
spent or it cannot move; it then picks a random susceptible farm next to its
end point, if any, and infects it with probability `spore_infection`. Returns
the number of farms infected.
"""
function coffee_rust_spores!(coffee::AbstractMatrix, canopy::AbstractMatrix, rng;
                             credit::Real=10.0, spore_infection::Real=0.5)
    n, m = size(coffee)
    starts = Tuple{Int,Int}[]
    nb = Tuple{Int,Int}[]
    for j in 1:m, i in 1:n
        coffee[i, j] == RUST_INFECTED || continue
        empty!(nb)
        for (di, dj) in MOORE_OFFSETS
            a, b = i + di, j + dj
            (1 <= a <= n && 1 <= b <= m) && coffee[a, b] == RUST_MATRIX && push!(nb, (a, b))
        end
        isempty(nb) || push!(starts, rand(rng, nb))
    end
    moves = Tuple{Int,Int}[]
    infected = 0
    for (i, j) in starts
        left = Float64(credit)
        while left > 0
            isempty(_spore_moves!(moves, coffee, canopy, i, j)) && break
            i, j = rand(rng, moves)
            left -= canopy[i, j]
        end
        empty!(nb)
        for (di, dj) in MOORE_OFFSETS
            a, b = i + di, j + dj
            (1 <= a <= n && 1 <= b <= m) && coffee[a, b] == RUST_SUSCEPTIBLE &&
                push!(nb, (a, b))
        end
        isempty(nb) && continue
        target = rand(rng, nb)
        if rand(rng) < spore_infection
            coffee[target...] = RUST_INFECTED
            infected += 1
        end
    end
    return infected
end

"""
    coffee_rust_spores_archived!(coffee, canopy, rng; credit=10.0, spore_infection=0.5,
                                 stop=99, limit=50) -> Int

One round of between-patch spread as the authors' archived code
(`LandscapeSim.py`) performs it, errors included, in place; written from a
description of that code. With 0-based indices `(r, c)` on an `n × n`
landscape:

1. every infected farm, in row-major order, with a neighbour that is matrix or
   off the landscape picks one such neighbour at random, and keeps the spore
   there only if `r < (n & c) < n` (Python's `a < b & c < d`, with `&` the
   bitwise and), so spores off the top edge can be kept;
2. each spore moves to random Moore neighbours that are matrix cells with
   canopy above zero, paying their canopy, while its credit is positive; it
   stops on entering a row or column above `stop`; a spore with no such
   neighbour moves by `(1, 1)` (Python 3.8's `random.choices` with zero
   weights), and entering coffee ends its walk;
3. then each spore with a susceptible neighbour of its own end point picks one
   at random, but the target is that offset from the *last* spore's end point;
   if both target indices are below `limit` (negative ones wrap around, as in
   NumPy) the target is infected with probability `spore_infection` if it is
   susceptible.

Returns the number of farms infected.
"""
function coffee_rust_spores_archived!(coffee::AbstractMatrix, canopy::AbstractMatrix, rng;
                                      credit::Real=10.0, spore_infection::Real=0.5,
                                      stop::Integer=99, limit::Integer=50)
    n = size(coffee, 1)
    size(coffee, 2) == n ||
        throw(ArgumentError("the archived spores need a square landscape"))
    inside(r, c) = 0 <= r < n && 0 <= c < n
    at(A, r, c) = A[mod(r, n) + 1, mod(c, n) + 1]   # 0-based, wrapping like NumPy
    spores = Tuple{Int,Int}[]
    free = Int[]
    for r in 0:(n - 1), c in 0:(n - 1)
        at(coffee, r, c) == RUST_INFECTED || continue
        empty!(free)
        for (k, (dr, dc)) in enumerate(MOORE_OFFSETS)
            a, b = r + dr, c + dc
            (!inside(a, b) || at(coffee, a, b) == RUST_MATRIX) && push!(free, k)
        end
        isempty(free) && continue
        dr, dc = MOORE_OFFSETS[rand(rng, free)]
        a, b = r + dr, c + dc
        a < (n & b) < n && push!(spores, (a, b))
    end
    isempty(spores) && return 0
    moves = Int[]
    for s in eachindex(spores)
        r, c = spores[s]
        left = Float64(credit)
        while left > 0
            empty!(moves)
            for (k, (dr, dc)) in enumerate(MOORE_OFFSETS)
                a, b = r + dr, c + dc
                inside(a, b) && at(coffee, a, b) == RUST_MATRIX && at(canopy, a, b) > 0 &&
                    push!(moves, k)
            end
            dr, dc = isempty(moves) ? MOORE_OFFSETS[end] : MOORE_OFFSETS[rand(rng, moves)]
            r, c = r + dr, c + dc
            (r > stop || c > stop) && break
            left -= at(coffee, r, c) == RUST_MATRIX ? at(canopy, r, c) : NaN
        end
        spores[s] = (r, c)
    end
    lr, lc = spores[end]
    infected = 0
    targets = Int[]
    for (r, c) in spores
        empty!(targets)
        for (k, (dr, dc)) in enumerate(MOORE_OFFSETS)
            a, b = r + dr, c + dc
            inside(a, b) && at(coffee, a, b) == RUST_SUSCEPTIBLE && push!(targets, k)
        end
        isempty(targets) && continue
        dr, dc = MOORE_OFFSETS[rand(rng, targets)]
        a, b = lr + dr, lc + dc
        (a < limit && b < limit) || continue
        if rand(rng) < spore_infection && at(coffee, a, b) == RUST_SUSCEPTIBLE
            coffee[mod(a, n) + 1, mod(b, n) + 1] = RUST_INFECTED
            infected += 1
        end
    end
    return infected
end

"""
    coffee_rust_spore_callback(; archived=false) -> DiscreteCallback

A callback that runs [`coffee_rust_spores!`](@ref) (or, with
`archived = true`, [`coffee_rust_spores_archived!`](@ref)) on the state after
every step, with the parameters `credit` and `spore_infection` of the problem
and a random number generator seeded from the problem's seed and the time.
"""
function coffee_rust_spore_callback(; archived::Bool=false)
    affect! = function (integ)
        rng = Random.Xoshiro(hash((integ.prob.seed, integ.t, :coffee_rust_spores)))
        spores! = archived ? coffee_rust_spores_archived! : coffee_rust_spores!
        spores!(integ.u.coffee, integ.u.canopy, rng; integ.p.credit,
                integ.p.spore_infection)
        return nothing
    end
    return DiscreteCallback((u, t, integ) -> true, affect!; save_positions=(false, false))
end

# ------------------------------------------------------------------ landscapes

"""
    classify_landscape(a, weights) -> Matrix{Int}

Classes `0, 1, ...` of the values of `a`, the lowest `weights[1] / sum(weights)`
of them in class 0 and so on, as NLMpy's `classifyArray` assigns them.
"""
function classify_landscape(a::AbstractMatrix, weights)
    s = sort(vec(a))
    cp = cumsum(weights) ./ sum(weights)
    bounds = [s[clamp(floor(Int, c * length(s)), 1, length(s))] for c in cp]
    bounds[end] = max(bounds[end], maximum(a))
    return [searchsortedfirst(bounds, x) - 1 for x in a]
end

# Fill every cell with the value of its nearest seed (Euclidean), then rescale to [0, 1].
function _nearest_fill(seeds::Vector{CartesianIndex{2}}, values, dims)
    out = zeros(dims)
    for I in CartesianIndices(dims)
        best, d = 1, typemax(Int)
        for (k, S) in enumerate(seeds)
            e = sum(abs2, Tuple(I - S))
            e < d && ((best, d) = (k, e))
        end
        out[I] = values[best]
    end
    lo, hi = extrema(out)
    return hi > lo ? (out .- lo) ./ (hi - lo) : zero(out)
end

"""
    random_cluster_landscape(rng, n, p) -> Matrix{Float64}

An `n × n` neutral landscape by the modified random cluster method, as NLMpy's
`randomClusterNN` with the Moore neighbourhood: a fraction `p` of cells is
switched on at random, each Moore-connected cluster of them gets a uniform
random value, every other cell takes the value of its nearest cluster cell, and
the result is rescaled to `[0, 1]`. Larger `p` gives larger patches.
"""
function random_cluster_landscape(rng, n::Integer, p::Real)
    0 < p < 1 || throw(ArgumentError("p must lie in (0, 1), got $p"))
    on = classify_landscape(rand(rng, n, n), [1 - p, p]) .== 1
    labels = zeros(Int, n, n)
    nlab = 0
    stack = CartesianIndex{2}[]
    for I in CartesianIndices(on)
        (on[I] && labels[I] == 0) || continue
        nlab += 1
        labels[I] = nlab
        push!(stack, I)
        while !isempty(stack)
            J = pop!(stack)
            for (di, dj) in MOORE_OFFSETS
                K = J + CartesianIndex(di, dj)
                checkbounds(Bool, on, K) && on[K] && labels[K] == 0 || continue
                labels[K] = nlab
                push!(stack, K)
            end
        end
    end
    value = rand(rng, nlab)
    seeds = findall(on)
    return _nearest_fill(seeds, [value[labels[S]] for S in seeds], (n, n))
end

"""
    random_element_landscape(rng, n, k) -> Matrix{Float64}

An `n × n` neutral landscape by nearest-neighbour interpolation of `k` random
elements, as NLMpy's `randomElementNN`: `k` distinct random cells get the
values `1, ..., k`, every other cell takes the value of its nearest one, and
the result is rescaled to `[0, 1]`. More elements give smaller patches.
"""
function random_element_landscape(rng, n::Integer, k::Integer)
    1 <= k <= n^2 || throw(ArgumentError("k must lie in 1:$(n^2), got $k"))
    cells = CartesianIndices((n, n))[Random.randperm(rng, n^2)[1:k]]
    return _nearest_fill(collect(cells), collect(1.0:k), (n, n))
end

"""
    coffee_rust_landscape(rng; n=120, clustering=0.2, deforestation=0.45, dispersion=3,
                          coffee=0.25) -> (coffee, canopy)

A simulated landscape of [Beasley2022](@cite): coffee (state 1) on the lowest
`coffee` share of a random cluster landscape with parameter `clustering`, and,
elsewhere, canopy cover from a random element landscape with
`2000 dispersion` elements, of which the lowest `deforestation` share is
deforested (canopy uniform on `[0, 0.3]`) and the rest forest (uniform on
`[0.6, 1]`). Coffee cells have canopy 0.
"""
function coffee_rust_landscape(rng; n::Integer=120, clustering::Real=0.2,
                               deforestation::Real=0.45, dispersion::Integer=3,
                               coffee::Real=0.25)
    0 < deforestation < 1 ||
        throw(ArgumentError("deforestation must lie in (0, 1), got $deforestation"))
    0 < coffee < 1 || throw(ArgumentError("coffee must lie in (0, 1), got $coffee"))
    farms = classify_landscape(random_cluster_landscape(rng, n, clustering),
                               [coffee, 1 - coffee]) .== 0
    k = min(2000 * dispersion, n^2)
    forest = classify_landscape(random_element_landscape(rng, n, k),
                                [deforestation, 1 - deforestation]) .== 1
    canopy = [farms[I] ? 0.0 : forest[I] ? 0.6 + 0.4rand(rng) : 0.3rand(rng)
              for I in CartesianIndices(farms)]
    return Int8.(farms), canopy
end

"""
    problem_coffee_rust(; n=120, buffer=10, clustering=0.2, deforestation=0.45, dispersion=3,
                        ca_rule=:archived, spores=true, tspan=1000, seed=0, pmap=(;),
                        kwargs...) -> CAProblem

The simulated-landscape setup of [Beasley2022](@cite): an `n × n` landscape
from [`coffee_rust_landscape`](@ref), drawn from `seed`, with one coffee farm
infected, for 1000 steps. `spores` is `:intended` (or `true`) for spores as
the paper describes them ([`coffee_rust_spores!`](@ref)), `:archived` for
spores as the authors' archived code moves them
([`coffee_rust_spores_archived!`](@ref)), or `:none` (or `false`). The first
infected farm is chosen at random inside the `buffer`, or, with
`start = :anywhere` (the default with archived spores, as in that code),
anywhere in the landscape.
"""
function problem_coffee_rust(; n::Integer=120, buffer::Integer=10, clustering::Real=0.2,
                             deforestation::Real=0.45, dispersion::Integer=3,
                             ca_rule::Symbol=:archived, spores::Union{Bool,Symbol}=true,
                             start::Symbol=spores === :archived ? :anywhere : :inner,
                             tspan=1000, seed::Integer=0, pmap=(;), kwargs...)
    n > 2buffer || throw(ArgumentError("the landscape must be larger than its buffers"))
    mode = spores === true ? :intended : spores === false ? :none : spores
    mode in (:intended, :archived, :none) ||
        throw(ArgumentError("spores must be :intended, :archived or :none, got :$spores"))
    start in (:inner, :anywhere) ||
        throw(ArgumentError("start must be :inner or :anywhere, got :$start"))
    rng = Random.Xoshiro(seed)
    coffee, canopy = coffee_rust_landscape(rng; n, clustering, deforestation, dispersion)
    region = start === :inner ?
             CartesianIndices(((buffer + 1):(n - buffer), (buffer + 1):(n - buffer))) :
             CartesianIndices(coffee)
    farms = [I for I in region if coffee[I] == RUST_SUSCEPTIBLE]
    isempty(farms) && throw(ArgumentError("the landscape has no coffee to infect"))
    coffee[rand(rng, farms)] = RUST_INFECTED
    extra = mode === :none ? (;) :
            (; callback=coffee_rust_spore_callback(; archived=mode === :archived))
    return CAProblem(build_coffee_rust(; ca_rule), (; coffee, canopy), tspan, pmap;
                     space=Grid(n, n), seed, extra..., kwargs...)
end

"""
    coffee_rust_prevalence(u; buffer=10) -> Float64

The share of coffee farms inside the `buffer` that are infected.
"""
function coffee_rust_prevalence(u; buffer::Integer=10)
    n, m = size(u.coffee)
    inner = view(u.coffee, (buffer + 1):(n - buffer), (buffer + 1):(m - buffer))
    farms = count(!=(RUST_MATRIX), inner)
    return farms == 0 ? NaN : count(==(RUST_INFECTED), inner) / farms
end

"""
    coffee_rust_rate(sol; buffer=10) -> Float64

The rate of spread of a run: new infections per step among the coffee farms
inside the `buffer`, from the first to the last saved state. (The source's
rate multiplies the inner prevalence by the coffee count of the whole
landscape; this one counts the inner farms only.)
"""
function coffee_rust_rate(sol; buffer::Integer=10)
    count_inner(u) = (n=size(u.coffee, 1); m=size(u.coffee, 2);
                      count(==(RUST_INFECTED),
                            view(u.coffee, (buffer + 1):(n - buffer),
                                 (buffer + 1):(m - buffer))))
    return (count_inner(sol.u[end]) - count_inner(sol.u[1])) / (sol.t[end] - sol.t[1])
end

_register!(:coffee_rust; build=build_coffee_rust, problem=problem_coffee_rust)
