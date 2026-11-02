# Zsoldos & Stepler (2024): a Wa-Tor style model of Huffaker's (1958) mite experiment. Prey and
# predatory mites occupy cells of a grid of orange patches, at most one per cell. Each step, the mites
# alive at its start act one at a time in a random order, in place: a predator eats a neighbouring
# prey if it can, otherwise a mite moves to a random empty neighbour; a mite that moves and has
# reached its fertility threshold leaves an offspring in the cell it left; a predator loses energy
# when it moves without eating, and dies at zero. Moves write two cells, which a local rule cannot,
# so the step is a callback stored with the problem and the model's own rule leaves cells unchanged.

const MITE_EMPTY = Int8(0)
const MITE_PREY = Int8(1)
const MITE_PREDATOR = Int8(2)

const VON_NEUMANN_OFFSETS = ((0, -1), (1, 0), (0, 1), (-1, 0))

"""
    HUFFAKER_POSTS

The 18 wooden posts of the source's Grid Design 4, as `(row, column)` cells of
its 30 × 36 grid.
"""
const HUFFAKER_POSTS = [(y + 1, x + 1)
                        for (x, y) in
                            ((4, 7), (4, 16), (4, 25), (7, 4), (7, 13), (7, 22), (16, 7),
                             (16, 16), (16, 25), (19, 4), (19, 13), (19, 22), (28, 7),
                             (28, 16), (28, 25), (31, 4), (31, 13), (31, 22))]

"""
    build_huffaker_mites() -> CAModel

The state of the mite model of [Zsoldos2024](@cite): layers `mite` (0 empty,
1 prey, 2 predator), `fertility` and `energy` (counters of the mite in a
cell) and `post` (whether the cell holds a wooden post). The model's rule
leaves every cell unchanged: the dynamics, in which mites act one at a time,
are the step callback [`huffaker_step_callback`](@ref), which
[`problem_huffaker_mites`](@ref) stores with the problem. The parameters
are `F_prey` and `F_predator` (fertility thresholds), `L_predator` (a
predator's initial energy), `loss` (energy lost by a predator moving without
eating), `gain` (energy gained by eating) and `prey_energy`.
"""
function build_huffaker_mites()
    rule = LocalRule((cell, nbrs, ctx) -> cell, VonNeumann(1))
    return CAModel(rule; name=:huffaker_mites,
                   states=[StateVariable(:mite, MITE_EMPTY; domain=Int8(0):Int8(2),
                                         description="0 empty, 1 prey, 2 predator"),
                           StateVariable(:fertility, Int16(0);
                                         description="steps since the mite last bred"),
                           StateVariable(:energy, Int16(0);
                                         description="the mite's life points"),
                           StateVariable(:post, false; description="a wooden post")],
                   parameters=[Parameter(:F_prey, 1;
                                         description="prey fertility threshold"),
                               Parameter(:F_predator, 10;
                                         description="predator fertility threshold"),
                               Parameter(:L_predator, 1;
                                         description="initial energy of a predator"),
                               Parameter(:loss, 1;
                                         description="energy a predator loses moving without eating"),
                               Parameter(:gain, 1;
                                         description="energy a predator gains by eating"),
                               Parameter(:prey_energy, 10;
                                         description="energy of a prey (never changes)")],
                   description="Zsoldos & Stepler (2024) predator and prey mites on oranges")
end

"""
    huffaker_step!(u, p, rng; posts=:code) -> u

One step of the mite model of [Zsoldos2024](@cite), in place on the layered
state `u`, following the source's code. The mites alive at the start act once
each, in a random order; each one's fertility rises by one; a predator with a
prey among its von Neumann neighbours eats one at random and moves into its
cell, gaining `gain`; otherwise the mite moves to a random empty neighbour, a
predator losing `loss`, or stays if there is none. A predator whose energy
reaches zero dies where it stands. A mite that moved and has reached its
fertility threshold resets it and leaves an offspring in the cell it left. No
cell lies beyond the edges.

Prey may also use the posts: with `posts = :code` (the source's code) every
empty post is among a prey's possible moves, wherever it stands; with
`posts = :text` (the source's description) a prey standing on a post may jump
to a random empty post instead of a neighbour, with the jump counted as one
more possible move; with `posts = :none` posts are ignored.
"""
function huffaker_step!(u, p, rng; posts::Symbol=:code)
    posts in (:code, :text, :none) ||
        throw(ArgumentError("posts must be :code, :text or :none, got :$posts"))
    mite, fert, energy, post = u.mite, u.fertility, u.energy, u.post
    nr, nc = size(mite)
    postcells = posts === :none ? CartesianIndex{2}[] : findall(post)
    # the mites alive at the start, as ids; `where[id]` is a mite's cell and `id_at` the reverse
    where = Random.shuffle!(rng, findall(!=(MITE_EMPTY), mite))
    id_at = zeros(Int, nr, nc)
    for (k, I) in enumerate(where)
        id_at[I] = k
    end
    alive = trues(length(where))
    cand = CartesianIndex{2}[]
    for k in eachindex(where)
        alive[k] || continue
        I = where[k]
        kind = mite[I]
        fert[I] += Int16(1)
        target = I
        moved = false
        if kind == MITE_PREDATOR
            empty!(cand)
            for (di, dj) in VON_NEUMANN_OFFSETS
                J = I + CartesianIndex(di, dj)
                checkbounds(Bool, mite, J) && mite[J] == MITE_PREY && push!(cand, J)
            end
            if !isempty(cand)
                target = rand(rng, cand)
                energy[I] += Int16(p.gain)
                victim = id_at[target]
                victim > 0 && (alive[victim] = false)
                mite[target] = MITE_EMPTY
                id_at[target] = 0
                moved = true
            end
        end
        if !moved
            empty!(cand)
            for (di, dj) in VON_NEUMANN_OFFSETS
                J = I + CartesianIndex(di, dj)
                checkbounds(Bool, mite, J) && mite[J] == MITE_EMPTY && push!(cand, J)
            end
            njump = 0
            if kind == MITE_PREY && posts === :code
                for J in postcells
                    J != I && mite[J] == MITE_EMPTY && !(J in cand) && push!(cand, J)
                end
            elseif kind == MITE_PREY && posts === :text && post[I]
                njump = any(J -> J != I && mite[J] == MITE_EMPTY, postcells) ? 1 : 0
            end
            if length(cand) + njump > 0
                c = rand(rng, 1:(length(cand) + njump))
                target = c <= length(cand) ? cand[c] :
                         rand(rng,
                              [J for J in postcells if J != I && mite[J] == MITE_EMPTY])
                kind == MITE_PREDATOR && (energy[I] -= Int16(p.loss))
                moved = true
            end
        end
        if kind == MITE_PREDATOR && energy[I] <= 0
            mite[I] = MITE_EMPTY
            fert[I] = 0
            energy[I] = 0
            id_at[I] = 0
            alive[k] = false
        elseif moved
            mite[target], fert[target], energy[target] = kind, fert[I], energy[I]
            id_at[target] = k
            where[k] = target
            threshold = kind == MITE_PREY ? p.F_prey : p.F_predator
            if fert[target] >= threshold
                fert[target] = 0
                fert[I] = 0
                energy[I] = kind == MITE_PREY ? p.prey_energy : p.L_predator
                id_at[I] = 0          # the newborn acts from the next step
            else
                mite[I] = MITE_EMPTY
                fert[I] = 0
                energy[I] = 0
                id_at[I] = 0
            end
        end
    end
    return u
end

"""
    huffaker_step_callback(; posts=:code) -> DiscreteCallback

A callback that runs [`huffaker_step!`](@ref) on the state after every step,
with the problem's parameters and a random number generator seeded from the
problem's seed and the time.
"""
function huffaker_step_callback(; posts::Symbol=:code)
    posts in (:code, :text, :none) ||
        throw(ArgumentError("posts must be :code, :text or :none, got :$posts"))
    affect! = function (integ)
        rng = Random.Xoshiro(hash((integ.prob.seed, integ.t, :huffaker_mites)))
        huffaker_step!(integ.u, integ.p, rng; posts)
        return nothing
    end
    return DiscreteCallback((u, t, integ) -> true, affect!; save_positions=(false, false))
end

"""
    problem_huffaker_mites(; experiment=2, prey=120, predators=27, posts=:code, tspan=500,
                           seed=0, pmap=(;), kwargs...) -> CAProblem

One of the four experiments of [Zsoldos2024](@cite), with `prey` prey and
`predators` predators on distinct random cells (drawn from `seed`) for 500
steps:

1. a 20 × 24 grid (four patches per orange), no posts;
2. a 30 × 36 grid (nine patches per orange), no posts;
3. as 2, with the 18 posts ([`HUFFAKER_POSTS`](@ref));
4. as 3, with the fertility thresholds, the predators' initial energy and their
   energy loss doubled (`F_prey = 2`, `F_predator = 20`, `L_predator = 2`,
   `loss = 2`).

The step callback, with the post rule `posts`, is stored with the problem.
"""
function problem_huffaker_mites(; experiment::Integer=2, prey::Integer=120,
                                predators::Integer=27, posts::Symbol=:code, tspan=500,
                                seed::Integer=0, pmap=(;), kwargs...)
    experiment in 1:4 ||
        throw(ArgumentError("experiment must be 1, 2, 3 or 4, got $experiment"))
    dims = experiment == 1 ? (20, 24) : (30, 36)
    prey + predators <= prod(dims) ||
        throw(ArgumentError("$(prey + predators) mites do not fit on $(dims[1]) × $(dims[2]) cells"))
    base = experiment == 4 ? (F_prey=2, F_predator=20, L_predator=2, loss=2) : (;)
    given = (; (Symbol(k) => v for (k, v) in pairs(pmap))...)
    p = merge(defaults(build_huffaker_mites()), base, given)
    rng = Random.Xoshiro(seed)
    mite = fill(MITE_EMPTY, dims)
    energy = zeros(Int16, dims)
    cells = Random.randperm(rng, prod(dims))
    mite[cells[1:prey]] .= MITE_PREY
    energy[cells[1:prey]] .= p.prey_energy
    mite[cells[(prey + 1):(prey + predators)]] .= MITE_PREDATOR
    energy[cells[(prey + 1):(prey + predators)]] .= p.L_predator
    post = fill(false, dims)
    experiment >= 3 && foreach(c -> post[c...] = true, HUFFAKER_POSTS)
    return CAProblem(build_huffaker_mites(), (; mite, energy, post), tspan,
                     merge(base, given);
                     space=Grid(dims...), seed, callback=huffaker_step_callback(; posts),
                     kwargs...)
end

"""
    huffaker_counts(sol) -> (prey, predators)

The numbers of prey and of predators in each saved state of a solution.
"""
function huffaker_counts(sol)
    return ([count(==(MITE_PREY), u.mite) for u in sol.u],
            [count(==(MITE_PREDATOR), u.mite) for u in sol.u])
end

"""
    count_peaks(x; window=10, prominence=0.2) -> Int

The number of peaks of the series `x` after a centred moving average over
`window` points: rises of at least `prominence` times the smoothed series'
maximum above the lowest point since the previous peak, each followed by a
fall of the same size.
"""
function count_peaks(x::AbstractVector; window::Integer=10, prominence::Real=0.2)
    n = length(x)
    n == 0 && return 0
    h = window ÷ 2
    s = [sum(x[max(1, i - h):min(n, i + h)]) / length(max(1, i - h):min(n, i + h))
         for i in 1:n]
    thr = prominence * maximum(s)
    thr > 0 || return 0
    peaks = 0
    lo = hi = s[1]
    rising = true
    for v in s
        if rising
            v < lo && hi - lo < thr && (lo = hi = v)
            v > hi && (hi = v)
            if hi - lo >= thr && hi - v >= thr
                peaks += 1
                rising = false
                lo = v
            end
        else
            v < lo && (lo = v)
            if v - lo >= thr
                rising = true
                hi = v
            end
        end
    end
    return peaks
end

"""
    huffaker_outcome(prey, predators; peaks=3, kwargs...) -> Symbol

The outcome of a run from its prey and predator counts, in the categories of
[Zsoldos2024](@cite), with an oscillation meaning at least `peaks` peaks of
the prey ([`count_peaks`](@ref), with `kwargs`): `:both_extinct_immediately`,
`:both_extinct` (after oscillating), `:prey_extinct`,
`:predator_extinct_immediately`, `:predator_extinct` (after oscillating) or
`:no_extinction`. The source classified its runs by eye.
"""
function huffaker_outcome(prey::AbstractVector, predators::AbstractVector;
                          peaks::Integer=3, kwargs...)
    oscillated = count_peaks(prey; kwargs...) >= peaks
    prey_out, pred_out = last(prey) == 0, last(predators) == 0
    if prey_out && pred_out
        return oscillated ? :both_extinct : :both_extinct_immediately
    elseif prey_out
        return :prey_extinct
    elseif pred_out
        return oscillated ? :predator_extinct : :predator_extinct_immediately
    end
    return :no_extinction
end

_register!(:huffaker_mites; build=build_huffaker_mites, problem=problem_huffaker_mites)
