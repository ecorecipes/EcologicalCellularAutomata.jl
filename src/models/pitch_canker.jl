# Möykkynen, Capretti & Pukkala (2015): the spread of pitch canker (Fusarium circinatum) through the pine
# and Douglas fir forests of Europe, on cells of 1 km in yearly steps. Infected seedling lots carried from
# entry points release spores where they are planted; infected trees spread the fungus for some years, by
# airborne spores to their own and nearby cells and by insects flying to other trees, and then die. The
# source's host and climate rasters and its program (PitchCanker) are not published, so problems run on a
# synthetic landscape.
#
# A year is one step: the rule ages the infected trees (last year's new infections start spreading, the
# oldest die) and adds the infections by spores from the cells within reach (Eqs. 2-4 of the source); a
# callback stored with the problem then adds those by arriving seedling lots and by insects, one tree at a
# time (Eq. 1). The state at time t is the end of year t.

# The Lanczos approximation (g = 7, nine terms) of the gamma function, for the scale of a Weibull flight
# distance with a given mean; Random is the zoo's only numerical dependency.
const _PC_LANCZOS = (0.99999999999980993, 676.5203681218851, -1259.1392167224028,
                     771.32342877765313, -176.61502916214059, 12.507343278686905,
                     -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7)

function _pc_gamma(x::Real)
    x < 0.5 && return π / (sinpi(x) * _pc_gamma(1 - x))
    x -= 1
    a = _PC_LANCZOS[1]
    for k in 2:9
        a += _PC_LANCZOS[k] / (x + k - 1)
    end
    s = x + 7.5
    return sqrt(2π) * s^(x + 0.5) * exp(-s) * a
end

# log Γ(x) for x ≥ 1: Lanczos below 10, Stirling's series above (absolute error below 1e-12).
function _pc_lgamma(x::Real)
    x < 10 && return log(_pc_gamma(x))
    return (x - 0.5) * log(x) - x + 0.5 * log(2π) + 1 / (12x) - 1 / (360x^3) + 1 / (1260x^5)
end

# A binomial(n, p) draw: inversion when n min(p, 1 - p) < 10, otherwise the transformed rejection with
# squeeze (BTRS) of Hörmann (1993), both exact.
function _pc_binomial(rng, n::Integer, p::Real)
    (n <= 0 || p <= 0) && return 0
    p >= 1 && return Int(n)
    p > 0.5 && return Int(n) - _pc_binomial(rng, n, 1 - p)
    q = 1 - p
    if n * p < 10
        s = p / q
        a = (n + 1) * s
        while true
            r = q^n
            u = rand(rng)
            x = 0
            while u > r
                u -= r
                x += 1
                x > n && break
                r *= a / x - s
            end
            x <= n && return x
        end
    end
    spq = sqrt(n * p * q)
    b = 1.15 + 2.53 * spq
    a = -0.0873 + 0.0248 * b + 0.01 * p
    c = n * p + 0.5
    vr = 0.92 - 4.2 / b
    α = (2.83 + 5.1 / b) * spq
    lpq = log(p / q)
    m = floor(Int, (n + 1) * p)
    h = _pc_lgamma(m + 1) + _pc_lgamma(n - m + 1)
    while true
        u = rand(rng) - 0.5
        v = rand(rng)
        us = 0.5 - abs(u)
        k = floor(Int, (2a / us + b) * u + c)
        (k < 0 || k > n) && continue
        us >= 0.07 && v <= vr && return k
        v = log(v * α / (a / us^2 + b))
        v <= h - _pc_lgamma(k + 1) - _pc_lgamma(n - k + 1) + (k - m) * lpq && return k
    end
end

# Most host trees in a cell: trees per hectare times the cell's area in hectares.
_pc_nmax(p) = p.trees_per_ha * p.cell^2 / 10_000

# The share of the way from the current to the future climate in year y (y = 1 is the first year).
_pc_climate_weight(p, y) = clamp((y - 1) / (p.climate_year - p.first_year), 0.0, 1.0)

# A negative exponential distance of mean μ, redrawn until it is at most M.
function _pc_truncexp(rng, μ, M)
    while true
        d = -μ * log1p(-rand(rng))
        d <= M && return d
    end
end

# An insect's flight distance: negative exponential, or Weibull of shape c with the same mean, redrawn
# until it is at most insect_max.
function _pc_flight(rng, p)
    p.insect_distribution === :exponential &&
        return _pc_truncexp(rng, p.insect_mean, p.insect_max)
    c = p.insect_shape
    b = p.insect_mean / _pc_gamma(1 + 1 / c)
    while true
        d = b * (-log1p(-rand(rng)))^(1 / c)
        d <= p.insect_max && return d
    end
end

# The distribution function of an insect's flight distance, truncated at insect_max.
function _pc_flight_cdf(r, p)
    r <= 0 && return 0.0
    F = if p.insect_distribution === :exponential
        x -> -expm1(-x / p.insect_mean)
    else
        b = p.insect_mean / _pc_gamma(1 + 1 / p.insect_shape)
        x -> -expm1(-(x / b)^p.insect_shape)
    end
    M = p.insect_max
    isfinite(M) || return F(r)
    return r >= M ? 1.0 : F(r) / F(M)
end

# The probability that a flight from the midpoint of a cell lands in each cell within reach,
# (1 / 2π) ∫ [F(r_out(θ)) - F(r_in(θ))] dθ over the directions θ, where [r_in, r_out] is the part of
# the ray in direction θ inside the cell (a convex square), by the midpoint rule on `steps` directions.
# The reach is insect_max, or without a maximum the distance beyond which a flight lands with
# probability below 1e-6 (such flights are lost). Returns the matrix and its radius R in cells
# (offset (a, b) at [a+R+1, b+R+1]).
function _pc_insect_kernel(p; steps::Integer=4096)
    c = p.cell
    reach = p.insect_max
    if !isfinite(reach)
        reach = p.insect_mean
        while 1 - _pc_flight_cdf(reach, p) > 1e-6
            reach *= 2
        end
    end
    R = ceil(Int, reach / c + 0.5)
    K = zeros(2R + 1, 2R + 1)
    for k in 1:steps
        θ = (k - 0.5) * 2π / steps
        dx, dy = cos(θ), sin(θ)
        for b in (-R):R, a in (-R):R
            lo, hi = 0.0, Inf
            for (d, o) in ((dx, a), (dy, b))
                l, h = (o - 0.5) * c, (o + 0.5) * c
                if abs(d) < 1e-15
                    (l <= 0 <= h) || (hi = -1.0)
                else
                    t1, t2 = minmax(l / d, h / d)
                    lo, hi = max(lo, t1), min(hi, t2)
                end
            end
            hi > lo || continue
            K[a + R + 1, b + R + 1] += _pc_flight_cdf(hi, p) - _pc_flight_cdf(lo, p)
        end
    end
    return K ./ steps, R
end

# The insects of all cells, pooled cell by cell. With the healthy trees as they stand when a cell's
# insects set out, each try lands at offset o with probability K(o) and infects there with probability
# q(o) = (N_healthy / N_max) EI, so an insect infects at o, on some try of `insect_tries`, with
# probability K(o) q(o) (1 - (1 - s)^T) / s, where s = Σ K q; the cell's insects then infect the cells
# around it by a multinomial draw (a binomial per cell), at most every healthy tree of each.
function _pc_insects_pooled!(u, rng, S, p, w, t)
    K, R = _pc_insect_kernel(p)
    nmax = _pc_nmax(p)
    T = p.insect_tries
    Q = zeros(size(K))
    for J in CartesianIndices(S)
        S[J] > 0 || continue
        insects = _pc_round(rng, S[J] * p.insects_per_tree)
        insects > 0 || continue
        s = 0.0
        for I in CartesianIndices(K)
            i, j = J[1] + I[1] - R - 1, J[2] + I[2] - R - 1
            q = 0.0
            if K[I] > 0 && checkbounds(Bool, u.hosts, i, j)
                ei = u.ei_now[i, j] + w * (u.ei_future[i, j] - u.ei_now[i, j])
                q = K[I] * u.hosts[i, j] / nmax * ei
            end
            Q[I] = q
            s += q
        end
        s > 0 || continue
        scale = (1 - (1 - s)^T) / s
        left, rest = insects, 1.0
        for I in CartesianIndices(K)
            π_o = Q[I] * scale
            π_o > 0 || continue
            x = _pc_binomial(rng, left, min(1.0, π_o / rest))
            rest -= π_o
            left -= x
            if x > 0
                i, j = J[1] + I[1] - R - 1, J[2] + I[2] - R - 1
                x = min(x, Int(u.hosts[i, j]))
                u.hosts[i, j] -= Int32(x)
                u.new[i, j] += Int32(x)
                x > 0 && u.infected_year[i, j] == 0 && (u.infected_year[i, j] = Int32(t))
            end
            left == 0 && break
        end
    end
    return u
end

# floor(x), plus one with probability x - floor(x).
_pc_round(rng, x) = floor(Int, x) + (rand(rng) < x - floor(x))

# The layer names of a model with Y spreading cohorts.
@generated function _pc_names(::Val{Y}) where {Y}
    return (:hosts, :coverage, :ei_now, :ei_future, :new, (Symbol(:i, k) for k in 1:Y)...,
            :dead, :infected_year)
end

# The trees spreading in the coming year from a cell at the end of a year (last year's infections and all
# but the oldest cohort): the cell's new infections and its cohorts 1 to Y - 1.
@generated function _pc_next_spreading(c, ::Val{Y}) where {Y}
    terms = [:(Int(c.$(Symbol(:i, k)))) for k in 1:(Y - 1)]
    return :(Int(c.new) + $(Expr(:call, :+, 0, terms...)))
end

# The cohorts after a year: this year's infections become cohort 1, cohort k becomes k + 1.
@generated function _pc_aged(c, ::Val{Y}) where {Y}
    terms = [k == 1 ? :(c.new) : :(c.$(Symbol(:i, k - 1))) for k in 1:Y]
    return :(tuple($(terms...)))
end

# The oldest cohort, which dies at the end of the year.
@generated _pc_oldest(c, ::Val{Y}) where {Y} = :(c.$(Symbol(:i, Y)))

"""
    build_pitch_canker(; spreading_years=5, spore_max=1300, cell=1000) -> CAModel

The pitch canker model of [Moykkynen2015](@cite), in yearly steps on square
cells of side `cell` metres. Each cell holds its healthy host trees (`hosts`),
its host coverage (`coverage`, 0 to 1), the ecoclimatic index of the current
and the future climate (`ei_now`, `ei_future`, 0 to 1), the trees infected
this year (`new`), those infected 1 to `spreading_years` years ago (`i1`,
`i2`, ...), which spread the fungus, the dead trees (`dead`) and the year of
its first infection (`infected_year`, 0 for none). With `N_max` the most host
trees in a cell (`trees_per_ha` times the cell's area in hectares) and `EI`
the index in year `y`, interpolated linearly from `ei_now` in the first year
(`first_year`) to `ei_future` in `climate_year`, each year, in every cell at
once:

1. the infected trees age: this year's infections join those spreading, and
   those infected `spreading_years` years ago die;
2. spores from every cell whose midpoint lies within `spore_max` metres
   (the cell itself included) give a density `Σ (N_spreading / N_max)
   exp(-d / spore_mean)` over those cells at distance `d` (Eqs. 2 and 3), and
   `N_healthy EI density` healthy trees are infected (Eq. 4): its whole part,
   and one more with probability its fractional part, at most every healthy
   tree.

Arrivals of infected seedling lots and the flights of insects are
[`pitch_canker_spread!`](@ref), a callback stored with the problem.
"""
function build_pitch_canker(; spreading_years::Integer=5, spore_max::Real=1300,
                            cell::Real=1000)
    spreading_years >= 1 ||
        throw(ArgumentError("spreading_years must be at least 1, got $spreading_years"))
    cell > 0 || throw(ArgumentError("cell must be positive, got $cell"))
    spore_max >= 0 || throw(ArgumentError("spore_max must not be negative, got $spore_max"))
    R = floor(Int, spore_max / cell)
    within = [(i, j) for j in (-R):R
              for i in (-R):R
              if (i, j) != (0, 0) && hypot(i, j) * cell <= spore_max]
    # with no cell but the donor's own within reach, one neighbour read with weight zero;
    # each name is assigned once, so the rule's closure captures concrete values
    offs = isempty(within) ? [(1, 0)] : within
    dists = isempty(within) ? (Inf,) : Tuple(hypot(i, j) * cell for (i, j) in within)
    coh = ntuple(k -> Symbol(:i, k), spreading_years)
    V = Val(spreading_years)
    rule = LocalRule(Offsets(offs...)) do c, nbrs, ctx
        p = ctx.p
        y = ctx.t + 1
        w = _pc_climate_weight(p, y)
        ei = c.ei_now + w * (c.ei_future - c.ei_now)
        nmax = _pc_nmax(p)
        density = _pc_next_spreading(c, V) / nmax
        for k in eachindex(dists)
            dists[k] <= p.spore_max || continue
            density += _pc_next_spreading(nbrs[k], V) / nmax * exp(-dists[k] / p.spore_mean)
        end
        λ = c.hosts * ei * density
        new = λ > 0 ? min(floor(Int, λ) + (uniform(ctx, 1) < λ - floor(λ)), Int(c.hosts)) :
              0
        year = new > 0 && c.infected_year == 0 ? Int32(y) : c.infected_year
        return NamedTuple{_pc_names(V)}((Int32(c.hosts - new), c.coverage, c.ei_now,
                                         c.ei_future, Int32(new), _pc_aged(c, V)...,
                                         Int32(c.dead + _pc_oldest(c, V)), year))
    end
    states = [StateVariable(:hosts, Int32(0); description="healthy host trees"),
              StateVariable(:coverage, 0.0; description="host coverage (0 to 1)"),
              StateVariable(:ei_now, 0.0;
                            description="ecoclimatic index of the current climate (0 to 1)"),
              StateVariable(:ei_future, 0.0;
                            description="ecoclimatic index of the future climate (0 to 1)"),
              StateVariable(:new, Int32(0); description="trees infected this year"),
              (StateVariable(coh[k], Int32(0);
                             description="trees infected $k year$(k == 1 ? "" : "s") ago, spreading")
               for k in 1:spreading_years)...,
              StateVariable(:dead, Int32(0); description="trees killed by the fungus"),
              StateVariable(:infected_year, Int32(0);
                            description="year of the cell's first infection (0 for none)")]
    parameters = [Parameter(:cell, Float64(cell); description="side of a cell (m)"),
                  Parameter(:trees_per_ha, 500.0;
                            description="most host trees per hectare at full coverage"),
                  Parameter(:spore_mean, 280.0; description="mean flight of spores (m)"),
                  Parameter(:spore_max, Float64(spore_max);
                            description="longest flight of spores (m)"),
                  Parameter(:dispatches, 100;
                            description="spore dispatches from each infected seedling lot"),
                  Parameter(:transport_max, 300_000.0;
                            description="longest transport of a seedling lot (m)"),
                  Parameter(:transport_tries, 10_000;
                            description="draws of a lot's destination before it is given up"),
                  Parameter(:insects_per_tree, 1.0;
                            description="infected insects leaving each spreading tree a year"),
                  Parameter(:insect_tries, 10;
                            description="flights an insect makes to find a tree to infect"),
                  Parameter(:insect_mean, 1000.0; description="mean flight of insects (m)"),
                  Parameter(:insect_max, 4500.0;
                            description="longest flight of insects (m); Inf for no limit"),
                  Parameter(:insect_distribution, :exponential;
                            description=":exponential or :weibull flight distances"),
                  Parameter(:insect_shape, 1.5;
                            description="shape of the Weibull flight distance"),
                  Parameter(:insect_method, :pooled;
                            description=":pooled (a multinomial draw per cell) or :each (every insect flown)"),
                  Parameter(:entries, NamedTuple[];
                            description="entry points (i, j, first, last, lots): lots a year from first to last"),
                  Parameter(:first_year, 2005;
                            description="calendar year of the first step"),
                  Parameter(:climate_year, 2100;
                            description="calendar year of the future climate")]
    return CAModel(rule; name=:pitch_canker, states, parameters,
                   description="Möykkynen et al. (2015) pitch canker in Europe")
end

# The cohorts of spreading trees of a state, i1, i2, ..., in order.
function _pc_cohorts(u)
    ks = [k for k in keys(u) if occursin(r"^i\d+$", String(k))]
    return sort(ks; by=k -> parse(Int, String(k)[2:end]))
end

# Infect one healthy tree, with probability (N_healthy / N_max) EI, in the cell holding the point (x, y)
# metres; points outside the lattice infect nothing. Returns whether a tree was infected.
function _pc_infect!(u, rng, x, y, p, w, t)
    c = p.cell
    i = floor(Int, x / c) + 1
    j = floor(Int, y / c) + 1
    checkbounds(Bool, u.hosts, i, j) || return false
    h = u.hosts[i, j]
    h > 0 || return false
    ei = u.ei_now[i, j] + w * (u.ei_future[i, j] - u.ei_now[i, j])
    rand(rng) < h / _pc_nmax(p) * ei || return false
    u.hosts[i, j] -= Int32(1)
    u.new[i, j] += Int32(1)
    u.infected_year[i, j] == 0 && (u.infected_year[i, j] = Int32(t))
    return true
end

# The point where a seedling lot from (x0, y0) is planted: a uniform distance up to transport_max in a
# uniform direction, accepted with probability the destination cell's host coverage and drawn again
# otherwise (points outside the lattice are rejected); nothing after transport_tries draws.
function _pc_transport(rng, coverage, x0, y0, p)
    c = p.cell
    for _ in 1:(p.transport_tries)
        r = p.transport_max * rand(rng)
        θ = 2π * rand(rng)
        x, y = x0 + r * cos(θ), y0 + r * sin(θ)
        i, j = floor(Int, x / c) + 1, floor(Int, y / c) + 1
        checkbounds(Bool, coverage, i, j) || continue
        rand(rng) < coverage[i, j] && return (x, y)
    end
    return nothing
end

"""
    pitch_canker_spread!(u, p, t, rng) -> u

The arrivals of infected seedling lots and the flights of infected insects of
[Moykkynen2015](@cite) in year `t`, applied to the state `u` in place after
the year's spores (see [`build_pitch_canker`](@ref)). Each infection takes
one healthy tree of the cell where it lands with probability `(N_healthy /
N_max) EI` (Eq. 1), lowering `N_healthy` at once; cells are squares of side
`cell` metres, and points outside the lattice infect nothing.

1. Every entry point `(i, j, first, last, lots)` in `entries` with `first ≤ t
   ≤ last` sends `lots` seedling lots from the midpoint of cell `(i, j)`, each
   a uniform distance up to `transport_max` in a uniform direction, accepted
   with probability the destination cell's host coverage and drawn again
   otherwise (given up after `transport_tries` draws). Each lot sends
   `dispatches` spore dispatches from where it is planted, each a negative
   exponential distance of mean `spore_mean` (drawn again above `spore_max`)
   in a uniform direction, which infect where they land.
2. The trees spreading this year (cohorts `i1` to the last) of each cell send
   `insects_per_tree` infected insects each (the total's whole part, and one
   more with probability its fractional part) from the cell's midpoint. Each
   insect flies a distance of mean `insect_mean`, negative exponential or
   Weibull of shape `insect_shape` (`insect_distribution`), drawn again above
   `insect_max`, in a uniform direction, and tries to infect where it lands;
   if it fails it flies again from the tree, `insect_tries` times in all.

Lots are taken in the order of `entries` and insects cell by cell in column
order, all drawing from `rng`. With `insect_method = :each` every insect is
flown as above. With `:pooled` (the default, far faster when many trees
spread) a cell's insects are drawn together: given the healthy trees when
they set out, an insect infects in the cell at offset `o`, on one of its
tries, with probability `K(o) q(o) (1 - (1 - s)^T) / s`, where `K(o)` is the
probability that a flight from the midpoint lands there, `q(o) = (N_healthy /
N_max) EI` there, `s = Σ K q` and `T = insect_tries`; the numbers infecting
each cell are a multinomial draw over those probabilities, at most the cell's
healthy trees. The two agree in distribution as long as one cell's insects
take few of any cell's healthy trees; without a maximum flight (`insect_max =
Inf`) the pooled draw loses the flights longer than the distance that fewer
than one in a million exceed.
"""
function pitch_canker_spread!(u, p, t::Integer, rng)
    w = _pc_climate_weight(p, t)
    c = p.cell
    for e in p.entries
        e.first <= t <= e.last || continue
        x0, y0 = (e.i - 0.5) * c, (e.j - 0.5) * c
        for _ in 1:(e.lots)
            dest = _pc_transport(rng, u.coverage, x0, y0, p)
            dest === nothing && continue
            for _ in 1:(p.dispatches)
                d = _pc_truncexp(rng, p.spore_mean, p.spore_max)
                θ = 2π * rand(rng)
                _pc_infect!(u, rng, dest[1] + d * cos(θ), dest[2] + d * sin(θ), p, w, t)
            end
        end
    end
    coh = _pc_cohorts(u)
    S = zeros(Int, size(u.hosts))
    for k in coh
        S .+= getfield(u, k)
    end
    p.insect_method === :pooled && return _pc_insects_pooled!(u, rng, S, p, w, t)
    p.insect_method === :each ||
        throw(ArgumentError("insect_method must be :pooled or :each, got :$(p.insect_method)"))
    for J in CartesianIndices(S)
        S[J] > 0 || continue
        insects = _pc_round(rng, S[J] * p.insects_per_tree)
        x0, y0 = (J[1] - 0.5) * c, (J[2] - 0.5) * c
        for _ in 1:insects
            for _ in 1:(p.insect_tries)
                d = _pc_flight(rng, p)
                θ = 2π * rand(rng)
                _pc_infect!(u, rng, x0 + d * cos(θ), y0 + d * sin(θ), p, w, t) && break
            end
        end
    end
    return u
end

"""
    pitch_canker_callback() -> DiscreteCallback

A callback that runs [`pitch_canker_spread!`](@ref) after each year, drawing
from a `Xoshiro` seeded from the problem's seed and the year.
"""
function pitch_canker_callback()
    affect! = function (integ)
        t = Int(integ.t)
        pitch_canker_spread!(integ.u, integ.p, t,
                             Random.Xoshiro(hash((integ.prob.seed, t, :pitch_canker))))
        return nothing
    end
    return DiscreteCallback((u, t, integ) -> true, affect!; save_positions=(false, false))
end

"""
    pitch_canker_landscape(rng; dims=(600, 200), forest=0.5, patchiness=0.3, coast=0.2,
                           inland=0.6) -> (; coverage, ei)

A synthetic coastal landscape of `dims` cells for
[`problem_pitch_canker`](@ref), the coast along the first column: host
`coverage` from a random cluster landscape (see
[`random_cluster_landscape`](@ref), with `patchiness` its share of seed
cells, drawn on a square and cropped), zero on the `1 - forest` share of
cells with the lowest values and rising from 0.05 to 0.95 on the rest; and an
ecoclimatic index `ei` that is 1 within `coast` of the width from the coast
and falls linearly to 0 at `inland` of the width, as the source's index is
highest on the coasts.
"""
function pitch_canker_landscape(rng; dims::Tuple{Integer,Integer}=(600, 200),
                                forest::Real=0.5, patchiness::Real=0.3, coast::Real=0.2,
                                inland::Real=0.6)
    m, n = dims
    min(m, n) >= 2 || throw(ArgumentError("dims must be at least 2 × 2, got $dims"))
    0 < forest <= 1 || throw(ArgumentError("forest must lie in (0, 1], got $forest"))
    0 <= coast < inland <= 1 ||
        throw(ArgumentError("need 0 ≤ coast < inland ≤ 1, got coast = $coast, inland = $inland"))
    v = random_cluster_landscape(rng, max(m, n), patchiness)[1:m, 1:n]
    cut = sort(vec(v))[clamp(round(Int, (1 - forest) * m * n), 1, m * n)]
    hi = maximum(v)
    coverage = [v[I] > cut || forest == 1 ?
                0.05 + 0.9 * (hi > cut ? (v[I] - cut) / (hi - cut) : 1.0) : 0.0
                for I in CartesianIndices(v)]
    ei = [clamp((inland - (j - 1) / (n - 1)) / (inland - coast), 0.0, 1.0)
          for _ in 1:m, j in 1:n]
    return (; coverage, ei)
end

"""
    pitch_canker_entries(m; years=20) -> Vector{NamedTuple}

Three entry points on the coast (column 2) of a landscape `m` cells long, as
the source's three (Table 2: Gijón and Bilbao, about 200 km apart, from year
1, and Brindisi from year 4, each 100 lots a year to year 20): at rows
`m/6` and `m/2` from year 1 and at row `5m/6` from year 4, each sending 100
lots a year until year `years`.
"""
function pitch_canker_entries(m::Integer; years::Integer=20)
    m >= 6 || throw(ArgumentError("m must be at least 6, got $m"))
    return [(; i=round(Int, m / 6), j=2, first=1, last=years, lots=100),
            (; i=round(Int, m / 2), j=2, first=1, last=years, lots=100),
            (; i=round(Int, 5m / 6), j=2, first=4, last=years, lots=100)]
end

"""
    problem_pitch_canker(; dims=(600, 200), years=20, spreading_years=5, coverage=nothing, ei=nothing,
                         ei_future=nothing, entries=nothing, seed=0, pmap=(;), kwargs...)
        -> CAProblem

The setup of [Moykkynen2015](@cite) on cells of 1 km: `years` yearly steps
from a landscape with no infection, the fungus arriving through `entries`
(see [`pitch_canker_spread!`](@ref); by default
[`pitch_canker_entries`](@ref) for the landscape's rows), with the source's parameters unless
`pmap` overrides them. The host `coverage` and the current index `ei` are
those of [`pitch_canker_landscape`](@ref) of size `dims` drawn from `seed`
unless given (as matrices of one size); `ei_future` is the index of the future climate (by default
the current one). Each cell starts with `round(coverage N_max)` healthy trees.
The lattice is a `Grid` with a constant zero boundary, so spores, lots and
insects leaving it are lost. `spreading_years` and `pmap.spore_max` and
`pmap.cell` set the model's cohorts and reach (see
[`build_pitch_canker`](@ref)); the callback is stored with the problem.
"""
function problem_pitch_canker(; dims::Tuple{Integer,Integer}=(600, 200), years::Integer=20,
                              spreading_years::Integer=5, coverage=nothing, ei=nothing,
                              ei_future=nothing, entries=nothing, seed::Integer=0,
                              pmap=(;), kwargs...)
    years >= 1 || throw(ArgumentError("years must be at least 1, got $years"))
    (coverage === nothing) == (ei === nothing) ||
        throw(ArgumentError("give both coverage and ei, or neither"))
    if coverage === nothing
        land = pitch_canker_landscape(Random.Xoshiro(hash((seed, 0, :pitch_canker))); dims)
        coverage, ei = land.coverage, land.ei
    end
    size(coverage) == size(ei) ||
        throw(ArgumentError("coverage and ei differ in size: $(size(coverage)) and $(size(ei))"))
    ei_future === nothing && (ei_future = ei)
    size(ei_future) == size(ei) ||
        throw(ArgumentError("ei_future is $(size(ei_future)), not $(size(ei))"))
    all(x -> 0 <= x <= 1, coverage) || throw(ArgumentError("coverage must lie in [0, 1]"))
    all(x -> 0 <= x <= 1, ei) && all(x -> 0 <= x <= 1, ei_future) ||
        throw(ArgumentError("ei and ei_future must lie in [0, 1]"))
    m1, m2 = size(coverage)
    entries === nothing && (entries = pitch_canker_entries(m1; years))
    for e in entries
        checkbounds(Bool, coverage, e.i, e.j) ||
            throw(ArgumentError("entry point ($(e.i), $(e.j)) lies outside the $m1 × $m2 landscape"))
    end
    given = merge((; entries=collect(entries)),
                  (; (Symbol(k) => v for (k, v) in pairs(pmap))...))
    model = build_pitch_canker(; spreading_years, spore_max=get(given, :spore_max, 1300.0),
                               cell=get(given, :cell, 1000.0))
    nmax = get(given, :trees_per_ha, 500.0) * get(given, :cell, 1000.0)^2 / 10_000
    z() = zeros(Int32, m1, m2)
    coh = ntuple(k -> Symbol(:i, k), spreading_years)
    u0 = NamedTuple{(:hosts, :coverage, :ei_now, :ei_future, :new, coh..., :dead,
                     :infected_year)}((Int32.(round.(coverage .* nmax)),
                                       Float64.(coverage), Float64.(ei),
                                       Float64.(ei_future), z(),
                                       (z() for _ in coh)..., z(), z()))
    return CAProblem(model, u0, years, given; space=Grid(m1, m2), seed,
                     callback=pitch_canker_callback(), kwargs...)
end

"""
    pitch_canker_infected(u) -> Matrix{Int}

The trees of each cell infected and still alive: this year's infections and
every cohort that spreads the fungus.
"""
function pitch_canker_infected(u)
    I = Int.(u.new)
    for k in _pc_cohorts(u)
        I .+= getfield(u, k)
    end
    return I
end

"""
    pitch_canker_infected_area(u) -> Float64

The source's "percentage of infected host area": the host area (coverage
times cell area) of the cells where any tree has been infected, alive or
dead, as a percentage of all the host area.
"""
function pitch_canker_infected_area(u)
    total = sum(u.coverage)
    total > 0 || return 0.0
    hit = pitch_canker_infected(u) .+ u.dead .> 0
    return 100 * sum(u.coverage[hit]) / total
end

_register!(:pitch_canker; build=build_pitch_canker, problem=problem_pitch_canker)
