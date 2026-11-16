# Cunniffe et al. (2014): Bahia bark scaling of citrus in the experimental grove of Laranjeira et
# al. (2006). A susceptible tree is infected at a rate beta times the sum, over symptomatic trees, of
# an exponential kernel normalised in two dimensions, (2 pi alpha^2)^-1 exp(-d / alpha); infected
# trees become symptomatic, and infectious, at rate rho; young trees can be neither infected nor
# infectious before a maturity delay delta. The trees are the vertices of a graph, each linked to
# every other, and each step is a month, so the continuous-time model of the source becomes a
# discrete-time one in which a rate R over a month gives the probability 1 - exp(-R).

const CITRUS_S = Int8(0)
const CITRUS_E = Int8(1)
const CITRUS_I = Int8(2)
# The survey's last month, and a month standing for "not within the survey".
const CITRUS_MONTHS = 82
const CITRUS_NEVER = typemax(Int32)

"""
    build_citrus_bark_scaling() -> CAModel

The spatial SEI model of Bahia bark scaling of citrus of [Cunniffe2014](@cite),
fitted there to the grove of [Laranjeira2006](@cite). Each cell is a tree on a
`GraphSpace` linking every tree with every other: susceptible (0),
infected but not yet symptomatic (1), or symptomatic and infectious (2). A
step is a month; in the step from month `t`:

- a susceptible tree is infected with probability `1 - exp(-λ a(t))`, where
  `λ = β Σ exp(-d / α) / (2π α²)` sums over the symptomatic trees at distances
  `d` (the covariate `distance` holds each tree's distances to its neighbours, in
  the order of the graph), and `a(t) = clamp(t + 1 - δ, 0, 1)` is the share of
  the month after the trees reach epidemiological maturity at month `δ`;
- an infected tree becomes symptomatic with probability `1 - exp(-ρ)`.

It is a `ProbabilisticRule`, so `path_logpdf` gives the exact
likelihood of a path observed every month, and [`citrus_loglik`](@ref) gives
the same quickly from the months of infection and symptoms.
"""
function build_citrus_bark_scaling()
    rule = ProbabilisticRule() do c, nbrs, ctx
        p = ctx.p
        if c == CITRUS_S
            a = clamp(ctx.t + 1 - p.δ, 0.0, 1.0)
            d = covariate(ctx, :distance)
            s = 0.0
            for k in eachindex(nbrs)
                nbrs[k] == CITRUS_I && (s += exp(-d[k] / p.α))
            end
            q = -expm1(-p.β * s / (2π * p.α^2) * a)
            return (CITRUS_S => 1 - q, CITRUS_E => q, CITRUS_I => 0.0)
        elseif c == CITRUS_E
            q = -expm1(-p.ρ)
            return (CITRUS_S => 0.0, CITRUS_E => 1 - q, CITRUS_I => q)
        end
        return (CITRUS_S => 0.0, CITRUS_E => 0.0, CITRUS_I => 1.0)
    end
    parameters = [Parameter(:α, 2.51; description="dispersal scale of the kernel (m)"),
                  Parameter(:β, 4.42; description="rate of infection (m² per month)"),
                  Parameter(:ρ, 0.170;
                            description="rate at which infected trees show symptoms (per month)"),
                  Parameter(:δ, 22.7;
                            description="months from planting to epidemiological maturity")]
    return CAModel(rule; name=:citrus_bark_scaling,
                   states=[StateVariable(:tree, CITRUS_S; domain=Int8(0):Int8(2),
                                         description="0 susceptible, 1 infected, 2 symptomatic")],
                   parameters, covariates=(:distance,),
                   description="Cunniffe et al. (2014) Bahia bark scaling of citrus")
end

"""
    citrus_layout(; sources=:figure) -> NamedTuple

The positions, in metres, of the trees in the grove of [Laranjeira2006](@cite):
`grove`, a 240 × 2 matrix of the young trees (16 rows of 15 at 2 m, row 1 at
`y = 0` nearest the source, tree `(r, c)` at row `15(r - 1) + c`), and
`sources`, the diseased adult trees at 6 m × 4 m, the nearest 5 m from row 1.
`sources = :figure` gives the 30 trees (6 rows of 5) drawn in Figure 2a of
[Cunniffe2014](@cite); `:text` the 25 (5 rows of 5) of its text.
"""
function citrus_layout(; sources::Symbol=:figure)
    sources in (:figure, :text) ||
        throw(ArgumentError("sources must be :figure or :text, got :$sources"))
    grove = [Float64[2(c - 1), 2(r - 1)] for r in 1:16 for c in 1:15]
    ys = sources === :figure ? (-5.0, -9.0, -13.0, -17.0, -21.0, -25.0) :
         (-5.0, -9.0, -13.0, -17.0, -21.0)
    src = [Float64[x, y] for y in ys for x in (2.0, 8.0, 14.0, 20.0, 26.0)]
    return (grove=permutedims(reduce(hcat, grove)), sources=permutedims(reduce(hcat, src)))
end

"""
    citrus_survey() -> NamedTuple

The survey of [Laranjeira2006](@cite), as published with [Cunniffe2014](@cite)
(CC BY 4.0; see `models/citrus_bark_scaling/data/NOTICE.md`): `rows` and `cols`
of the 240 trees, `first_symptomatic` (the month since planting of the
inspection at which each tree was first seen with symptoms, `missing` for the
three never seen with them), and the `inspections` (months since planting),
with their `dates` and the numbers `newly_symptomatic`.
"""
function citrus_survey()
    dir = joinpath(MODELS_DIR, "citrus_bark_scaling", "data")
    trees = [split(l, ',') for l in Iterators.drop(eachline(joinpath(dir, "trees.csv")), 1)]
    insp = [split(l, ',')
            for l in Iterators.drop(eachline(joinpath(dir, "inspections.csv")), 1)]
    first = Union{Missing,Int}[isempty(t[3]) ? missing : parse(Int, t[3]) for t in trees]
    return (rows=[parse(Int, t[1]) for t in trees], cols=[parse(Int, t[2]) for t in trees],
            first_symptomatic=first, inspections=[parse(Int, r[2]) for r in insp],
            dates=[String(r[1]) for r in insp],
            newly_symptomatic=[parse(Int, r[3]) for r in insp])
end

function _citrus_graph(positions)
    n = size(positions, 1)
    neighbours = [[j for j in 1:n if j != i] for i in 1:n]
    distance = [[hypot(positions[i, 1] - positions[j, 1],
                       positions[i, 2] - positions[j, 2])
                 for j in neighbours[i]] for i in 1:n]
    return GraphSpace(neighbours), distance
end

"""
    problem_citrus_bark_scaling(; sources=:figure, tspan=82, seed=0, pmap=(;)) -> CAProblem

The grove of [Laranjeira2006](@cite) from planting (month 0): 240 susceptible
young trees and the symptomatic adult trees of [`citrus_layout`](@ref), the
first 240 vertices being the grove, for `tspan` months (82, the last survey).
"""
function problem_citrus_bark_scaling(; sources::Symbol=:figure, tspan=CITRUS_MONTHS,
                                     seed::Integer=0, pmap=(;))
    layout = citrus_layout(; sources)
    positions = vcat(layout.grove, layout.sources)
    space, distance = _citrus_graph(positions)
    u0 = [fill(CITRUS_S, size(layout.grove, 1)); fill(CITRUS_I, size(layout.sources, 1))]
    return CAProblem(build_citrus_bark_scaling(), (tree=u0,), tspan, pmap; space, seed,
                     covariates=(; distance))
end

"""
    citrus_symptomatic(u) -> Int

The number of symptomatic young trees in a state (the first 240 vertices).
"""
citrus_symptomatic(u) = count(==(CITRUS_I), view(u, 1:240))

# ------------------------------------------------------------------ fitting

"""
    CitrusFit

What fitting the citrus model to the survey needs, from [`citrus_fit`](@ref):
the `distances` from each young tree to every tree (young trees first, then the
sources), the number of `sources`, and for each young tree the months between
which it first showed symptoms: after `lower` and by `upper` (`upper` is
`CITRUS_NEVER` for a tree never seen with symptoms, whose `lower` is the last
survey).
"""
struct CitrusFit
    distances::Matrix{Float64}
    sources::Int
    lower::Vector{Int}
    upper::Vector{Int}
    months::Int
end

"""
    citrus_fit(; sources=:figure) -> CitrusFit

The data for fitting the citrus model to the survey of [Laranjeira2006](@cite).
A tree first seen with symptoms at an inspection became symptomatic after the
previous inspection (or after planting, for the first) and by that one.
"""
function citrus_fit(; sources::Symbol=:figure)
    layout = citrus_layout(; sources)
    all = vcat(layout.grove, layout.sources)
    n = size(layout.grove, 1)
    D = [hypot(all[i, 1] - all[j, 1], all[i, 2] - all[j, 2])
         for i in 1:n, j in 1:size(all, 1)]
    s = citrus_survey()
    insp = s.inspections
    lower = Int[]
    upper = Int[]
    for m in s.first_symptomatic
        if ismissing(m)
            push!(lower, last(insp))
            push!(upper, CITRUS_NEVER)
        else
            k = findfirst(==(m), insp)
            push!(lower, k == 1 ? 0 : insp[k - 1])
            push!(upper, m)
        end
    end
    return CitrusFit(D, size(layout.sources, 1), lower, upper, last(insp))
end

"""
    CitrusLatent

The unobserved history of the young trees: for each, the month it was first
`infected` and the month it was first `symptomatic` (both `CITRUS_NEVER` if not
within the survey). Its state at month `t` is susceptible before `infected`,
infected from then until `symptomatic`, and symptomatic afterwards.
"""
mutable struct CitrusLatent
    infected::Vector{Int}
    symptomatic::Vector{Int}
end
Base.copy(z::CitrusLatent) = CitrusLatent(copy(z.infected), copy(z.symptomatic))

"""
    citrus_latent(fit::CitrusFit; δ=22.7, lag=4) -> CitrusLatent

A history consistent with the survey to start a sampler from: each tree seen
with symptoms became symptomatic at the inspection where it was seen and was
infected `lag` months earlier, but not before month `δ`; trees never seen with
symptoms stay susceptible.
"""
function citrus_latent(fit::CitrusFit; δ::Real=22.7, lag::Integer=4)
    tmin = floor(Int, δ) + 1
    sym = copy(fit.upper)
    inf = [u == CITRUS_NEVER ? CITRUS_NEVER : max(u - lag, tmin) for u in sym]
    all(i -> inf[i] == CITRUS_NEVER || inf[i] < sym[i], eachindex(inf)) ||
        throw(ArgumentError("δ = $δ leaves no month for infection before the first symptoms"))
    return CitrusLatent(inf, sym)
end

_citrus_maturity(t, δ) = clamp(t + 1 - δ, 0.0, 1.0)

# Cumulative maturity A[m + 1] = sum of a(t) for t = 0, ..., m - 1.
function _citrus_cumulative(δ, T)
    A = zeros(T + 1)
    for m in 1:T
        A[m + 1] = A[m] + _citrus_maturity(m - 1, δ)
    end
    return A
end

function _citrus_consistent(z::CitrusLatent, fit::CitrusFit, i)
    tE, tI = z.infected[i], z.symptomatic[i]
    if fit.upper[i] == CITRUS_NEVER
        tI == CITRUS_NEVER || return false
        return tE == CITRUS_NEVER || 1 <= tE <= fit.months
    end
    return fit.lower[i] < tI <= fit.upper[i] && 1 <= tE < tI
end

"""
    citrus_loglik(θ, z::CitrusLatent, fit::CitrusFit) -> Float64

The log-probability of the history `z` of the young trees under the citrus
model with parameters `θ` (`α`, `β`, `ρ`, `δ`), month by month to the last
survey, and `-Inf` if `z` disagrees with the survey: the complete-data
likelihood for data augmentation. It equals `path_logpdf` of the problem on
the path [`citrus_path`](@ref) builds, computed from the months of infection
and symptoms in `O(n²)` operations for `n` trees instead of `O(n² T)`.
"""
function citrus_loglik(θ, z::CitrusLatent, fit::CitrusFit)
    n = length(z.infected)
    all(i -> _citrus_consistent(z, fit, i), 1:n) || return -Inf
    T = fit.months
    A = _citrus_cumulative(θ.δ, T)
    c = θ.β / (2π * θ.α^2)
    stay, move = -θ.ρ, log(-expm1(-θ.ρ))
    ll = 0.0
    for i in 1:n
        tE, tI = z.infected[i], z.symptomatic[i]
        m = tE == CITRUS_NEVER ? T : tE - 1
        src = 0.0
        for s in 1:(fit.sources)
            src += exp(-fit.distances[i, n + s] / θ.α)
        end
        H = src * A[m + 1]
        λ = src
        for j in 1:n
            j == i && continue
            tj = z.symptomatic[j]
            tj < m || continue
            k = exp(-fit.distances[i, j] / θ.α)
            H += k * (A[m + 1] - A[tj + 1])
            λ += k
        end
        ll -= c * H
        tE == CITRUS_NEVER && continue
        # tree j counts in the infection step m when symptomatic by month m: tj <= m
        for j in 1:n
            j != i && z.symptomatic[j] == m && (λ += exp(-fit.distances[i, j] / θ.α))
        end
        ll += log(-expm1(-c * λ * _citrus_maturity(m, θ.δ)))
        ll += tI == CITRUS_NEVER ? (T - tE) * stay : (tI - 1 - tE) * stay + move
    end
    return ll
end

"""
    citrus_path(z::CitrusLatent, prob) -> Vector

The states of `prob` (from [`problem_citrus_bark_scaling`](@ref)) at months
`0, 1, ...` to the end of its time span under the history `z` of the young
trees, the sources symptomatic throughout: the path whose `path_logpdf` is
[`citrus_loglik`](@ref).
"""
function citrus_path(z::CitrusLatent, prob)
    n = length(z.infected)
    u0 = prob.u0
    return map(prob.tspan[1]:prob.tspan[2]) do t
        u = copy(u0)
        for i in 1:n
            u[i] = t >= z.symptomatic[i] ? CITRUS_I :
                   t >= z.infected[i] ? CITRUS_E : CITRUS_S
        end
        return u
    end
end

"""
    citrus_latent_step!(rng, z::CitrusLatent, θ, fit::CitrusFit) -> Float64

One sweep of Metropolis-Hastings updates of the history `z`, tree by tree,
leaving the distribution of the history given `θ` and the survey invariant; it
returns the share of moves accepted. For a tree seen with symptoms, the month
of symptoms is proposed uniformly between its inspections and the month of
infection uniformly from maturity to the month before; for a tree never seen
with symptoms, susceptible or infected at a uniform month, each with
probability one half. Each move's change in [`citrus_loglik`](@ref) is
computed from the terms that involve the tree, in `O(n)` operations.
"""
function citrus_latent_step!(rng, z::CitrusLatent, θ, fit::CitrusFit)
    n = length(z.infected)
    T = fit.months
    cache = _citrus_cache(z, θ, fit)
    tmin = floor(Int, θ.δ) + 1
    accepted = 0
    for k in 1:n
        tE, tI = z.infected[k], z.symptomatic[k]
        if fit.upper[k] == CITRUS_NEVER
            tE′ = (rand(rng) < 0.5 || tmin > T) ? CITRUS_NEVER : rand(rng, tmin:T)
            tI′ = CITRUS_NEVER
            logq = 0.0
        else
            tI′ = rand(rng, (fit.lower[k] + 1):(fit.upper[k]))
            tI′ - 1 >= tmin || continue
            tE′ = rand(rng, tmin:(tI′ - 1))
            logq = log(tI′ - tmin) - log(tI - tmin)
        end
        if tE′ == tE && tI′ == tI
            accepted += 1
            continue
        end
        Δ, lamk′ = _citrus_move_delta(z, k, tE′, tI′, cache, θ, T)
        if log(rand(rng)) < Δ + logq
            _citrus_accept!(z, k, tE′, tI′, lamk′, cache)
            accepted += 1
        end
    end
    return accepted / n
end

# Kernel weights, source pressure, cumulative maturity and each infected tree's hazard in its
# month of infection, for the current parameters and history.
function _citrus_cache(z, θ, fit)
    n = length(z.infected)
    K = [exp(-fit.distances[i, j] / θ.α) for i in 1:n, j in 1:n]
    src = [sum(s -> exp(-fit.distances[i, n + s] / θ.α), 1:(fit.sources); init=0.0)
           for i in 1:n]
    lam = [z.infected[i] == CITRUS_NEVER ? 0.0 :
           _citrus_hazard(z, i, z.infected[i] - 1, K, src) for i in 1:n]
    return (K=K, src=src, A=_citrus_cumulative(θ.δ, fit.months),
            c=θ.β / (2π * θ.α^2), lam=lam)
end

# The change in citrus_loglik when tree k's history becomes (tE′, tI′): its own terms, and its
# part in every other tree's susceptible period and month of infection.
function _citrus_move_delta(z, k, tE′, tI′, cache, θ, T)
    (; K, src, A, c, lam) = cache
    tE, tI = z.infected[k], z.symptomatic[k]
    own′, lamk′ = _citrus_own(z, k, tE′, tI′, K, src, A, c, θ, T)
    own, _ = _citrus_own(z, k, tE, tI, K, src, A, c, θ, T)
    Δ = own′ - own
    for i in eachindex(z.infected)
        i == k && continue
        tEi = z.infected[i]
        mi = tEi == CITRUS_NEVER ? T : tEi - 1
        Δ -= c * K[i, k] * (_citrus_overlap(A, tI′, mi) - _citrus_overlap(A, tI, mi))
        if tEi != CITRUS_NEVER
            din = (tI′ <= mi) - (tI <= mi)
            if din != 0
                a = _citrus_maturity(mi, θ.δ)
                Δ += log(-expm1(-c * (lam[i] + din * K[i, k]) * a)) -
                     log(-expm1(-c * lam[i] * a))
            end
        end
    end
    return Δ, lamk′
end

function _citrus_accept!(z, k, tE′, tI′, lamk′, cache)
    tI = z.symptomatic[k]
    for i in eachindex(z.infected)
        i == k && continue
        tEi = z.infected[i]
        tEi == CITRUS_NEVER && continue
        cache.lam[i] += ((tI′ <= tEi - 1) - (tI <= tEi - 1)) * cache.K[i, k]
    end
    z.infected[k] = tE′
    z.symptomatic[k] = tI′
    cache.lam[k] = lamk′
    return z
end

_citrus_overlap(A, t, m) = t < m ? A[m + 1] - A[t + 1] : 0.0

# Hazard on tree i (before the factor c and maturity) in the step from month m.
function _citrus_hazard(z, i, m, K, src)
    λ = src[i]
    for j in eachindex(z.symptomatic)
        j != i && z.symptomatic[j] <= m && (λ += K[i, j])
    end
    return λ
end

# Tree k's own terms of the log-likelihood with history (tEk, tIk), and its hazard in the month
# of infection (zero if not infected).
function _citrus_own(z, k, tEk, tIk, K, src, A, c, θ, T)
    n = length(z.infected)
    m = tEk == CITRUS_NEVER ? T : tEk - 1
    H = src[k] * A[m + 1]
    for j in 1:n
        j == k && continue
        H += K[k, j] * _citrus_overlap(A, z.symptomatic[j], m)
    end
    total = -c * H
    tEk == CITRUS_NEVER && return total, 0.0
    λ = _citrus_hazard(z, k, m, K, src)
    total += log(-expm1(-c * λ * _citrus_maturity(m, θ.δ)))
    total += tIk == CITRUS_NEVER ? (T - tEk) * -θ.ρ :
             (tIk - 1 - tEk) * -θ.ρ + log(-expm1(-θ.ρ))
    return total, λ
end

_register!(:citrus_bark_scaling; build=build_citrus_bark_scaling,
           problem=problem_citrus_bark_scaling)
