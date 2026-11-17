# Adrakey, Gibson, Eveillard, Malembic-Maher & Fabre (2023): flavescence dorée in three vineyard plots at
# Faleyras, a yearly SIR with roguing. A healthy vine is infected in a year with probability
# 1 - exp(-(epsilon + beta sum_j c_j K(r_ij; alpha) / alpha^2)), summing over the vines infectious that year
# (c_j a factor for Cabernet Sauvignon sources); an infected vine is cryptic for the rest of the year,
# symptomatic and infectious the next, and rogued the winter after; healthy vines also disappear for other
# reasons with probability q a year. The kernel is not normalised. Vines stand on a fine raster, and the
# sum over sources is a kernel field (the core normalises kernel weights to one, so the rule multiplies
# by the normalising constant).

const FD_NONE = Int8(0)
const FD_S = Int8(1)
const FD_E = Int8(2)
const FD_I = Int8(3)
const FD_R = Int8(4)
const FD_KERNELS = (:gaussian, :exponential, :inverse_square, :power4, :power6, :power8)

"""
    fd_kernel(shape, r, α) -> Float64

The dispersal kernels of [Adrakey2023](@cite), unnormalised, at distance `r`
(m) with scale `α` (m): `:gaussian` `exp(-r²/α²)`, `:exponential`
`exp(-r/α)`, `:inverse_square` `(α/r)²`, and the power laws `:power4`,
`:power6`, `:power8`, `(1 + r/α)^-b` with `b` 4, 6 or 8.
"""
function fd_kernel(shape::Symbol, r::Real, α::Real)
    shape === :gaussian && return exp(-(r / α)^2)
    shape === :exponential && return exp(-r / α)
    shape === :inverse_square && return (α / r)^2
    shape === :power4 && return (1 + r / α)^-4
    shape === :power6 && return (1 + r / α)^-6
    shape === :power8 && return (1 + r / α)^-8
    return throw(ArgumentError("kernel must be one of $FD_KERNELS, got :$shape"))
end

"""
    build_flavescence_doree(; kernel=:power6, α=5.0, cultivar_effect=1.0, radius=60.0,
                            resolution=0.5) -> CAModel

The yearly SIR with roguing of [Adrakey2023](@cite) on a raster of
`resolution` m. Layers: `vine` (0 no vine, 1 healthy, 2 infected and
cryptic, 3 symptomatic and infectious, 4 removed) and `cultivar` (1 Merlot, 2
Cabernet Sauvignon). In each year's step a healthy vine is removed for other
reasons with probability `q`, and otherwise infected with probability
`1 - exp(-(ε + β Z field))`; an infected vine becomes symptomatic, and a
symptomatic one is rogued. A step is a season: a vine infected in the season
before becomes symptomatic, and infects others, in this one (the vector acquires
and transmits within a season), and is rogued in the winter after it. So the
field `pressure` sums, over the vines infected the season before (state 2,
symptomatic this season) within `radius` m, the `kernel`'s weight ([`fd_kernel`](@ref) with scale
`α`), times `cultivar_effect` for a Cabernet Sauvignon source (the source's
cultivar effect is 1.42); `Z` (a parameter, set by
[`problem_flavescence_doree`](@ref)) is the sum of the kernel's weights over
the raster divided by `α²`, which undoes the core's normalisation, so that the
force of infection is `ε + β Σ c_j K(r_ij; α) / α²` as in the source. `α` and
the kernel are fixed when the model is built, since the field's weights depend
on them.
"""
function build_flavescence_doree(; kernel::Symbol=:power6, α::Real=5.0,
                                 cultivar_effect::Real=1.0, radius::Real=60.0,
                                 resolution::Real=0.5)
    kernel in FD_KERNELS ||
        throw(ArgumentError("kernel must be one of $FD_KERNELS, got :$kernel"))
    rule = ProbabilisticRule() do c, nbrs, ctx
        p = ctx.p
        if c.vine == FD_S
            qi = -expm1(-(p.ε + p.β * p.Z * field(ctx, :pressure)))
            return ((vine=FD_S, cultivar=c.cultivar) => (1 - p.q) * (1 - qi),
                    (vine=FD_E, cultivar=c.cultivar) => (1 - p.q) * qi,
                    (vine=FD_R, cultivar=c.cultivar) => p.q)
        elseif c.vine == FD_E
            return ((vine=FD_I, cultivar=c.cultivar) => 1.0,)
        elseif c.vine == FD_I
            return ((vine=FD_R, cultivar=c.cultivar) => 1.0,)
        end
        return (c => 1.0,)
    end
    w = _fd_weights(kernel, α, radius, resolution)
    source = c -> c.vine == FD_E ? (c.cultivar == 2 ? Float64(cultivar_effect) : 1.0) : 0.0
    states = [StateVariable(:vine, FD_NONE; domain=Int8(0):Int8(4),
                            description="0 none, 1 healthy, 2 infected, 3 symptomatic, 4 removed"),
              StateVariable(:cultivar, Int8(0);
                            description="1 Merlot, 2 Cabernet Sauvignon")]
    parameters = [Parameter(:ε, 0.0012; description="primary infection rate (per year)"),
                  Parameter(:β, 1.0; description="secondary infection rate"),
                  Parameter(:q, 0.005;
                            description="yearly probability a healthy vine is removed for other reasons"),
                  Parameter(:Z, sum(w) / α^2;
                            description="normalising constant of the kernel field over α²"),
                  Parameter(:first_year, 2013; description="the year of time 0"),
                  Parameter(:resolution, Float64(resolution);
                            description="raster cell (m)"),
                  Parameter(:kernel, kernel;
                            description="the kernel's shape (fixed at build)"),
                  Parameter(:radius, Float64(radius);
                            description="the kernel's radius (m)"),
                  Parameter(:cultivar_effect, Float64(cultivar_effect);
                            description="the factor of a Cabernet Sauvignon source"),
                  Parameter(:cells, CartesianIndex{2}[];
                            description="the raster cell of each vine of fd_vines")]
    return CAModel(rule; name=:flavescence_doree, states, parameters,
                   fields=(pressure=KernelField(EmpiricalKernel(w), source;
                                                include_centre=false, method=:fft),),
                   description="Adrakey et al. (2023) flavescence dorée, $kernel kernel")
end

function _fd_weights(kernel, α, radius, h)
    R = ceil(Int, radius / h)
    w = zeros(2R + 1, 2R + 1)
    for i in (-R):R, j in (-R):R
        r = h * hypot(i, j)
        (r == 0 || r > radius) && continue
        w[i + R + 1, j + R + 1] = fd_kernel(kernel, r, α)
    end
    return w
end

"""
    fd_vines() -> NamedTuple

The vines of [Adrakey2023](@cite) (Etalab 2.0; see
`models/flavescence_doree/data/NOTICE.md`): `x`, `y` (m), `plot` (`"CS"`,
`"M1"`, `"M2"`), `cultivar` (`:merlot`, `:cabernet_sauvignon`) and `status`
(`:absent_2018`, `:symptomatic_2018`, `:absent_2019`, `:symptomatic_2019`,
`:healthy`).
"""
function fd_vines()
    path = joinpath(MODELS_DIR, "flavescence_doree", "data", "vines.csv")
    rows = [split(l, ',') for l in Iterators.drop(eachline(path), 1)]
    return (x=[parse(Float64, r[1]) for r in rows], y=[parse(Float64, r[2]) for r in rows],
            plot=[String(r[3]) for r in rows], cultivar=[Symbol(r[4]) for r in rows],
            status=[Symbol(r[5]) for r in rows])
end

"""
    problem_flavescence_doree(; first_year=2013, last_year=2019, kernel=:power6, α=5.0,
                              cultivar_effect=1.0, radius=60.0, resolution=0.5, seed=0,
                              pmap=(;)) -> CAProblem

The three plots of [Adrakey2023](@cite) on a periodic raster of `resolution` m
with an empty margin of `radius` m (so that nothing wraps around), every vine
healthy in `first_year`, one step a year to `last_year`; time `t` is year
`first_year + t`. `vines` (in the problem's parameters) maps each vine of
[`fd_vines`](@ref) to its raster cell.
"""
function problem_flavescence_doree(; first_year::Integer=2013, last_year::Integer=2019,
                                   kernel::Symbol=:power6, α::Real=5.0,
                                   cultivar_effect::Real=1.0, radius::Real=60.0,
                                   resolution::Real=0.5, seed::Integer=0, pmap=(;))
    last_year > first_year || throw(ArgumentError("last_year must follow first_year"))
    v = fd_vines()
    h = resolution
    pad = ceil(Int, radius / h) + 1
    ix = [round(Int, x / h) + pad + 1 for x in v.x]
    iy = [round(Int, y / h) + pad + 1 for y in v.y]
    m, n = maximum(ix) + pad, maximum(iy) + pad
    allunique(zip(ix, iy)) ||
        throw(ArgumentError("two vines share a cell at resolution $h m; use a finer one"))
    vine = fill(FD_NONE, m, n)
    cultivar = zeros(Int8, m, n)
    for k in eachindex(ix)
        vine[ix[k], iy[k]] = FD_S
        cultivar[ix[k], iy[k]] = v.cultivar[k] === :cabernet_sauvignon ? Int8(2) : Int8(1)
    end
    model = build_flavescence_doree(; kernel, α, cultivar_effect, radius, resolution)
    p = merge(defaults(model), (; first_year, cells=CartesianIndex.(ix, iy)),
              (; (Symbol(k) => x for (k, x) in pairs(pmap))...))
    return CAProblem(model, (; vine, cultivar), last_year - first_year, p;
                     space=Torus(m, n),
                     seed)
end

"""
    fd_symptomatic(u, cells) -> Vector{Bool}

Whether each vine (at its raster cell in `cells`) is symptomatic in the state
`u`: what a survey that year records.
"""
fd_symptomatic(u, cells) = [u.vine[c] == FD_I for c in cells]

_register!(:flavescence_doree; build=build_flavescence_doree,
           problem=problem_flavescence_doree)

# ------------------------------------------------------------------ fitting

"""
    FDHistory

A history of the vines for the yearly model, by the years of the steps in which
things happen: `infected[i]` (the vine is infected in the step from that year,
symptomatic two years on and rogued after; 0 for never), `removed[i]` (removed
for another reason in the step from that year; 0 for never), `open_end[i]`
(the last step is summed over: the vine was healthy at the last survey, so
susceptible or newly infected), and `observed_to[i]` (the last year whose step
is observed; later steps are censored, as for a plot uprooted after a survey).
"""
mutable struct FDHistory
    infected::Vector{Int}
    removed::Vector{Int}
    open_end::Vector{Bool}
    observed_to::Vector{Int}
end
function Base.copy(z::FDHistory)
    return FDHistory(copy(z.infected), copy(z.removed), copy(z.open_end),
                     copy(z.observed_to))
end

"""
    FDFit

What [`fd_loglik`](@ref) needs: the vines' raster `cells` and the raster's
`resolution`, the `kernel` and its `radius`, each vine's `factor` as a source
(the cultivar effect), the `first_year` and `last_year` of the steps, which vines
can be sources (`source`), for each of those its neighbours within the radius
with their squared distances in cells (`pairs`), and the latent vines'
`candidates`. Build one with [`fd_fit`](@ref).
"""
struct FDFit
    cells::Vector{CartesianIndex{2}}
    resolution::Float64
    kernel::Symbol
    radius::Float64
    factor::Vector{Float64}
    first_year::Int
    last_year::Int
    source::Vector{Bool}
    pairs::Vector{Vector{Tuple{Int,Int32}}}
    candidates::Vector{Vector{Tuple{Int,Int}}}
end

"""
    fd_fit(cells, cultivar; resolution=0.5, kernel=:power6, radius=60.0, cultivar_effect=1.0,
           first_year, last_year, source=trues(length(cells)), candidates=[]) -> FDFit

The data structure for fitting: vines at raster `cells` with `cultivar` (1
Merlot, 2 Cabernet Sauvignon); `source` marks the vines that can ever be
infectious (others' neighbour lists are skipped); `candidates[i]`, if not
empty, lists the histories `(infected, removed)` the latent step may give vine
`i`. [`fd_survey_fit`](@ref) builds it for the survey.
"""
function fd_fit(cells, cultivar; resolution::Real=0.5, kernel::Symbol=:power6,
                radius::Real=60.0, cultivar_effect::Real=1.0, first_year::Integer,
                last_year::Integer, source=trues(length(cells)),
                candidates=[Tuple{Int,Int}[] for _ in cells])
    kernel in FD_KERNELS ||
        throw(ArgumentError("kernel must be one of $FD_KERNELS, got :$kernel"))
    n = length(cells)
    pairs = [Tuple{Int,Int32}[] for _ in 1:n]
    for j in 1:n
        source[j] || continue
        for i in 1:n
            i == j && continue
            d2 = sum(abs2, Tuple(cells[i] - cells[j]))
            resolution * sqrt(d2) <= radius && push!(pairs[j], (i, Int32(d2)))
        end
    end
    factor = [c == 2 ? Float64(cultivar_effect) : 1.0 for c in cultivar]
    return FDFit(collect(cells), Float64(resolution), kernel, Float64(radius), factor,
                 Int(first_year), Int(last_year), collect(Bool, source), pairs,
                 collect(candidates))
end

# The kernel over α² at every squared raster distance within the radius: distances on the raster
# are resolution * sqrt(d2), so a table indexed by d2 replaces a power per pair.
function _fd_table(θ, fit::FDFit)
    m = floor(Int, (fit.radius / fit.resolution)^2)
    a2 = θ.α^2
    return [d2 == 0 ? 0.0 : fd_kernel(fit.kernel, fit.resolution * sqrt(d2), θ.α) / a2
            for d2 in 0:m]
end

# The infection pressure on every vine in the step from every year from the vines symptomatic in
# that step's season (infected in the step before): S[i, y] for years first_year:last_year - 1,
# without ε and β.
function _fd_pressure(θ, z::FDHistory, fit::FDFit, tab=_fd_table(θ, fit))
    Y0, Y1 = fit.first_year, fit.last_year
    S = zeros(length(fit.cells), Y1 - Y0)
    for j in eachindex(z.infected)
        y = z.infected[j]
        y == 0 && continue
        y += 1
        Y0 <= y <= Y1 - 1 || continue
        _fd_add!(S, j, y - Y0 + 1, 1.0, fit, tab)
    end
    return S
end

function _fd_add!(S, j, col, sign, fit, tab)
    f = sign * fit.factor[j]
    @inbounds for (i, d2) in fit.pairs[j]
        S[i, col] += f * tab[d2 + 1]
    end
    return S
end

# Vine i's term for the step from year y given its pressure s in that year (0 if it is not then
# healthy or the step is censored).
@inline function _fd_step(θ, z, fit, i, y, s)
    inf, rem = z.infected[i], z.removed[i]
    (y > z.observed_to[i] || (inf != 0 && y > inf) || (rem != 0 && y > rem)) && return 0.0
    y == rem && return log(θ.q)
    l1q = log1p(-θ.q)
    y == inf && return l1q + log(-expm1(-(θ.ε + θ.β * s)))
    (z.open_end[i] && y == fit.last_year - 1) && return l1q
    return l1q - (θ.ε + θ.β * s)
end

function _fd_terms(θ, z, fit, S, i, ys)
    ll = 0.0
    for y in ys
        ll += _fd_step(θ, z, fit, i, y, S[i, y - fit.first_year + 1])
    end
    return ll
end

"""
    fd_loglik(θ, z::FDHistory, fit::FDFit) -> Float64

The log-probability of the history `z` under the yearly model with parameters
`θ` (`ε`, `β`, `q`, `α`): over every vine and every observed step in which it is
healthy, `log q` if it is then removed for another reason, otherwise
`log(1 - q)` plus `-λ` (it stays healthy) or `log(1 - exp(-λ))` (it is infected),
with `λ = ε + β Σ c_j K(r_ij; α) / α²` over the vines symptomatic, and
infectious, in that season (infected in the step before); a vine healthy at the last survey contributes
`log(1 - q)` for the last step. With every step observed and no open ends it is
the core's `path_logpdf` of the corresponding path.
"""
function fd_loglik(θ, z::FDHistory, fit::FDFit)
    S = _fd_pressure(θ, z, fit)
    ys = fit.first_year:(fit.last_year - 1)
    return sum(i -> _fd_terms(θ, z, fit, S, i, ys), eachindex(z.infected); init=0.0)
end

"""
    fd_path(z::FDHistory, prob) -> Vector

The states of `prob` (vines at `prob.p.cells`) in each year from `first_year`
under the history `z`, every vine healthy at the start.
"""
function fd_path(z::FDHistory, prob)
    Y0 = prob.p.first_year
    return map(0:(prob.tspan[2])) do t
        y = Y0 + t
        vine = copy(prob.u0.vine)
        for (i, c) in enumerate(prob.p.cells)
            inf, rem = z.infected[i], z.removed[i]
            vine[c] = rem != 0 && y > rem ? FD_R :
                      inf == 0 || y <= inf ? FD_S :
                      y == inf + 1 ? FD_E : y == inf + 2 ? FD_I : FD_R
        end
        return (vine=vine, cultivar=prob.u0.cultivar)
    end
end

"""
    fd_latent_step!(rng, z::FDHistory, θ, fit::FDFit) -> Float64

One sweep of Metropolis-Hastings updates of the vines with `candidates`: each
proposes one of its candidate histories uniformly, accepted with the change in
[`fd_loglik`](@ref), computed from the vine's own terms and its neighbours' in
the years whose pressure changes. Returns the share of moves accepted.
"""
function fd_latent_step!(rng, z::FDHistory, θ, fit::FDFit)
    tab = _fd_table(θ, fit)
    S = _fd_pressure(θ, z, fit, tab)
    ys = fit.first_year:(fit.last_year - 1)
    tried = accepted = 0
    for k in eachindex(fit.candidates)
        cands = fit.candidates[k]
        isempty(cands) && continue
        tried += 1
        old = (z.infected[k], z.removed[k])
        new = cands[rand(rng, eachindex(cands))]
        new == old && (accepted += 1; continue)
        Δ = _fd_move!(z, S, k, old, new, θ, fit, ys, tab)
        if log(rand(rng)) < Δ
            accepted += 1
        else
            _fd_move!(z, S, k, new, old, θ, fit, ys, tab)
        end
    end
    return tried == 0 ? 1.0 : accepted / tried
end

# Move vine k's history from `old` to `new`, updating the pressures, and return the change in
# the log-likelihood: k's own terms, and its neighbours' terms in the (at most two) steps in which
# k is infectious under the old or the new history.
function _fd_move!(z, S, k, old, new, θ, fit, ys, tab=_fd_table(θ, fit))
    Y0 = fit.first_year
    own_before = _fd_terms(θ, z, fit, S, k, ys)
    z.infected[k], z.removed[k] = new
    Δ = _fd_terms(θ, z, fit, S, k, ys) - own_before
    yo = old[1] == 0 ? 0 : old[1] + 1
    yn = new[1] == 0 ? 0 : new[1] + 1
    yo == yn && return Δ
    f = fit.factor[k]
    for (y, sign) in ((yo, -1.0), (yn, 1.0))
        (y != 0 && Y0 <= y <= fit.last_year - 1) || continue
        col = y - Y0 + 1
        @inbounds for (i, d2) in fit.pairs[k]
            s0 = S[i, col]
            s1 = s0 + sign * f * tab[d2 + 1]
            Δ += _fd_step(θ, z, fit, i, y, s1) - _fd_step(θ, z, fit, i, y, s0)
            S[i, col] = s1
        end
    end
    return Δ
end

"""
    fd_survey_fit(prob) -> (fit, history)

The fitting data for the survey of [Adrakey2023](@cite) on the raster, kernel,
radius and cultivar effect of `prob` ([`problem_flavescence_doree`](@ref))
(steps from `prob.p.first_year` to 2019): vines symptomatic in 2018 or 2019 were
infected two years before; vines absent from a Merlot plot in 2019 were removed
for other reasons in 2018; the Cabernet Sauvignon plot was uprooted after the
2018 survey, so its vines' last step is censored; vines healthy in 2019 have an
open last step; and each vine absent in 2018 is latent, either rogued after
showing symptoms in 2015, 2016 or 2017, or removed for another reason in a year
before 2018 (its candidates). The history starts each latent vine as removed
for another reason in the first year.
"""
function fd_survey_fit(prob)
    v = fd_vines()
    Y0 = prob.p.first_year
    Y0 <= 2013 || throw(ArgumentError("the survey needs steps from 2013 or earlier"))
    n = length(v.x)
    cultivar = [c === :cabernet_sauvignon ? 2 : 1 for c in v.cultivar]
    inf, rem = zeros(Int, n), zeros(Int, n)
    open_end, observed_to = falses(n), fill(2018, n)
    candidates = [Tuple{Int,Int}[] for _ in 1:n]
    for i in 1:n
        s = v.status[i]
        if s === :symptomatic_2018
            inf[i] = 2016
        elseif s === :symptomatic_2019
            inf[i] = 2017
        elseif s === :absent_2019
            v.plot[i] == "CS" ? (observed_to[i] = 2017) : (rem[i] = 2018)
        elseif s === :healthy
            open_end[i] = true
        else
            candidates[i] = [[(y, 0) for y in (2013, 2014, 2015)];
                             [(0, y) for y in Y0:2017]]
            rem[i] = Y0
        end
    end
    source = (v.status .=== :symptomatic_2018) .| (v.status .=== :absent_2018)
    fit = fd_fit(prob.p.cells, cultivar; resolution=prob.p.resolution, kernel=prob.p.kernel,
                 radius=prob.p.radius, cultivar_effect=prob.p.cultivar_effect,
                 first_year=Y0,
                 last_year=2019, source, candidates)
    return fit, FDHistory(inf, rem, open_end, observed_to)
end
