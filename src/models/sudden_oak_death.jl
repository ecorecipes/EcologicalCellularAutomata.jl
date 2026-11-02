# Meentemeyer, Cunniffe, Cook, Filipe, Hunter, Rizzo & Gilligan (2011): sudden oak death (Phytophthora
# ramorum) spreading through California's forests in weekly steps, on cells of 250 m holding up to 100
# host units, susceptible or infected. Infected units release spores when the forest type's season and the
# weather allow; spores travel by a mixture of a short-range and a long-range Cauchy kernel, and infect
# susceptible units where they land, again as the weather allows. The source's statewide host and weather
# rasters were never published, so problems run on a synthetic landscape with synthetic weather.
#
# The spores landing in each cell are the sum, over all sources, of Poisson emissions thinned by the
# kernel, so they are Poisson with a mean that is a convolution: two kernel fields (short and long range)
# computed on a torus twice the landscape's size whose padding is masked, so spores leaving are lost.

const SOD_NONE = Int8(0)
const SOD_R = Int8(1)
const SOD_M = Int8(2)

"""
    sod_kernel(scale; radius, cell=250, samples=16) -> Matrix{Float64}

The probability that a spore released at the centre of a cell lands in each
cell within `radius` cells, for a Cauchy dispersal kernel of `scale` metres
on cells of side `cell` metres: the spore travels in a uniform direction a
distance from the half-Cauchy distribution `2 / (π scale (1 + (d/scale)²))`,
so that the kernel `(1 + (d/scale)²)^-1` of [Meentemeyer2011](@cite) is the
density per unit distance. The centre cell is integrated exactly; the others
by the midpoint rule on `samples × samples` points. Spores landing beyond the
radius are lost, so the probabilities sum to less than one.
"""
function sod_kernel(scale::Real; radius::Integer, cell::Real=250, samples::Integer=16)
    scale > 0 || throw(ArgumentError("scale must be positive, got $scale"))
    radius >= 1 || throw(ArgumentError("radius must be at least 1, got $radius"))
    g(r) = 1 / (π^2 * scale * r * (1 + (r / scale)^2))   # density per unit area
    P = zeros(2radius + 1, 2radius + 1)
    h = cell / samples
    for j in (-radius):radius, i in (-radius):radius
        (i, j) == (0, 0) && continue
        s = 0.0
        for b in 1:samples, a in 1:samples
            x = (i - 0.5) * cell + (a - 0.5) * h
            y = (j - 0.5) * cell + (b - 0.5) * h
            s += g(hypot(x, y))
        end
        P[i + radius + 1, j + radius + 1] = s * h^2
    end
    # the centre: in direction θ the cell's edge lies at (cell/2)/max(|cos θ|, |sin θ|)
    n = 4000
    P[radius + 1, radius + 1] = sum(1:n) do k
        θ = (k - 0.5) * (π / 2) / n
        return 2 / π * atan(cell / 2 / max(cos(θ), sin(θ)) / scale)
    end / n
    return P
end

# Whether forest of type `f` takes part in week `w` (1-52): redwood-tanoak January to September,
# mixed evergreen February to September.
function _sod_active(f, w)
    month = cld(12w, 52)
    return f == SOD_R ? month <= 9 : f == SOD_M ? 2 <= month <= 9 : false
end

# The week of the year (1-52) of step `t` (the step from t - 1 to t, starting from t = 1).
_sod_week(t) = mod1(t, 52)

# min(Poisson(λ), cap) by inversion of one uniform.
function _sod_poisson(u, λ, cap)
    λ <= 0 && return 0
    λ > 2cap + 50 && return cap
    k, q = 0, exp(-λ)
    F = q
    while u > F && k < cap
        k += 1
        q *= λ / k
        F += q
    end
    return k
end

"""
    build_sudden_oak_death(; short=sod_kernel(20.57; radius=4),
                           long=sod_kernel(9504; radius=99), method=:fft) -> CAModel

The sudden oak death model of [Meentemeyer2011](@cite), in weekly steps. Each
cell holds `S` susceptible and `I` infected host units (at most `nmax`
together), a forest type (`forest`: 0 none, 1 redwood-tanoak, 2 mixed
evergreen) and whether it releases spores this week (`sporulating`). In week
`t`, with weather index `w_t` (moisture times temperature suitability, from
`weather`), every infected unit in a sporulating cell releases a Poisson
number of spore units at rate `beta w_t`; each travels by the short-range
kernel `short` with probability `gamma` and by the long-range kernel `long`
otherwise (with the long-range share multiplied by `long_factor` from week
`management_week` on); where it lands, in a cell whose forest is in season,
it infects one of the `S` susceptible units with probability
`w_t S / nmax`. New infections in a cell are therefore Poisson, capped at
`S`. Redwood-tanoak forest takes part from January to September and mixed
evergreen from February to September. `method = :fft` computes the long-range
field by fast Fourier transform.
"""
function build_sudden_oak_death(; short::AbstractMatrix=sod_kernel(20.57; radius=4),
                                long::AbstractMatrix=sod_kernel(9504; radius=99),
                                method::Symbol=:fft)
    totals = (sum(short), sum(long))
    emit = c -> c.sporulating ? Float64(c.I) : 0.0
    fields = (; short=KernelField(EmpiricalKernel(short), emit),
              long=KernelField(EmpiricalKernel(long), emit; method))
    rule = LocalRule(Offsets((1, 0))) do c, nbrs, ctx
        p = ctx.p
        c.forest == SOD_NONE && return c
        t = ctx.t + 1
        week = _sod_week(t)
        wt = p.weather[mod1(t, length(p.weather))]
        S, I = c.S, c.I
        if S > 0 && wt > 0 && _sod_active(c.forest, week)
            longshare = (1 - p.gamma) * (t >= p.management_week ? p.long_factor : 1.0)
            spores = p.beta * wt *
                     (p.gamma * totals[1] * field(ctx, :short) +
                      longshare * totals[2] * field(ctx, :long))
            new = _sod_poisson(uniform(ctx, 1), spores * wt * S / p.nmax, Int(S))
            S -= Int8(new)
            I += Int8(new)
        end
        return (; S, I, forest=c.forest,
                sporulating=I > 0 && _sod_active(c.forest, _sod_week(t + 1)))
    end
    states = [StateVariable(:S, Int8(0); description="susceptible host units"),
              StateVariable(:I, Int8(0); description="infected host units"),
              StateVariable(:forest, SOD_NONE; domain=Int8(0):Int8(2),
                            description="0 no host, 1 redwood-tanoak, 2 mixed evergreen"),
              StateVariable(:sporulating, false;
                            description="whether the cell's infected units release spores this week")]
    parameters = [Parameter(:beta, 4.4;
                            description="spore units released per infected unit a week"),
                  Parameter(:gamma, 0.9947;
                            description="share of spore units dispersing at short range"),
                  Parameter(:nmax, 100; description="most host units in a cell"),
                  Parameter(:weather, [1.0];
                            description="weekly weather index (moisture times temperature suitability), from week 1"),
                  Parameter(:long_factor, 1.0;
                            description="factor on the long-range share from management_week on"),
                  Parameter(:management_week, typemax(Int);
                            description="week from which long-range dispersal is reduced")]
    return CAModel(rule; name=:sudden_oak_death, states, parameters, fields,
                   description="Meentemeyer et al. (2011) sudden oak death in California")
end

"""
    sod_temperature_index(x) -> Float64

The source's suitability of a week of mean temperature `x` °C for sporulation,
fitted to Davidson et al. (2005): `-0.066 + 0.056 x - 0.0036 (x - 15)² -
0.0003 (x - 15)³`, at least 0.
"""
function sod_temperature_index(x)
    return max(0.0,
               -0.066 + 0.056x - 0.0036 * (x - 15)^2 - 0.0003 * (x - 15)^3)
end

"""
    sod_weather(rng; years=19, wetness=1, mean_temperature=13, temperature_range=4)
        -> Matrix{Float64}

Synthetic weekly weather for a coastal Californian forest, `52 × years`: each
week's index is its moisture (the share of its 7 days with rain) times its
temperature suitability ([`sod_temperature_index`](@ref)). Temperatures
follow a seasonal cycle about `mean_temperature`, warmest in August, with
weekly noise; rain falls mostly from November to April, more in wet years
(each year's wetness is `wetness` times a lognormal factor).
"""
function sod_weather(rng; years::Integer=19, wetness::Real=1, mean_temperature::Real=13,
                     temperature_range::Real=4)
    W = zeros(52, years)
    for y in 1:years
        wet = wetness * exp(0.5 * randn(rng))
        for w in 1:52
            x = mean_temperature - temperature_range * cospi(2 * (w - 6) / 52) +
                1.5 * randn(rng)
            p = clamp(wet * (0.05 + 0.4 * max(0.0, cospi(2 * (w - 2) / 52))), 0.0, 1.0)
            days = count(_ -> rand(rng) < p, 1:7)
            W[w, y] = days / 7 * sod_temperature_index(x)
        end
    end
    return W
end

"""
    problem_sudden_oak_death(; n=100, weather=:random, years=1990:2030, observed=1990:2008,
                             sites=3, wetness=0.5, management=nothing, seed=0, pmap=(;),
                             kwargs...)
        -> CAProblem

The setup of [Meentemeyer2011](@cite) on a synthetic landscape of `n × n`
cells of 250 m: host units from a random cluster landscape (half the cells
hold none), redwood-tanoak forest on about a third of the host cells and
mixed evergreen on the rest, and `sites` host cells near the centre (within
`n / 8` cells of it, or the nearest host cells if there are too few there)
each with one infected unit at the start of `years`. The weather of the `observed` years is
drawn by [`sod_weather`](@ref) with `wetness`; later years repeat observed years at random
(`weather = :random`), from the better half by weather index
(`:favourable`) or from the worse half (`:unfavourable`), as in the source.
`management = (year, factor)` multiplies the long-range share of dispersal by
`factor` from that year on. All is drawn from `seed`. The landscape lies in
the corner of a torus twice its size whose other cells are masked, so spores
leaving it are lost; the state at time `t` is the end of week `t`.
"""
function problem_sudden_oak_death(; n::Integer=100, weather::Symbol=:random,
                                  years::AbstractUnitRange=1990:2030,
                                  observed::AbstractUnitRange=1990:2008, sites::Integer=3,
                                  wetness::Real=0.5, management=nothing, seed::Integer=0,
                                  pmap=(;), kwargs...)
    weather in (:random, :favourable, :unfavourable) ||
        throw(ArgumentError("weather must be :random, :favourable or :unfavourable, got :$weather"))
    first(observed) == first(years) && last(observed) <= last(years) ||
        throw(ArgumentError("the observed years must start the run, got $observed in $years"))
    n >= 10 || throw(ArgumentError("n must be at least 10, got $n"))
    rng = Random.Xoshiro(hash((seed, 0, :sudden_oak_death)))
    W = sod_weather(rng; years=length(observed), wetness)
    ranked = sortperm(vec(sum(W; dims=1)))
    half = length(ranked) ÷ 2
    pool = weather === :random ? collect(1:length(observed)) :
           weather === :favourable ? ranked[(end - half + 1):end] : ranked[1:half]
    future = [rand(rng, pool) for _ in 1:(length(years) - length(observed))]
    series = vec(hcat(W, W[:, future]))
    v = random_cluster_landscape(rng, n, 0.3)
    hosts = classify_landscape(v, [0.5, 0.5]) .== 1
    kind = classify_landscape(random_element_landscape(rng, n, 200), [1 / 3, 2 / 3])
    S = zeros(Int8, 2n, 2n)
    forest = zeros(Int8, 2n, 2n)
    for j in 1:n, i in 1:n
        hosts[i, j] || continue
        S[i, j] = clamp(round(Int, 100 * v[i, j]), 1, 100)
        forest[i, j] = kind[i, j] == 0 ? SOD_R : SOD_M
    end
    away(I) = maximum(abs.(Tuple(I) .- (n + 1) / 2))
    host = [I for I in CartesianIndices((n, n)) if forest[I] != SOD_NONE]
    length(host) >= sites ||
        throw(ArgumentError("the landscape has fewer than $sites host cells"))
    reach = max(n / 8, sort(away.(host))[sites])   # widened to the nearest hosts if needed
    centre = [I for I in host if away(I) <= reach]
    I0 = zeros(Int8, 2n, 2n)
    for c in centre[Random.randperm(rng, length(centre))[1:sites]]
        S[c] -= Int8(1)
        I0[c] = Int8(1)
    end
    sporulating = [I0[c] > 0 && _sod_active(forest[c], 1) for c in CartesianIndices(I0)]
    mgmt = management === nothing ? (;) :
           (; management_week=52 * (management[1] - first(years)) + 1,
            long_factor=Float64(management[2]))
    given = merge((; weather=series), mgmt,
                  (; (Symbol(k) => v for (k, v) in pairs(pmap))...))
    model = build_sudden_oak_death(; long=sod_kernel(9504; radius=n - 1))
    mask = fill(false, 2n, 2n)
    mask[1:n, 1:n] .= true
    return CAProblem(model, (; S, I=I0, forest, sporulating), 52 * length(years), given;
                     space=Torus(2n, 2n), seed, mask, kwargs...)
end

"""
    sod_infected_area(u; n=size(u.I, 1) ÷ 2) -> Float64

The source's infected area, in km²: the cells of the `n × n` landscape with
any infected unit, at 6.25 ha each.
"""
function sod_infected_area(u; n::Integer=size(u.I, 1) ÷ 2)
    return 0.0625 * count(>(0), view(u.I, 1:n, 1:n))
end

_register!(:sudden_oak_death; build=build_sudden_oak_death,
           problem=problem_sudden_oak_death)
