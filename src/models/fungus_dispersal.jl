# Ferreira, Moral, Ferreira & Godoy (2013): fungus dispersal among permissive and resistant hosts.
# A healthy host is infected with probability F (1 - prod_i (1 - f(d_i))^c) over every infectious
# host i, with f(d) = exp(-d / alpha) / (2 pi alpha^2). The product is a kernel field: with weights
# w(d) = -log(1 - f(d)) normalised by their total W, the field is sum_i w(d_i) / W, and the
# probability is F (1 - exp(-c W field)).

const FUNGUS_HEALTHY = Int8(0)
const FUNGUS_LATENT = Int8(1)
const FUNGUS_INFECTIOUS = Int8(2)

"""
    fungus_kernel(range; radius=59) -> (weights, total)

The dispersal weights of the fungus model: a `(2radius+1)²` matrix with
`-log(1 - f(d))` at distance `d` from the centre, where
`f(d) = exp(-d / range) / (2π range²)` is the probability that one propagule
lands on a host at distance `d`, and their sum `total`.
"""
function fungus_kernel(range::Real; radius::Integer=59)
    range > 0 || throw(ArgumentError("the propagule range must be positive, got $range"))
    radius >= 1 || throw(ArgumentError("the kernel radius must be at least 1, got $radius"))
    f(d) = exp(-d / range) / (2π * range^2)
    w = [-log1p(-f(hypot(i, j))) for i in (-radius):radius, j in (-radius):radius]
    return w, sum(w)
end

"""
    build_fungus_dispersal(; range=3.0, radius=59, seasonality=:none) -> CAModel

The fungus dispersal model of [Ferreira2013](@cite). Each cell is a host,
permissive or resistant (the layer `resistant`), and healthy, latent or
infectious (the layer `host`: 0, 1, 2), with `clock` the days since infection.
In each daily step, synchronously:

- a healthy host becomes latent with probability
  `F (1 - ∏ᵢ (1 - f(dᵢ))^c)`, the product over every infectious host `i` at
  distance `dᵢ`, where `f(d) = exp(-d / range) / (2π range²)`, `F` is `F_permissive`
  or `F_resistant`, and `c` is the daily propagule load of each infectious host;
- a latent host becomes infectious `latency` days after it was infected;
- an infectious host stays infectious.

`range` is the mean distance a propagule travels, in cells, and `radius` is
the farthest distance the kernel reaches; the default reaches across the
source's 60 × 60 lattice, so dispersal is global, as in the source. Both are
fixed when the model is built. `seasonality` sets the daily load: `:none`, a
constant `load`; `:temperate`, `load_max * exp(-(day - peak)² / width²)`;
`:tropical`, `load_max * amplitude * (sin(2π day / period) + 1)`, where `day = t + 1`.
"""
function build_fungus_dispersal(; range::Real=3.0, radius::Integer=59,
                                seasonality::Symbol=:none)
    seasonality in (:none, :temperate, :tropical) ||
        throw(ArgumentError("seasonality must be :none, :temperate or :tropical, got :$seasonality"))
    w, total = fungus_kernel(range; radius)
    infectious_field = KernelField(EmpiricalKernel(w), c -> c.host == FUNGUS_INFECTIOUS)
    loadf = if seasonality === :none
        (p, day) -> p.load
    elseif seasonality === :temperate
        (p, day) -> p.load_max * exp(-(day - p.peak)^2 / p.width^2)
    else
        (p, day) -> p.load_max * p.amplitude * (sin(2π * day / p.period) + 1)
    end
    rule = LocalRule(VonNeumann(1)) do cell, nbrs, ctx
        p = ctx.p
        if cell.host == FUNGUS_HEALTHY
            c = loadf(p, ctx.t + 1)
            F = cell.resistant ? p.F_resistant : p.F_permissive
            prob = F * -expm1(-c * total * field(ctx, :infectious))
            bernoulli(ctx, prob, 1) || return cell
            # infected now: latent with one day elapsed at the next step
            return (host=p.latency <= 1 ? FUNGUS_INFECTIOUS : FUNGUS_LATENT,
                    resistant=cell.resistant, clock=Int16(1))
        elseif cell.host == FUNGUS_LATENT
            clock = cell.clock + Int16(1)
            return (host=clock >= p.latency ? FUNGUS_INFECTIOUS : FUNGUS_LATENT,
                    resistant=cell.resistant, clock)
        else
            return cell
        end
    end
    states = [StateVariable(:host, FUNGUS_HEALTHY; domain=Int8(0):Int8(2),
                            description="0 healthy, 1 latent, 2 infectious"),
              StateVariable(:resistant, false;
                            description="resistant (true) or permissive"),
              StateVariable(:clock, Int16(0); description="days since infection")]
    parameters = [Parameter(:F_permissive, 0.003;
                            description="permissiveness of a permissive host"),
                  Parameter(:F_resistant, 0.001;
                            description="permissiveness of a resistant host"),
                  Parameter(:latency, 10; description="latent period (days)"),
                  Parameter(:load, 40.0;
                            description="daily propagule load (seasonality :none)"),
                  Parameter(:load_max, 270.0; description="peak daily load (seasonal)"),
                  Parameter(:peak, 60.0; description="day of peak load (:temperate)"),
                  Parameter(:width, 10.0;
                            description="width of the season in days (:temperate)"),
                  Parameter(:amplitude, 0.15; description="relative amplitude (:tropical)"),
                  Parameter(:period, 20.0; description="period in days (:tropical)")]
    return CAModel(rule; name=:fungus_dispersal, states, parameters,
                   fields=(; infectious=infectious_field),
                   description="Ferreira et al. (2013) fungus dispersal, range $range, seasonality $seasonality")
end

"""
    fungus_hosts(n; layout=:permissive, fraction=0.5, band_width=15, rng) -> Matrix{Bool}

The resistant hosts of an `n × n` field: `:permissive` (none), `:resistant`
(all), `:random` (each host resistant with probability `fraction`, drawn
from `rng`) or `:bands` (alternating vertical bands `band_width` columns wide,
permissive first).
"""
function fungus_hosts(n::Integer; layout::Symbol=:permissive, fraction::Real=0.5,
                      band_width::Integer=15, rng=Random.Xoshiro(0))
    layout === :permissive && return fill(false, n, n)
    layout === :resistant && return fill(true, n, n)
    layout === :random && return rand(rng, n, n) .< fraction
    layout === :bands &&
        return [isodd(fld(j - 1, band_width)) for i in 1:n, j in 1:n]
    return throw(ArgumentError("layout must be :permissive, :resistant, :random or :bands, got :$layout"))
end

"""
    problem_fungus_dispersal(; n=60, layout=:permissive, fraction=0.5, band_width=15,
                             range=3.0, radius=n-1, seasonality=:none, tspan=120,
                             seed=0, pmap=(;)) -> CAProblem

The setup of [Ferreira2013](@cite): an `n × n` field with no wrap-around, the
hosts laid out by [`fungus_hosts`](@ref) (a `:random` layout is drawn from
`seed`), nine infectious hosts in a 3 × 3 block at the centre, and 120 daily
steps.
"""
function problem_fungus_dispersal(; n::Integer=60, layout::Symbol=:permissive,
                                  fraction::Real=0.5, band_width::Integer=15,
                                  range::Real=3.0, radius::Integer=n - 1,
                                  seasonality::Symbol=:none, tspan=120, seed::Integer=0,
                                  pmap=(;))
    n >= 3 || throw(ArgumentError("the field must be at least 3 × 3, got $n"))
    model = build_fungus_dispersal(; range, radius, seasonality)
    resistant = fungus_hosts(n; layout, fraction, band_width, rng=Random.Xoshiro(seed))
    host = fill(FUNGUS_HEALTHY, n, n)
    c = n ÷ 2
    host[c:(c + 2), c:(c + 2)] .= FUNGUS_INFECTIOUS
    return CAProblem(model, (; host, resistant), tspan, pmap; space=Grid(n, n), seed)
end

"""
    fungus_infected(u) -> Int

The number of infected hosts, latent or infectious, in a state of the fungus
model.
"""
fungus_infected(u) = count(>=(FUNGUS_LATENT), u.host)

_register!(:fungus_dispersal; build=build_fungus_dispersal,
           problem=problem_fungus_dispersal)
