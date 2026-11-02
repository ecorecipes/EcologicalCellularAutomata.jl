# Kleczkowski, Gilligan & Bailey (1997): primary and secondary infection of radish by Rhizoctonia
# solani on a 5 x 10 lattice of plants. A susceptible plant is infected by the inoculum units in its
# 20 mm square (primary infection, rate r_p(x, t) per unit at distance x) and by its infectious von
# Neumann neighbours (secondary infection, rate r_s(t)); both rates decay with absolute time. Latent
# plants become infectious at rate g. Rates add as independent hazards, and a step of length dt
# turns a total rate R into the probability 1 - exp(-R dt).

const RHIZOCTONIA_S = Int8(0)
const RHIZOCTONIA_L = Int8(1)
const RHIZOCTONIA_I = Int8(2)

"""
    build_rhizoctonia_contact(; secondary=:per_neighbour) -> CAModel

The primary and secondary infection model of [Kleczkowski1997](@cite), for
damping-off of radish by *Rhizoctonia solani* in microcosms. Each cell is a
plant, susceptible (0), latent (1) or infectious (2). A step is `dt` days of
absolute time `τ_p + t * dt`. In each step, synchronously:

- a susceptible plant becomes latent with probability `1 - exp(-R dt)`, where
  `R` is the sum of the primary rate `C_p (1 + η_p x) exp(-σ_p x) exp(-d_p t)`
  over the inoculum units in its square, at distances `x`, and the secondary
  rate `C_s exp(-σ_s x_d²) exp(-d_s t)` times the number of infectious
  neighbours, with `C_p = a_p d_p exp(d_p τ_p)` and `C_s = a_s d_s exp(d_s τ_s)`;
- a latent plant becomes infectious with probability `1 - exp(-g dt)`.

The neighbours are the von Neumann neighbours, with none beyond the edges.
The primary part of `R` reads the covariate `pathozone`, the sum of
`(1 + η_p x) exp(-σ_p x)` over the units in each plant's square, which
[`rhizoctonia_pathozone`](@ref) computes. With `secondary = :any` the secondary
rate is counted once if any neighbour is infectious, the literal reading of
the source's rule, instead of once per infectious neighbour.
"""
function build_rhizoctonia_contact(; secondary::Symbol=:per_neighbour)
    secondary in (:per_neighbour, :any) ||
        throw(ArgumentError("secondary must be :per_neighbour or :any, got :$secondary"))
    per_neighbour = secondary === :per_neighbour
    rule = LocalRule(VonNeumann(1)) do cell, nbrs, ctx
        p = ctx.p
        t = p.τ_p + ctx.t * p.dt
        if cell == RHIZOCTONIA_S
            Cp = p.a_p * p.d_p * exp(p.d_p * p.τ_p)
            Cs = p.a_s * p.d_s * exp(p.d_s * p.τ_s)
            k = count(==(RHIZOCTONIA_I), nbrs)
            k = per_neighbour ? k : min(k, 1)
            rate = Cp * exp(-p.d_p * t) * covariate(ctx, :pathozone) +
                   k * Cs * exp(-p.σ_s * p.x_d^2) * exp(-p.d_s * t)
            return bernoulli(ctx, -expm1(-rate * p.dt), 1) ? RHIZOCTONIA_L : RHIZOCTONIA_S
        elseif cell == RHIZOCTONIA_L
            return bernoulli(ctx, -expm1(-p.dt / p.latent_period), 2) ? RHIZOCTONIA_I :
                   RHIZOCTONIA_L
        else
            return cell
        end
    end
    parameters = [Parameter(:a_p, 1.09; description="amplitude of primary infection"),
                  Parameter(:σ_p, 0.277;
                            description="spatial decay of primary infection (1/mm)"),
                  Parameter(:η_p, 0.49;
                            description="critical-exponential correction (1/mm)"),
                  Parameter(:d_p, 0.55;
                            description="temporal decay of primary infection (1/d)"),
                  Parameter(:τ_p, 4.898;
                            description="time offset of primary infection (d)"),
                  Parameter(:a_s, 14.9; description="amplitude of secondary infection"),
                  Parameter(:σ_s, 0.00099;
                            description="spatial decay of secondary infection (1/mm²)"),
                  Parameter(:d_s, 0.36;
                            description="temporal decay of secondary infection (1/d)"),
                  Parameter(:τ_s, 3.73;
                            description="time offset of secondary infection (d)"),
                  Parameter(:latent_period, 3.89; description="mean latent period 1/g (d)"),
                  Parameter(:x_d, 20.0; description="distance between plants (mm)"),
                  Parameter(:dt, 0.1; description="length of a step (d)")]
    return CAModel(rule; name=:rhizoctonia_contact,
                   states=[StateVariable(:plant, RHIZOCTONIA_S; domain=Int8(0):Int8(2),
                                         description="0 susceptible, 1 latent, 2 infectious")],
                   parameters, covariates=(:pathozone,),
                   description="Kleczkowski et al. (1997) Rhizoctonia solani primary and secondary infection, secondary $secondary")
end

"""
    rhizoctonia_inoculum(rng; units=10, rows=5, cols=10, spacing=20.0) -> Vector{NTuple{2,Float64}}

`units` inoculum positions, in mm, drawn uniformly at random from the
`rows * spacing` by `cols * spacing` rectangle that holds the plants.
"""
function rhizoctonia_inoculum(rng; units::Integer=10, rows::Integer=5, cols::Integer=10,
                              spacing::Real=20.0)
    return [(rows * spacing * rand(rng), cols * spacing * rand(rng)) for _ in 1:units]
end

"""
    rhizoctonia_pathozone(positions; rows=5, cols=10, spacing=20.0, σ_p=0.277, η_p=0.49)
        -> Matrix{Float64}

The covariate `pathozone` of the Rhizoctonia model: at each plant, the sum of
`(1 + η_p x) exp(-σ_p x)` over the inoculum units in its square, at distances
`x` from it. Plant `(i, j)` stands at `((i - 1/2) spacing, (j - 1/2) spacing)`,
so the plants' squares tile the rectangle and every unit lies in one square.
"""
function rhizoctonia_pathozone(positions; rows::Integer=5, cols::Integer=10,
                               spacing::Real=20.0, σ_p::Real=0.277, η_p::Real=0.49)
    z = zeros(rows, cols)
    for (x, y) in positions
        (0 <= x <= rows * spacing && 0 <= y <= cols * spacing) ||
            throw(ArgumentError("inoculum at ($x, $y) mm lies outside the plants' rectangle"))
        i = clamp(floor(Int, x / spacing) + 1, 1, rows)
        j = clamp(floor(Int, y / spacing) + 1, 1, cols)
        d = hypot(x - (i - 0.5) * spacing, y - (j - 0.5) * spacing)
        z[i, j] += (1 + η_p * d) * exp(-σ_p * d)
    end
    return z
end

"""
    problem_rhizoctonia_contact(; units=10, rows=5, cols=10, secondary=:per_neighbour,
                                tspan=152, seed=0, pmap=(;)) -> CAProblem

The setup of [Kleczkowski1997](@cite): a 5 × 10 lattice of susceptible plants
20 mm apart, and 10 inoculum units at uniformly random positions drawn from
`seed`, from day `τ_p` for 152 steps of 0.1 d (to day 20.1).
"""
function problem_rhizoctonia_contact(; units::Integer=10, rows::Integer=5,
                                     cols::Integer=10, secondary::Symbol=:per_neighbour,
                                     tspan=152, seed::Integer=0, pmap=(;))
    model = build_rhizoctonia_contact(; secondary)
    p = merge(defaults(model), (; (Symbol(k) => v for (k, v) in pairs(pmap))...))
    positions = rhizoctonia_inoculum(Random.Xoshiro(seed); units, rows, cols,
                                     spacing=p.x_d)
    pathozone = rhizoctonia_pathozone(positions; rows, cols, spacing=p.x_d, p.σ_p, p.η_p)
    return CAProblem(model, (;), tspan, pmap; space=Grid(rows, cols), seed,
                     covariates=(; pathozone))
end

"""
    rhizoctonia_days(sol; τ_p=4.898, dt=0.1) -> Vector{Float64}

The days since planting of the saved states of a solution of the Rhizoctonia
model, whose steps are `dt` days long from day `τ_p`.
"""
rhizoctonia_days(sol; τ_p::Real=4.898, dt::Real=0.1) = τ_p .+ dt .* sol.t

"""
    rhizoctonia_infected(u) -> Int

The number of infected plants, latent or infectious, in a state of the
Rhizoctonia model: what the source's experiments count, since latent and
infectious plants both show damping-off.
"""
rhizoctonia_infected(u) = count(!=(RHIZOCTONIA_S), u)

_register!(:rhizoctonia_contact; build=build_rhizoctonia_contact,
           problem=problem_rhizoctonia_contact)
