# Martins et al. (2000): citrus variegated chlorosis in orange groves. Each tree is healthy or
# diseased and stressed or not; each diseased tree sends n_v(t) vector flights a month to trees at
# symmetric Levy-distributed distances; a targeted healthy tree is infected unless stressed; an
# untargeted healthy tree redraws its stress; one sinusoid drives both the flights and the stress.
# The source picks targets from each diseased tree ("push"); here each healthy tree computes the
# probability of being targeted by any diseased tree ("pull"), 1 - prod_x (1 - K(y - x))^n_v, a
# kernel field with weights -log(1 - K). The kernel reaches across the grove, so the field is
# computed by FFT on a torus twice the grove's size whose other three quarters are masked.

const CVC_GROVES = (SJ01=(p_inf=0.0002, p_min=0.001), SJ67=(p_inf=0.04, p_min=0.0),
                    SJ71=(p_inf=0.06, p_min=0.35), SJ75=(p_inf=0.009, p_min=0.0))

"""
    cvc_flight_kernel(; stability=0.68, location=4.0, radius=199, samples=10^7, seed=1)
        -> Matrix{Float64}

The probability that one vector flight of the citrus model of
[Martins2000](@cite) lands on each offset within `radius`: the flight length
is `|X|` with `X` a symmetric stable variate of index `stability`, location
`location` and scale 1, the direction is uniform, and the landing point is
rounded to the nearest tree. Estimated from `samples` flights drawn with the
Chambers–Mallows–Stuck method from a fixed `seed`; flights beyond `radius` are
lost.
"""
function cvc_flight_kernel(; stability::Real=0.68, location::Real=4.0, radius::Integer=199,
                           samples::Integer=10^7, seed::Integer=1)
    0 < stability < 2 ||
        throw(ArgumentError("stability must lie in (0, 2), got $stability"))
    stability == 1 &&
        throw(ArgumentError("stability 1 (the Cauchy law) needs another sampler; use 0.999"))
    rng = Random.Xoshiro(seed)
    a = Float64(stability)
    K = zeros(2radius + 1, 2radius + 1)
    for _ in 1:samples
        V = π * (rand(rng) - 0.5)
        W = -log(rand(rng))
        X = location + sin(a * V) / cos(V)^(1 / a) * (cos(V - a * V) / W)^((1 - a) / a)
        r = abs(X)
        θ = 2π * rand(rng)
        i, j = round(Int, r * cos(θ)), round(Int, r * sin(θ))
        (abs(i) <= radius && abs(j) <= radius) || continue
        K[i + radius + 1, j + radius + 1] += 1
    end
    return K ./ samples
end

"""
    build_cvc_citrus(; kernel=cvc_flight_kernel(), method=:fft) -> CAModel

The citrus variegated chlorosis model of [Martins2000](@cite), in monthly
steps. Each tree is healthy or diseased (`diseased`) and stressed or not
(`stressed`; a diseased tree is stressed). With `n_v(t) = n0 + trunc(n0 sin(2πt/T + φ))`
flights from each diseased tree landing on offset `o` with probability
`kernel[o]` ([`cvc_flight_kernel`](@ref)), a healthy tree is targeted with
probability `1 - ∏ₓ (1 - K(y - x))^n_v(t)` over the diseased trees `x`. A
targeted tree is infected unless it is stressed, in which case it stays as it
was; an untargeted healthy tree is stressed at the next step with probability
`1 - p_nst(t)`, where `p_nst(t) = p_min + (1 - 2p_min)(1 + sin(2πt/T + φ)) / 2`.
Diseased trees stay diseased.

The kernel is fixed when the model is built. `method = :fft` computes the
targeting field by fast Fourier transform, which needs a periodic space:
[`problem_cvc_citrus`](@ref) embeds the grove in a torus twice its size and
masks the rest, so that nothing lies beyond the grove, as in the source.
"""
function build_cvc_citrus(; kernel::AbstractMatrix=cvc_flight_kernel(),
                          method::Symbol=:fft)
    all(<(1), kernel) || throw(ArgumentError("kernel probabilities must be below 1"))
    w = -log1p.(-kernel)
    total = sum(w)
    total > 0 || throw(ArgumentError("the kernel has no positive probability"))
    targeting = KernelField(EmpiricalKernel(w), c -> c.diseased; method)
    rule = LocalRule(VonNeumann(1)) do cell, nbrs, ctx
        cell.diseased && return cell
        p = ctx.p
        phase = 2π * ctx.t / p.period + p.phase
        nv = p.n0 + trunc(p.n0 * sin(phase))
        if bernoulli(ctx, -expm1(-nv * total * field(ctx, :targeting)), 1)
            return cell.stressed ? cell : (diseased=true, stressed=true)
        end
        pnst = p.p_min + (1 - 2p.p_min) * (1 + sin(phase)) / 2
        return (diseased=false, stressed=bernoulli(ctx, 1 - pnst, 2))
    end
    return CAModel(rule; name=:cvc_citrus,
                   states=[StateVariable(:diseased, false;
                                         description="diseased (true) or healthy"),
                           StateVariable(:stressed, false;
                                         description="under hydric or nutritional stress")],
                   parameters=[Parameter(:n0, 5.0;
                                         description="mean number of vector flights per diseased tree per month"),
                               Parameter(:p_min, 0.001;
                                         description="minimum probability that a tree is not stressed"),
                               Parameter(:period, 12.0;
                                         description="seasonal period (months)"),
                               Parameter(:phase, -π / 3;
                                         description="seasonal phase (radians)")],
                   fields=(; targeting),
                   description="Martins et al. (2000) citrus variegated chlorosis")
end

"""
    problem_cvc_citrus(; grove=:SJ01, L=200, kernel=nothing, tspan=20, seed=0, pmap=(;),
                       kwargs...) -> CAProblem

The setup of [Martins2000](@cite) for one of its four groves (`:SJ01`, `:SJ67`,
`:SJ71`, `:SJ75`, which set `p_inf` and `p_min` as in its Table 2): an `L × L`
grove whose trees are diseased with probability `p_inf` (and then stressed)
and otherwise stressed with probability `1 - p_nst(0)`, drawn from `seed`,
for 20 months. The grove is the corner `1:L, 1:L` of a `2L × 2L` torus whose
other cells are masked healthy trees, so flights leaving the grove are lost.
`kernel` defaults to [`cvc_flight_kernel`](@ref) with radius `L - 1`.
"""
function problem_cvc_citrus(; grove::Symbol=:SJ01, L::Integer=200, kernel=nothing, tspan=20,
                            seed::Integer=0, pmap=(;), kwargs...)
    haskey(CVC_GROVES, grove) ||
        throw(ArgumentError("grove must be one of $(keys(CVC_GROVES)), got :$grove"))
    g = CVC_GROVES[grove]
    kernel = kernel === nothing ? cvc_flight_kernel(; radius=L - 1) : kernel
    size(kernel, 1) <= 2L - 1 ||
        throw(ArgumentError("a kernel wider than 2L - 1 would wrap around the torus"))
    model = build_cvc_citrus(; kernel)
    given = (; (Symbol(k) => v for (k, v) in pairs(pmap))...)
    p = merge(defaults(model), (; p_min=g.p_min), given)
    pnst0 = p.p_min + (1 - 2p.p_min) * (1 + sin(p.phase)) / 2
    rng = Random.Xoshiro(seed)
    diseased = fill(false, 2L, 2L)
    stressed = fill(false, 2L, 2L)
    for j in 1:L, i in 1:L
        diseased[i, j] = rand(rng) < g.p_inf
        stressed[i, j] = diseased[i, j] || rand(rng) < 1 - pnst0
    end
    mask = fill(false, 2L, 2L)
    mask[1:L, 1:L] .= true
    return CAProblem(model, (; diseased, stressed), tspan, merge((; p_min=g.p_min), given);
                     space=Torus(2L, 2L), seed, mask, kwargs...)
end

"""
    cvc_incidence(u; L=size(u.diseased, 1) ÷ 2) -> Float64

The share of diseased trees in the grove, the corner `1:L, 1:L` of the state.
"""
function cvc_incidence(u; L::Integer=size(u.diseased, 1) ÷ 2)
    return count(view(u.diseased, 1:L, 1:L)) / L^2
end

"""
    cvc_field_curve(grove, t) -> Float64

The generalised logistic fitted to the observed incidence of a grove by
[Martins2000](@cite) (its Table 1), at month `t`:
`p₁ / (1 + exp(-(p₂ + p₃ t + p₄ t² + p₅ t³)))`. The coefficients are printed
to few digits, so late values are uncertain.
"""
function cvc_field_curve(grove::Symbol, t::Real)
    c = (SJ01=(0.223, -6.000, 1.647, -0.158, 0.005),
         SJ67=(0.908, -6.414, 1.665, -0.153, 0.005),
         SJ71=(1.000, -4.854, 1.307, -0.126, 0.004),
         SJ75=(0.526, -6.638, 1.396, -0.120, 0.004))
    haskey(c, grove) || throw(ArgumentError("grove must be one of $(keys(c)), got :$grove"))
    p1, p2, p3, p4, p5 = c[grove]
    return p1 / (1 + exp(-(p2 + p3 * t + p4 * t^2 + p5 * t^3)))
end

_register!(:cvc_citrus; build=build_cvc_citrus, problem=problem_cvc_citrus)
