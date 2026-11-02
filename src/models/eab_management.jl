# Mercader, Siegert, Liebhold & McCullough (2011, Can. J. For. Res.): three ways to slow the local spread
# of the emerald ash borer from an outlier site, simulated with the coupled map lattice of the
# emerald_ash_borer model: removing ash, girdling trap trees that are destroyed with their larvae, and
# protecting trees with a systemic insecticide, in cells within a radius of the origin, from a given year.
# Dispersal may be biased towards cells with girdled or stressed trees. The model follows the source's
# code (its Supplemental Appendix S1): the rule is the year's egg laying, treatments, feeding and
# development; the stress bias, emergence and dispersal are a callback stored with the problem.

const EAB_TREATMENTS = (:none, :ash_reduction, :girdling, :insecticide)

"""
    build_eab_management() -> CAModel

The management model of [Mercader2011b](@cite), the emerald ash borer model
of [Mercader2011](@cite) (see [`build_emerald_ash_borer`](@ref)) with a
treatment applied to the cells within a radius of the origin, as the source's
code computes it. Besides the base model's layers, each cell holds whether it
is `treated` and its phloem in `large` trees (m²), which alone can be
treated. Each year `t`, in every cell at once, the adults lay
`N = min(r0 adults, phloem / C)` larvae; from the year `start` on,

1. under girdling, a treated cell girdles a tree of `girdle_area` m² (no more
   than its large-tree phloem), which takes a share
   `y g / (y g + phloem - g)` of the larvae, where
   `y = max(girdle_a (x + girdle_shift)^girdle_b, girdle_min)` is the ratio of
   the density of larvae on girdled to other trees at the density
   `x = N / phloem` (or `girdle_a` with `attraction_decline = false`); the tree
   is removed with them;
2. under insecticide, the larvae laid in a treated cell's large-tree phloem
   die (`x` times its area);
3. the larvae left eat `C` m² of phloem each, the insecticide-protected share
   excepted, as in the source's code, which reduces the feeding by the
   unprotected share a second time;
4. under ash reduction, in the year `start` only, a treated cell's remaining
   large-tree phloem is removed after the larvae have fed;

and the larvae take one or two years to develop, as in the base model. The
large-tree phloem loses its share of what is eaten. The stress bias,
emergence and dispersal are [`eab_manage!`](@ref), a callback stored with the
problem.
"""
function build_eab_management()
    rule = LocalRule(Offsets((1, 0))) do c, nbrs, ctx
        p = ctx.p
        t = ctx.t
        odd = isodd(t)   # the source's iteration t + 2
        P = c.phloem
        Lg = c.large
        new = max(min(p.r0 * c.adults, P / p.C), 0.0)
        if t >= p.start
            trap = p.treatment === :girdling && c.treated ? p.girdle_area : 0.0
            protect = p.treatment === :insecticide && c.treated ? Lg : 0.0
            x = new / (P + 0.0001)
            y = p.attraction_decline ?
                clamp(p.girdle_a * (x + p.girdle_shift)^p.girdle_b, p.girdle_min,
                      p.girdle_max) :
                p.girdle_a
            g = min(Lg, trap)
            divisor = g * y + (P - g)
            divisor <= 0 && (divisor = 0.0001)
            new = max(new - new * (g * y / divisor) - x * protect, 0.0)
            unprotected = max((P - protect + 0.0001) / (P + 0.0001), 0.0)
            phloem = max(P - unprotected * p.C * new - g, 0.0)
            large = protect > 0 ? max(Lg - g - Lg / (P + 0.0001) * p.C * new, protect) :
                    Lg - g - (Lg - g + 0.0001) / (P + 0.01) * p.C * new
            large = max(large, 0.0)
        else
            phloem = max(P - p.C * new, 0.0)
            large = max(Lg - (Lg + 0.01) / (P + 0.01) * p.C * new, 0.0)
        end
        if p.treatment === :ash_reduction && t == p.start
            cut = c.treated ? large : 0.0
            phloem = max(phloem - cut, 0.0)
            large = large - cut < 0 ? 0.0 : phloem   # as the source's code resets it
        end
        p2 = clamp(-0.194 * log((new + 0.001) / (phloem + 0.0001)) + 0.986, 0.0, 1.0)
        return (; phloem, adults=0.0, larvae1=(1 - p2) * new,
                larvae2_odd=odd ? p2 * new : c.larvae2_odd,
                larvae2_even=odd ? c.larvae2_even : p2 * new, border=c.border,
                phloem0=c.phloem0, treated=c.treated, large)
    end
    base = build_emerald_ash_borer()
    states = [base.states...,
              StateVariable(:treated, false;
                            description="the cell lies in the treated area"),
              StateVariable(:large, 0.0; description="phloem in large trees (m²)")]
    parameters = [filter(q -> q.name ∉ (:behaviour, :bias, :window, :min_source),
                         base.parameters)...,
                  Parameter(:behaviour, :avoidance; description=":standard or :avoidance"),
                  Parameter(:window, 40;
                            description="half-width of the dispersal window (cells)"),
                  Parameter(:min_source, 0.9;
                            description="adults a cell must hold to disperse them; smaller numbers are lost"),
                  Parameter(:treatment, :none;
                            description=":none, :ash_reduction, :girdling or :insecticide"),
                  Parameter(:radius, 300.0;
                            description="radius of the treated area around the origin (m)"),
                  Parameter(:share, 0.5;
                            description="share of phloem in large trees, which are treated"),
                  Parameter(:large_min, 6.32;
                            description="large-tree phloem (m²) below which a cell has none"),
                  Parameter(:start, 0; description="year in which treatment starts"),
                  Parameter(:girdle_area, 10.32;
                            description="phloem of a girdled tree of 27.9 cm DBH (m²)"),
                  Parameter(:girdle_a, 18.25;
                            description="coefficient of the girdled-to-other density ratio"),
                  Parameter(:girdle_shift, 1.0;
                            description="added to the density in the girdled-to-other density ratio"),
                  Parameter(:girdle_b, -0.533;
                            description="exponent of the girdled-to-other density ratio"),
                  Parameter(:girdle_min, 1.0;
                            description="smallest girdled-to-other density ratio"),
                  Parameter(:girdle_max, Inf;
                            description="largest girdled-to-other density ratio"),
                  Parameter(:attraction_decline, true;
                            description="girdled trees attract less as larval density rises"),
                  Parameter(:stress_bias, false;
                            description="dispersal biased towards girdled and stressed cells"),
                  Parameter(:stress_factor, 7.0;
                            description="attraction of a girdled or fully stressed cell"),
                  Parameter(:stress_threshold, 1 / 3;
                            description="share of phloem lost at which a cell is fully stressed")]
    return CAModel(rule; name=:eab_management, states, parameters,
                   description="Mercader et al. (2011) management of emerald ash borer outliers")
end

"""
    eab_attraction(u, p) -> Matrix{Float64}

The attraction of each cell to dispersing adults under the stress bias of
[Mercader2011b](@cite), as the source's code computes it: under girdling,
`stress_factor` for every treated cell, from the first year; otherwise
`stress_factor × lost / (stress_threshold × max(phloem0, 3))`, where `lost`
is the phloem the cell has lost since the start (`phloem0`) by any means,
kept within 1 and `stress_factor`.
"""
function eab_attraction(u, p)
    A = ones(size(u.phloem))
    for I in eachindex(A)
        if p.treatment === :girdling && u.treated[I]
            A[I] = max(p.stress_factor, 1.0)
        else
            lost = u.phloem0[I] - u.phloem[I]
            s = p.stress_factor * lost / (p.stress_threshold * max(u.phloem0[I], 3.0))
            A[I] = max(min(s, p.stress_factor), 1.0)
        end
    end
    return A
end

"""
    eab_manage!(u, p, t) -> u

The emergence and dispersal of [Mercader2011b](@cite) after year `t`, in
place: the adults emerge and disperse as in [`eab_disperse!`](@ref), weighted
by [`eab_attraction`](@ref) if `stress_bias` is set; then, under ash
reduction in the year `start`, the adults landing in a treated cell lose the
share on its large trees, which are removed.
"""
function eab_manage!(u, p, t::Integer)
    attraction = p.stress_bias ? eab_attraction(u, p) : nothing
    eab_disperse!(u, merge(p, (; bias=0.0)), t; attraction)
    if p.treatment === :ash_reduction && t == p.start
        for I in eachindex(u.adults)
            u.treated[I] || continue
            u.adults[I] -= u.large[I] * u.adults[I] / (u.phloem[I] + 0.0001)
        end
    end
    return u
end

"""
    eab_management_callback() -> DiscreteCallback

A callback that runs [`eab_manage!`](@ref) after each year.
"""
function eab_management_callback()
    affect! = function (integ)
        eab_manage!(integ.u, integ.p, Int(integ.t))
        return nothing
    end
    return DiscreteCallback((u, t, integ) -> true, affect!; save_positions=(false, false))
end

"""
    problem_eab_management(; treatment=:none, radius=300, share=0.5, start=0,
                           stress_bias=false, n=100, phloem=200, founders=100, years=15,
                           pmap=(;), kwargs...) -> CAProblem

The setup of [Mercader2011b](@cite): a landscape of `n × n` cells of 50 m
with `phloem` m² of ash phloem each (or a matrix of it), a `share` of it in
large trees (none in a cell where that is at most `large_min` m²), and
`founders` adults at the origin, the cell `(n ÷ 2 + 1, n ÷ 2 + 1)`, which
emerge and disperse in the first year, for `years` years, with the
`treatment` applied to the cells within `radius` m of the origin from year
`start` (0 is the founders' first year). As in the source's code, there is
no border: dispersal is normalised within the landscape. The model is
deterministic.
"""
function problem_eab_management(; treatment::Symbol=:none, radius::Real=300,
                                share::Real=0.5, start::Integer=0, stress_bias::Bool=false,
                                n::Integer=100, phloem=200.0, founders::Real=100,
                                years::Integer=15, seed::Integer=0, pmap=(;), kwargs...)
    treatment in EAB_TREATMENTS ||
        throw(ArgumentError("treatment must be one of $(EAB_TREATMENTS), got :$treatment"))
    0 <= share <= 1 || throw(ArgumentError("share must lie in [0, 1], got $share"))
    start >= 0 || throw(ArgumentError("start must not be negative, got $start"))
    n >= 1 || throw(ArgumentError("n must be positive, got $n"))
    P = phloem isa Real ? fill(Float64(phloem), n, n) : Float64.(phloem)
    size(P) == (n, n) || throw(ArgumentError("phloem must be $n × $n, got $(size(P))"))
    model = build_eab_management()
    given = merge((; treatment, radius=Float64(radius), share=Float64(share), start,
                   stress_bias), (; (Symbol(k) => v for (k, v) in pairs(pmap))...))
    p = merge(defaults(model), given)
    p.behaviour in (:standard, :avoidance) ||
        throw(ArgumentError("behaviour must be :standard or :avoidance, got :$(p.behaviour)"))
    c = n ÷ 2 + 1
    larvae1 = zeros(n, n)
    larvae1[c, c] = founders
    large = [x > p.large_min ? x : 0.0 for x in p.share .* P]
    treated = [hypot(I[1] - c, I[2] - c) * p.cell <= p.radius for I in CartesianIndices(P)]
    u0 = (; phloem=copy(P), adults=zeros(n, n), larvae1, larvae2_odd=zeros(n, n),
          larvae2_even=zeros(n, n), border=fill(false, n, n), phloem0=copy(P), treated,
          large)
    eab_manage!(u0, p, 0)
    return CAProblem(model, u0, years, given; space=Grid(n, n), seed,
                     callback=eab_management_callback(), kwargs...)
end

"""
    eab_beyond(u, radius; cell=50) -> Float64

The number of larvae (next years' beetles) in cells more than `radius` m from
the origin, the source's measure of spread beyond the treated area.
"""
function eab_beyond(u, radius::Real; cell::Real=50)
    N = size(u.phloem, 1)
    c = N ÷ 2 + 1
    total = 0.0
    for I in CartesianIndices(u.phloem)
        (u.border[I] || hypot(I[1] - c, I[2] - c) * cell <= radius) && continue
        total += u.larvae1[I] + u.larvae2_odd[I] + u.larvae2_even[I]
    end
    return total
end

_register!(:eab_management; build=build_eab_management, problem=problem_eab_management)
