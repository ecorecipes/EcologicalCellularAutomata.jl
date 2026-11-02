# Mercader, Siegert, Liebhold & McCullough (2011, Can. J. For. Res.): three ways to slow the local spread
# of the emerald ash borer from an outlier site, simulated with the coupled map lattice of the
# emerald_ash_borer model: removing ash, girdling trap trees that are destroyed with their larvae, and
# protecting trees with a systemic insecticide, in cells within a radius of the origin, from a given year.
# Dispersal may be biased towards cells with girdled or stressed trees. The rule is the base model's
# growth, feeding and development with the treatments' losses; treatments, emergence and dispersal are a
# callback stored with the problem.

const EAB_TREATMENTS = (:none, :ash_reduction, :girdling, :insecticide)

"""
    build_eab_management() -> CAModel

The management model of [Mercader2011b](@cite), the emerald ash borer model
of [Mercader2011](@cite) (see [`build_emerald_ash_borer`](@ref)) with a
treatment applied to cells within a radius of the origin. Besides the base
model's layers, each cell holds whether it is `treated`, the phloem
`protected` by insecticide, the phloem in large trees left to girdle
(`stock`), the phloem `girdled` this year and its phloem at the start less
what was removed or girdled (`ref`), against which feeding alone can count as
stress. Each year, in every cell
at once:

1. the adults lay `N = min(r0 adults, phloem / C)` larvae;
2. if a tree of `girdled` m² was girdled, it holds a share
   `y g / (y g + phloem - g)` of them, where `y = a x^b` is the ratio of the
   density of larvae on girdled to other trees at a density `x = N / phloem`
   (or `a` with `attraction_decline = false`); the tree is destroyed with
   them;
3. larvae laid in `protected` phloem die;
4. the rest eat `C` m² of unprotected phloem each and take one or two years
   to develop, as in the base model.

The treatments, emergence and dispersal are [`eab_manage!`](@ref), a
callback stored with the problem.
"""
function build_eab_management()
    rule = LocalRule(Offsets((1, 0))) do c, nbrs, ctx
        p = ctx.p
        odd = isodd(ctx.t)   # the source's iteration ctx.t + 2
        P = c.phloem
        new = max(min(p.r0 * c.adults, P / p.C), 0.0)
        g = min(c.girdled, P)
        if g > 0 && new > 0
            y = p.attraction_decline ?
                min(p.girdle_a * (new / P)^p.girdle_b, p.girdle_max) :
                p.girdle_max
            new *= (P - g) / (y * g + P - g)
        end
        P -= g
        P > 0 && (new *= 1 - min(c.protected, P) / P)
        phloem = max(P - p.C * new, 0.0)
        p2 = clamp(-0.194 * log((new + 0.001) / (phloem + 0.0001)) + 0.986, 0.0, 1.0)
        return (; phloem, adults=0.0, larvae1=(1 - p2) * new,
                larvae2_odd=odd ? p2 * new : c.larvae2_odd,
                larvae2_even=odd ? c.larvae2_even : p2 * new, border=c.border,
                phloem0=c.phloem0, treated=c.treated, protected=c.protected, stock=c.stock,
                girdled=0.0, ref=c.ref - g)
    end
    base = build_emerald_ash_borer()
    states = [base.states...,
              StateVariable(:treated, false;
                            description="the cell lies in the treated area"),
              StateVariable(:protected, 0.0;
                            description="phloem protected by insecticide (m²)"),
              StateVariable(:stock, 0.0;
                            description="phloem in large trees left to girdle (m²)"),
              StateVariable(:girdled, 0.0; description="phloem girdled this year (m²)"),
              StateVariable(:ref, 0.0;
                            description="phloem against which feeding counts as stress (m²)")]
    parameters = [filter(q -> q.name ∉ (:behaviour, :bias), base.parameters)...,
                  Parameter(:behaviour, :avoidance; description=":standard or :avoidance"),
                  Parameter(:treatment, :none;
                            description=":none, :ash_reduction, :girdling or :insecticide"),
                  Parameter(:radius, 300.0;
                            description="radius of the treated area around the origin (m)"),
                  Parameter(:share, 0.5;
                            description="share of phloem in large trees, which are treated"),
                  Parameter(:start, 0; description="year in which treatment starts"),
                  Parameter(:girdle_area, 10.32;
                            description="phloem of a girdled tree of 27.9 cm DBH (m²)"),
                  Parameter(:girdle_a, 20.468;
                            description="coefficient of the girdled-to-other density ratio"),
                  Parameter(:girdle_b, -0.559;
                            description="exponent of the girdled-to-other density ratio"),
                  Parameter(:girdle_max, 20.468;
                            description="largest girdled-to-other density ratio"),
                  Parameter(:attraction_decline, true;
                            description="girdled trees attract less as larval density rises"),
                  Parameter(:stress_bias, false;
                            description="dispersal biased towards girdled and stressed cells"),
                  Parameter(:stress_factor, 7.0;
                            description="attraction of a girdled or fully stressed cell"),
                  Parameter(:stress_threshold, 0.33;
                            description="share of phloem lost at which a cell is fully stressed"),
                  Parameter(:stress_from, :lost;
                            description=":lost (any loss of phloem is stress) or :eaten (only feeding)")]
    return CAModel(rule; name=:eab_management, states, parameters,
                   description="Mercader et al. (2011) management of emerald ash borer outliers")
end

"""
    eab_attraction(u, p) -> Matrix{Float64}

The attraction of each cell to dispersing adults under the stress bias of
[Mercader2011b](@cite): `stress_factor` for a cell with a girdled tree, and
otherwise rising linearly from 1, when it has lost none of its phloem, to
`stress_factor`, when it has lost `stress_threshold` of it. With
`stress_from = :lost` the loss is of its phloem at the start (`phloem0`),
however it was lost; with `:eaten`, only phloem eaten by larvae counts (the
loss is from `ref`, the phloem at the start less what was removed or
girdled).
"""
function eab_attraction(u, p)
    A = ones(size(u.phloem))
    for I in eachindex(A)
        if u.girdled[I] > 0
            A[I] = p.stress_factor
        else
            base = p.stress_from === :lost ? u.phloem0[I] : u.ref[I]
            base > 0 || continue
            lost = clamp((base - u.phloem[I]) / base, 0.0, 1.0)
            A[I] = 1 + (p.stress_factor - 1) * min(lost / p.stress_threshold, 1.0)
        end
    end
    return A
end

"""
    eab_manage!(u, p, t; larvae=true) -> u

The treatment, emergence and dispersal of [Mercader2011b](@cite) after year
`t`, in place. In the year `start`, every cell within `radius` m of the origin
is marked as treated and

- `:ash_reduction`: loses a `share` of its phloem and of its larvae, once;
- `:girdling`: gets a stock of `share × phloem0` m² of large trees (no more
  than its phloem) to girdle;
- `:insecticide`: has `share × phloem0` m² protected (no more than its
  phloem), and the share of its larvae in that phloem dies.

From then on, under `:girdling`, a treated cell girdles one tree of
`girdle_area` m² a year from its stock while it lasts. The adults then
emerge and disperse as in [`eab_disperse!`](@ref), weighted by
[`eab_attraction`](@ref) if `stress_bias` is set. With `larvae = false` no
larvae are removed (for treatments before the founders are released).
"""
function eab_manage!(u, p, t::Integer; larvae::Bool=true)
    if t == p.start && p.treatment !== :none
        N = size(u.phloem, 1)
        c = (N + 1) ÷ 2
        for I in CartesianIndices(u.phloem)
            u.border[I] && continue
            hypot(I[1] - c, I[2] - c) * p.cell <= p.radius || continue
            u.treated[I] = true
            if p.treatment === :ash_reduction
                removed = p.share * u.phloem[I]
                u.phloem[I] -= removed
                u.ref[I] -= removed
                if larvae
                    for L in (u.larvae1, u.larvae2_odd, u.larvae2_even)
                        L[I] *= 1 - p.share
                    end
                end
            elseif p.treatment === :girdling
                u.stock[I] = min(p.share * u.phloem0[I], u.phloem[I])
            elseif p.treatment === :insecticide
                u.protected[I] = min(p.share * u.phloem0[I], u.phloem[I])
                if larvae && u.phloem[I] > 0
                    q = u.protected[I] / u.phloem[I]
                    for L in (u.larvae1, u.larvae2_odd, u.larvae2_even)
                        L[I] *= 1 - q
                    end
                end
            end
        end
    end
    if p.treatment === :girdling && t >= p.start
        for I in eachindex(u.girdled)
            u.treated[I] || continue
            g = min(p.girdle_area, u.stock[I], u.phloem[I])
            u.girdled[I] = g
            u.stock[I] -= g
        end
    end
    attraction = p.stress_bias ? eab_attraction(u, p) : nothing
    return eab_disperse!(u, merge(p, (; bias=0.0)), t; attraction)
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
                           stress_bias=false, n=101, founders=100, years=15,
                           pmap=(;), kwargs...) -> CAProblem

The setup of [Mercader2011b](@cite): a homogeneous landscape of `n × n` cells
of 50 m with 200 m² of ash phloem each, surrounded by the base model's border
of 2 cells, and `founders` adults at its centre that emerge and disperse in
the first year, for `years` years, with the `treatment` applied within
`radius` m of the centre from year `start` (0 is before the founders
disperse). The model is deterministic.
"""
function problem_eab_management(; treatment::Symbol=:none, radius::Real=300,
                                share::Real=0.5, start::Integer=0, stress_bias::Bool=false,
                                n::Integer=101, founders::Real=100, years::Integer=15,
                                seed::Integer=0, pmap=(;), kwargs...)
    treatment in EAB_TREATMENTS ||
        throw(ArgumentError("treatment must be one of $(EAB_TREATMENTS), got :$treatment"))
    0 <= share <= 1 || throw(ArgumentError("share must lie in [0, 1], got $share"))
    start >= 0 || throw(ArgumentError("start must not be negative, got $start"))
    isodd(n) || throw(ArgumentError("n must be odd, so the landscape has a centre cell"))
    N = n + 4
    phloem = fill(200.0, N, N)
    border = fill(true, N, N)
    border[3:(n + 2), 3:(n + 2)] .= false
    larvae1 = zeros(N, N)
    larvae1[(N + 1) ÷ 2, (N + 1) ÷ 2] = founders
    model = build_eab_management()
    given = merge((; treatment, radius=Float64(radius), share=Float64(share), start,
                   stress_bias), (; (Symbol(k) => v for (k, v) in pairs(pmap))...))
    p = merge(defaults(model), given)
    p.behaviour in (:standard, :avoidance) ||
        throw(ArgumentError("behaviour must be :standard or :avoidance, got :$(p.behaviour)"))
    u0 = (; phloem, adults=zeros(N, N), larvae1, larvae2_odd=zeros(N, N),
          larvae2_even=zeros(N, N), border, phloem0=copy(phloem), treated=fill(false, N, N),
          protected=zeros(N, N), stock=zeros(N, N), girdled=zeros(N, N), ref=copy(phloem))
    eab_manage!(u0, p, 0; larvae=false)
    return CAProblem(model, u0, years, given; space=Grid(N, N), seed,
                     callback=eab_management_callback(), kwargs...)
end

"""
    eab_beyond(u, radius; cell=50) -> Float64

The number of larvae (next years' beetles) in cells more than `radius` m from
the origin, the source's measure of spread beyond the treated area.
"""
function eab_beyond(u, radius::Real; cell::Real=50)
    N = size(u.phloem, 1)
    c = (N + 1) ÷ 2
    total = 0.0
    for I in CartesianIndices(u.phloem)
        (u.border[I] || hypot(I[1] - c, I[2] - c) * cell <= radius) && continue
        total += u.larvae1[I] + u.larvae2_odd[I] + u.larvae2_even[I]
    end
    return total
end

_register!(:eab_management; build=build_eab_management, problem=problem_eab_management)
