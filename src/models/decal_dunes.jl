# DECAL, the Discrete Ecogeomorphic Aeolian Landscapes model of Baas & Nield (2007) and Nield & Baas
# (2008): Werner's sand-slab automaton with vegetation. Slabs of sand stand in stacks on a grid; each
# iteration polls every cell once in random order, and a polled cell outside the lee of upwind
# topography loses its top slab with probability pe - ρ; the slab moves one cell downwind at a time
# until it is deposited, with probability pd + ρ (1 - pd) (pd 0.6 on sand, 0.4 on a bare cell),
# certainly in a shadow zone. Avalanches keep slopes within the angle of repose, 30 degrees, or 40
# where ρ >= 0.3. Each cell holds an effectiveness ρ for each plant type, and their sum, within
# [0, 1], acts on transport; once a year each type's ρ changes by its growth function of the cell's
# net sedimentation balance. Transport is sequential, so each iteration is a callback.

"""
    DecalGrowth(vertices; lower, upper)

A DECAL growth function: the yearly change in a plant type's vegetation
effectiveness ρ as a function of the cell's net sedimentation balance over the
year (m), piecewise linear through five `vertices` `(balance, change)` (A to E
of [Nield2008](@cite), Figure 5), constant beyond A and E, with the change
normalised by the type's physiological range `[lower, upper]`, to which ρ is
held. [`decal_growth`](@ref) evaluates it.
"""
struct DecalGrowth
    vertices::NTuple{5,NTuple{2,Float64}}
    lower::Float64
    upper::Float64
    function DecalGrowth(vertices; lower::Real, upper::Real)
        v = Tuple(Tuple(Float64.(x)) for x in vertices)
        length(v) == 5 || throw(ArgumentError("a growth function needs 5 vertices"))
        issorted(first.(v)) ||
            throw(ArgumentError("the vertices' sedimentation balances must increase"))
        lower < upper || throw(ArgumentError("lower must be below upper"))
        return new(v, Float64(lower), Float64(upper))
    end
end

"""
    DecalGrowth(g::DecalGrowth; peak) -> DecalGrowth

`g` with the normalised peak growth (vertex C) set to `peak`, as in the source's
phase diagram over peak growth rates.
"""
function DecalGrowth(g::DecalGrowth; peak::Real)
    v = collect(g.vertices)
    v[3] = (v[3][1], Float64(peak))
    return DecalGrowth(v; lower=g.lower, upper=g.upper)
end

"""
    decal_growth(g::DecalGrowth, balance) -> Float64

The yearly change in vegetation effectiveness for a net sedimentation balance
`balance` (m, positive for deposition): the normalised change interpolated
between the vertices, times the width `upper - lower` of the physiological
range.
"""
function decal_growth(g::DecalGrowth, b::Real)
    v = g.vertices
    y = if b <= v[1][1]
        v[1][2]
    elseif b >= v[5][1]
        v[5][2]
    else
        k = findlast(x -> x[1] <= b, v)
        (x0, y0), (x1, y1) = v[k], v[k + 1]
        x1 == x0 ? y1 : y0 + (y1 - y0) * (b - x0) / (x1 - x0)
    end
    return (g.upper - g.lower) * y
end

"""
    DECAL_GRASS, DECAL_SHRUB, DECAL_MESQUITE

The growth functions of [Nield2008](@cite), Figure 5, read from the figure: a
pioneer grass that needs fresh sand (range [0, 1], declining without
deposition, peak 0.48 at 0.6 m of deposition, dying under more than 2 m), a stabilizing shrub that prefers a
quiet surface (range [-0.5, 1.5], peak 0.04 with no change), and a
burial-tolerant shrub such as mesquite (range [-1, 3], no change without
sedimentation, peak 0.25 at 0.5 m). The figure's axes were calibrated in
pixels and each vertex read to about 0.02.
"""
const DECAL_GRASS = DecalGrowth(((-1.39, -1.0), (0.09, 0.0), (0.6, 0.48), (2.0, 0.0),
                                 (2.21, -1.0)); lower=0, upper=1)
const DECAL_SHRUB = DecalGrowth(((-1.29, -0.5), (-0.64, 0.0), (0.0, 0.04), (0.21, 0.0),
                                 (2.82, -0.5)); lower=-0.5, upper=1.5)
const DECAL_MESQUITE = DecalGrowth(((-0.5, -0.25), (0.0, 0.0), (0.5, 0.25), (2.76, 0.0),
                                    (3.5, -0.08)); lower=-1, upper=3)
@doc (@doc DECAL_GRASS) DECAL_SHRUB
@doc (@doc DECAL_GRASS) DECAL_MESQUITE

const _DECAL_TYPES = (:grass, :shrub, :mesquite)

"""
    build_decal_dunes() -> CAModel

DECAL [Baas2007](@cite) [Nield2008](@cite) on a grid whose second index runs
downwind. Layers: `sand` (slabs above the erosion base), the vegetation
effectiveness of each plant type (`grass`, `shrub`, `mesquite`) and `balance`
(slabs gained this year). The model's local rule leaves cells unchanged; the
dynamics are [`decal_callback`](@ref), stored with the problem, which runs
one iteration ([`decal_step!`](@ref)) per step.
"""
function build_decal_dunes()
    rule = LocalRule((c, nbrs, ctx) -> c, VonNeumann(1))
    states = [StateVariable(:sand, Int32(0); description="slabs above the erosion base"),
              StateVariable(:grass, 0.0; description="effectiveness of the pioneer grass"),
              StateVariable(:shrub, 0.0;
                            description="effectiveness of the stabilizing shrub"),
              StateVariable(:mesquite, 0.0;
                            description="effectiveness of the burial-tolerant shrub"),
              StateVariable(:balance, Int32(0);
                            description="slabs gained (lost if negative) this year")]
    parameters = [Parameter(:cell, 1.0; description="cell width (m)"),
                  Parameter(:slab, 0.1; description="slab height (m)"),
                  Parameter(:pe, 1.0; description="erosion probability on bare sand"),
                  Parameter(:pd_sand, 0.6; description="deposition probability on sand"),
                  Parameter(:pd_bare, 0.4;
                            description="deposition probability on a cell with no sand"),
                  Parameter(:shadow, 15.0; description="shadow angle (degrees)"),
                  Parameter(:repose_bare, 30.0;
                            description="angle of repose, bare (degrees)"),
                  Parameter(:repose_veg, 40.0;
                            description="angle of repose where ρ >= rho_repose (degrees)"),
                  Parameter(:rho_repose, 0.3;
                            description="effectiveness that steepens the angle of repose"),
                  Parameter(:iterations, 80; description="iterations per year"),
                  Parameter(:types, (:grass, :shrub);
                            description="the plant types that grow"),
                  Parameter(:grass_growth, DECAL_GRASS;
                            description="growth function of the grass"),
                  Parameter(:shrub_growth, DECAL_SHRUB;
                            description="growth function of the shrub"),
                  Parameter(:mesquite_growth, DECAL_MESQUITE;
                            description="growth function of the mesquite"),
                  Parameter(:open, false;
                            description="slabs leave at the downwind edge instead of wrapping")]
    return CAModel(rule; name=:decal_dunes, states, parameters,
                   description="DECAL (Baas & Nield 2007; Nield & Baas 2008)")
end

_decal_rho(u, r, x) = clamp(u.grass[r, x] + u.shrub[r, x] + u.mesquite[r, x], 0.0, 1.0)

"""
    decal_shadow(sand, r, x, p, top=maximum(sand)) -> Bool

Whether cell `(r, x)` lies in a shadow zone: some cell `k` cells upwind in its
row (wrapping around a periodic grid) stands more than `k tan(shadow)` above it,
in slabs of `p.slab` m on cells `p.cell` m wide. `top`, an upper bound on the
heights, limits how far upwind to look.
"""
function decal_shadow(sand, r, x, p, top=maximum(sand))
    n = size(sand, 2)
    s = tand(p.shadow) * p.cell / p.slab
    h = sand[r, x]
    kmax = min(n - 1, floor(Int, (top - h) / s))
    for k in 1:kmax
        y = x - k
        if y < 1
            p.open && break
            y += n
        end
        sand[r, y] - h > k * s && return true
    end
    return false
end

@inline _decal_wrap(i, n) = mod1(i, n)

function _decal_limit(u, r, x, p)
    θ = _decal_rho(u, r, x) >= p.rho_repose ? p.repose_veg : p.repose_bare
    return tand(θ) * p.cell / p.slab
end

"""
    decal_avalanche!(u, p, rng, start) -> Int

Avalanching from the cells in `start`: a cell whose drop to its steepest lower
von Neumann neighbour (wrapping around the grid, or not across the downwind
edge of an open grid) exceeds the angle of repose (`tan(θ) cell / slab`
slabs) passes that neighbour one slab, ties broken at random, and the two
cells and the first's neighbours are checked again, until no slope exceeds it.
Returns the number of slabs moved.
"""
function decal_avalanche!(u, p, rng, start)
    sand = u.sand
    m, n = size(sand)
    work = collect(start)
    moved = 0
    while !isempty(work)
        (r, x) = pop!(work)
        sand[r, x] > 0 || continue
        lim = _decal_limit(u, r, x, p)
        best = 0
        ties = 0
        tr, tx = r, x
        for (dr, dx) in ((-1, 0), (1, 0), (0, -1), (0, 1))
            xx = x + dx
            (p.open && !(1 <= xx <= n)) && continue
            rr, xx = _decal_wrap(r + dr, m), _decal_wrap(xx, n)
            drop = Int(sand[r, x]) - Int(sand[rr, xx])
            drop > lim || continue
            if drop > best
                best, ties, tr, tx = drop, 1, rr, xx
            elseif drop == best
                ties += 1
                rand(rng) < 1 / ties && ((tr, tx) = (rr, xx))
            end
        end
        ties == 0 && continue
        sand[r, x] -= 1
        sand[tr, tx] += 1
        u.balance[r, x] -= 1
        u.balance[tr, tx] += 1
        moved += 1
        push!(work, (tr, tx), (r, x))
        for (dr, dx) in ((-1, 0), (1, 0), (0, -1), (0, 1))
            push!(work, (_decal_wrap(r + dr, m), _decal_wrap(x + dx, n)))
        end
    end
    return moved
end

"""
    decal_poll!(u, p, rng, r, x, top=Ref(maximum(u.sand))) -> Int

One polling event at cell `(r, x)`: if it has sand and lies outside a shadow
zone, its top slab is eroded with probability `pe - ρ`; the slab moves one
cell downwind at a time and is deposited at a cell in a shadow zone, or with
probability `pd + ρ (1 - pd)` (`pd_sand` on a cell with sand, `pd_bare`
otherwise), and avalanching follows the erosion and the deposition. Returns
the number of cells the slab travelled (0 if none was eroded, `-1` if it left
an open grid). `top` holds an upper bound on the heights, raised by
depositions (avalanches never raise a cell above the one it came from).
"""
function decal_poll!(u, p, rng, r, x, top=Ref(maximum(u.sand)))
    sand = u.sand
    m, n = size(sand)
    sand[r, x] > 0 || return 0
    decal_shadow(sand, r, x, p, top[]) && return 0
    rand(rng) < p.pe - _decal_rho(u, r, x) || return 0
    sand[r, x] -= 1
    u.balance[r, x] -= 1
    decal_avalanche!(u, p, rng,
                     ((r, x), (_decal_wrap(r - 1, m), x), (_decal_wrap(r + 1, m), x),
                      (r, _decal_wrap(x - 1, n)), (r, _decal_wrap(x + 1, n))))
    y = x
    k = 0
    while true
        k += 1
        y += 1
        if y > n
            p.open && return -1
            y = 1
        end
        decal_shadow(sand, r, y, p, top[]) && break
        pd = sand[r, y] > 0 ? p.pd_sand : p.pd_bare
        ρ = _decal_rho(u, r, y)
        rand(rng) < pd + ρ * (1 - pd) && break
    end
    sand[r, y] += 1
    u.balance[r, y] += 1
    top[] = max(top[], sand[r, y])
    decal_avalanche!(u, p, rng, ((r, y),))
    return k
end

"""
    decal_step!(u, p, t, rng) -> NamedTuple

Iteration `t`: every cell polled once, in an order drawn from `rng`
([`decal_poll!`](@ref)); after the last iteration of a year (`t` a multiple of
`p.iterations`), each growing type's effectiveness changes by its growth
function of the year's balance (in m), within its physiological range, the
balance is reset and the whole grid avalanches. Returns the numbers of slabs
eroded, avalanched in the year's end, and lost at an open edge.
"""
function decal_step!(u, p, t::Integer, rng)
    m, n = size(u.sand)
    eroded = 0
    lost = 0
    top = Ref(maximum(u.sand))
    for i in Random.randperm(rng, m * n)
        r, x = Tuple(CartesianIndices((m, n))[i])
        k = decal_poll!(u, p, rng, r, x, top)
        k != 0 && (eroded += 1)
        k == -1 && (lost += 1)
    end
    avalanched = 0
    if t % p.iterations == 0
        for type in _DECAL_TYPES
            type in p.types || continue
            g = getproperty(p, Symbol(type, :_growth))
            ρ = getproperty(u, type)
            for I in CartesianIndices(ρ)
                ρ[I] = clamp(ρ[I] + decal_growth(g, u.balance[I] * p.slab), g.lower,
                             g.upper)
            end
        end
        fill!(u.balance, 0)
        avalanched = decal_avalanche!(u, p, rng, Tuple.(vec(CartesianIndices((m, n)))))
        fill!(u.balance, 0)
    end
    return (; eroded, avalanched, lost)
end

"""
    decal_callback() -> DiscreteCallback

A callback that runs [`decal_step!`](@ref) after each step, with a random
number generator seeded from the problem's seed and the iteration.
"""
function decal_callback()
    affect! = function (integ)
        t = Int(integ.t)
        decal_step!(integ.u, integ.p, t, Random.Xoshiro(hash((integ.prob.seed, t, :decal))))
        return nothing
    end
    return DiscreteCallback((u, t, integ) -> true, affect!; save_positions=(false, false))
end

"""
    problem_decal_dunes(; scenario=:nebkha, cell=nothing, size=100.0, years=40, flux=5.0,
                        seed=0, pmap=(;), kwargs...) -> CAProblem

The setups of [Nield2008](@cite) and [Baas2007](@cite):

- `:nebkha`: a flat bare bed 0.3 m deep on a periodic `size` × `size` m grid
  (cells 1.25 m by default), mesquite only, slabs 0.15 m, and as many
  iterations a year as give a flux of `flux` m³ per m per year,
  `round(flux pd_sand / (slab cell))` (16 for 1.25 m cells): the same physical
  conditions at any cell size;
- `:parabolic`: a flat bed `depth` m deep (2), 1 m cells, slabs 0.1 m, 80 iterations a
  year, with shrub and grass effectiveness drawn uniformly from `[0, shrub_max]`
  and `[0, grass_max]` (0.9 and 0.1, as in the source's Figure 11a) except in
  five bare blowouts 6 m across, on a grid 100 m across and 300 m downwind whose
  downwind edge is open;
- `:bare`: Werner's bare-sand model on a periodic grid with seven mounds,
  slabs 0.1 cell widths, `pd_sand` 0.6.

`pmap` overrides parameters; `years` sets the time span.
"""
function problem_decal_dunes(; scenario::Symbol=:nebkha, cell=nothing, size::Real=100.0,
                             years::Integer=40, flux::Real=5.0, shrub_max::Real=0.9,
                             grass_max::Real=0.1, depth::Real=2.0, seed::Integer=0,
                             pmap=(;),
                             kwargs...)
    scenario in (:nebkha, :parabolic, :bare) ||
        throw(ArgumentError("scenario must be :nebkha, :parabolic or :bare, got :$scenario"))
    rng = Random.Xoshiro(hash((seed, 0, :decal)))
    given = (; (Symbol(k) => v for (k, v) in pairs(pmap))...)
    if scenario === :nebkha
        c = cell === nothing ? 1.25 : Float64(cell)
        slab = get(given, :slab, 0.15)
        pd = get(given, :pd_sand, 0.6)
        n = round(Int, size / c)
        iterations = max(1, round(Int, flux * pd / (slab * c)))
        p = merge((; cell=c, slab, iterations, types=(:mesquite,)), given)
        sand = fill(Int32(round(0.3 / slab)), n, n)
        grass = zeros(n, n)
        shrub = zeros(n, n)
        u0 = (; sand, grass, shrub, mesquite=zeros(n, n))
    elseif scenario === :parabolic
        c = cell === nothing ? 1.0 : Float64(cell)
        slab = get(given, :slab, 0.1)
        m, n = round(Int, 100 / c), round(Int, 300 / c)
        p = merge((; cell=c, slab, iterations=80, types=(:grass, :shrub), open=true), given)
        sand = fill(Int32(round(depth / slab)), m, n)
        shrub = [shrub_max * rand(rng) for _ in 1:m, _ in 1:n]
        grass = [grass_max * rand(rng) for _ in 1:m, _ in 1:n]
        w = round(Int, 6 / c)
        for k in 1:5
            r0 = round(Int, (k - 0.5) * m / 5 - w / 2) + 1
            x0 = round(Int, 20 / c)
            shrub[r0:(r0 + w - 1), x0:(x0 + w - 1)] .= 0
            grass[r0:(r0 + w - 1), x0:(x0 + w - 1)] .= 0
        end
        u0 = (; sand, grass, shrub, mesquite=zeros(m, n))
    else
        c = cell === nothing ? 1.0 : Float64(cell)
        n = round(Int, size / c)
        p = merge((; cell=c, slab=0.1 * c, iterations=typemax(Int32), types=()), given)
        sand = zeros(Int32, n, n)
        for _ in 1:7
            r0, x0 = rand(rng, 1:n), rand(rng, 1:n)
            for dr in -6:6, dx in -6:6
                d = hypot(dr, dx)
                d <= 6 && (sand[mod1(r0 + dr, n), mod1(x0 + dx, n)] += Int32(round(6 - d)))
            end
        end
        u0 = (; sand, grass=zeros(n, n), shrub=zeros(n, n), mesquite=zeros(n, n))
    end
    model = build_decal_dunes()
    return CAProblem(model, u0, years * p.iterations, p; space=Grid(Base.size(u0.sand)...),
                     seed, callback=decal_callback(), kwargs...)
end

"""
    decal_mounds(u, p; above=0.0) -> NamedTuple

The mounds of a state: the groups of cells (von Neumann neighbours, wrapping
around the grid) whose sand stands more than `above` m above the erosion base;
with the default, the patches of sand on a bare base, which are the nebkhas of
a nebkha field. Returns their `count`, the `median_length` and `median_width`
(m; each group's extent along and across the wind), their `median_height`
(m; the highest cell of each), the `lengths` of all of them, and the
`relief` of the surface (m).
"""
function decal_mounds(u, p; above::Real=0.0)
    sand = u.sand
    m, n = size(sand)
    mask = sand .* p.slab .> above
    seen = falses(m, n)
    lengths, widths, heights = Float64[], Float64[], Float64[]
    for I in CartesianIndices(mask)
        (mask[I] && !seen[I]) || continue
        stack = [I]
        seen[I] = true
        rows, cols = Set{Int}(), Set{Int}()
        top = 0
        while !isempty(stack)
            J = pop!(stack)
            push!(rows, J[1])
            push!(cols, J[2])
            top = max(top, Int(sand[J]))
            for (dr, dx) in ((-1, 0), (1, 0), (0, -1), (0, 1))
                K = CartesianIndex(mod1(J[1] + dr, m), mod1(J[2] + dx, n))
                if mask[K] && !seen[K]
                    seen[K] = true
                    push!(stack, K)
                end
            end
        end
        push!(lengths, length(cols) * p.cell)
        push!(widths, length(rows) * p.cell)
        push!(heights, top * p.slab)
    end
    med(v) = isempty(v) ? 0.0 : sort(v)[(length(v) + 1) ÷ 2]
    z = sand .* p.slab
    return (count=length(lengths), median_length=med(lengths), median_width=med(widths),
            median_height=med(heights), lengths=lengths, relief=maximum(z) - minimum(z))
end

_register!(:decal_dunes; build=build_decal_dunes, problem=problem_decal_dunes)
