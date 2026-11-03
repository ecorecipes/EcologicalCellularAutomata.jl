# Teixeira, Horstman & Wijnberg (2023): the growth of a coastal foredune by wind-blown sand, in the
# DuBeVeg model (Keijsers et al. 2016), a descendant of Werner's (1995) sand-slab automaton with
# vegetation (Nield & Baas 2008). The surface is a stack of slabs on each 1 m cell. Each aeolian
# iteration, slabs are picked up and hop downwind (landward) until they land, more readily where
# vegetation traps them; avalanches keep slopes below the angle of repose. Every second iteration the
# sea resets cells below the water level towards a reference beach profile, and once a year two kinds of
# plant grow or die back with the sand they gain or lose. The source's change is to hop a saltation
# length of 5 cells with thinner slabs, letting vegetation trap slabs part-way through a hop.
#
# Slabs travel any distance along a row, so the dynamics are a callback stored with the problem; the
# model's rule leaves cells as they are.

# Piecewise-linear growth of a plant with vertices `v` (m/yr of burial) and peak growth `peak`: -1 at
# and beyond the outer vertices, 0 at the second and fourth, `peak` at the third.
function _dune_growth(Δ, v, peak)
    a, b, c, d, e = v
    (Δ <= a || Δ >= e) && return -1.0
    Δ < b && return -1 + (Δ - a) / (b - a)
    Δ < c && return peak * (Δ - b) / (c - b)
    Δ < d && return peak * (d - Δ) / (d - c)
    return -(Δ - d) / (e - d)
end

"""
    dune_runup(η, slope) -> Float64

The wave run-up (m) above a water level `η` (m) on a beach of gradient `slope`,
by Stockdon et al.'s relation as DuBeVeg applies it, with the wave height and
length correlated with the water level: `H = max(0, 2.931 η - 2.637)`,
`L = max(0, 46.74 η - 30.59)`.
"""
function dune_runup(η, slope)
    H = max(0.0, -2.637 + 2.931η)
    L = max(0.0, -30.59 + 46.74η)
    (H > 0 && L > 0) || return 0.0
    ξ = slope / sqrt(H / L)
    return ξ < 0.3 ? 0.043 * sqrt(H * L) :
           1.1 * (0.35 * slope * sqrt(H * L) + sqrt(H * L * (0.563 * slope^2 + 0.004)) / 2)
end

"""
    build_dune_formation() -> CAModel

The beach-dune model of [Teixeira2023](@cite) (DuBeVeg). Each cell of 1 m holds
an elevation in slabs of height `h` (`topo`), the cover of a pioneer and of a
stabilizing plant (`pioneer`, `stabilizer`; their sum, at most 1, is the
vegetation effectiveness `ρ`), the slabs gained or lost this year
(`balance`), the elevation of the reference beach surface (`zref`), the
probability that the sea has reset the cell this year (`sea`), and the
elevation of a non-erodible surface (`hard`: a seawall or a promenade; none by
default), below which neither wind, avalanches nor the sea remove sand and on
which no plants grow until sand covers it. Rows
run cross-shore, from the sea at the first column to the land at the last;
the wind blows landward. The dynamics are [`dune_step!`](@ref), a callback
stored with the problem.

Beyond the source, `speedup > 0` lets the topography steer the wind
([`dune_wind`](@ref)): the erosion probability scales with the shear stress
along the wind, and flying slabs turn alongshore with the shear across it.
The default, 0, is the source's wind, uniform and straight.
"""
function build_dune_formation()
    rule = LocalRule(Offsets((1, 0))) do c, nbrs, ctx
        return c
    end
    states = [StateVariable(:topo, Int32(0); description="elevation (slabs)"),
              StateVariable(:pioneer, 0.0; description="cover of the pioneer plant"),
              StateVariable(:stabilizer, 0.0; description="cover of the stabilizing plant"),
              StateVariable(:balance, Int32(0);
                            description="slabs gained (lost if negative) this year"),
              StateVariable(:zref, Int32(0);
                            description="elevation of the reference beach surface (slabs)"),
              StateVariable(:sea, 0.0;
                            description="probability the sea has reset the cell this year"),
              StateVariable(:hard, typemin(Int32);
                            description="elevation of a non-erodible surface, such as a seawall (slabs)")]
    parameters = [Parameter(:h, 0.02; description="slab height (m)"),
                  Parameter(:L, 5; description="jumping length (cells)"),
                  Parameter(:rho_star, 0.3;
                            description="vegetation that traps slabs part-way through a hop"),
                  Parameter(:pe0, 0.5; description="basic erosion probability"),
                  Parameter(:pd0, 0.1; description="basic deposition probability"),
                  Parameter(:shadow, 15.0; description="shadow angle (degrees)"),
                  Parameter(:speedup, 0.0;
                            description="strength A of the wind's response to the topography (4 in Kroy et al. 2002); 0 for a wind the topography does not steer"),
                  Parameter(:upwind_shift, 0.25;
                            description="B, which moves the strongest wind upwind of a crest"),
                  Parameter(:repose_bare, 30.0; description="angle of repose of bare sand"),
                  Parameter(:repose_veg, 35.0;
                            description="angle of repose of vegetated sand (ρ at least 0.3)"),
                  Parameter(:G, 0.8; description="groundwater depth factor"),
                  Parameter(:Rv, 0.5;
                            description="resistance of vegetation to wave erosion"),
                  Parameter(:WE, 1.0; description="maximum erosive strength of waves"),
                  Parameter(:Wdiss, 0.012; description="wave dissipation strength"),
                  Parameter(:Pinun, 0.0; description="probability of reset by inundation"),
                  Parameter(:wave_min, 0.1;
                            description="least wave reset probability in exposed wet cells"),
                  Parameter(:depth_min, 0.01;
                            description="least water depth for dissipation (m)"),
                  Parameter(:iterations, 50; description="aeolian iterations a year"),
                  Parameter(:hydro_every, 2;
                            description="aeolian iterations between hydrodynamic steps"),
                  Parameter(:water_levels, [1.0];
                            description="highest water level of each hydrodynamic step (m)"),
                  Parameter(:slope, 0.01; description="beach gradient for wave run-up"),
                  Parameter(:references, Vector{Int32}[];
                            description="reference surface (slabs, sea to land) at the end of each year, replacing zref in every row; empty to keep it fixed"),
                  Parameter(:pioneer_vertices, (-1.4, 0.2, 0.6, 2.0, 2.2);
                            description="burial rates (m/yr) shaping pioneer growth"),
                  Parameter(:pioneer_peak, 0.2;
                            description="greatest yearly growth of pioneers"),
                  Parameter(:stabilizer_vertices, (-1.4, -0.65, 0.0, 0.2, 2.8);
                            description="burial rates (m/yr) shaping stabilizer growth"),
                  Parameter(:stabilizer_peak, 0.05;
                            description="greatest yearly growth of stabilizers"),
                  Parameter(:p_lateral, 0.2;
                            description="yearly probability a bare cell next to plants is colonised"),
                  Parameter(:p_establish, 0.05;
                            description="yearly probability a bare cell is colonised anyway")]
    return CAModel(rule; name=:dune_formation, states, parameters,
                   description="Teixeira et al. (2023) beach-dune evolution with saltation hops")
end

_dune_rho(u, I) = min(1.0, u.pioneer[I] + u.stabilizer[I])

# The non-erodible surface at a cell; a state without a `hard` layer has none.
_dune_hard(u, I) = hasproperty(u, :hard) ? u.hard[I] : typemin(Int32)

# Cells in the lee of upwind topography: some upwind cell i cells away stands more than
# i tan(shadow) above it.
function _dune_shadow!(S, topo, p)
    s = tand(p.shadow) / p.h
    for r in axes(topo, 1)
        M = -Inf
        for x in axes(topo, 2)
            S[r, x] = M > topo[r, x]
            M = max(M, Float64(topo[r, x])) - s
        end
    end
    return S
end

"""
    dune_wind(u, p) -> NamedTuple

The surface shear stress of the wind over the state `u`, relative to its value
over a flat bed, by the linear theory of turbulent flow over low hills
(Jackson & Hunt; Weng et al.) in the form of [Kroy2002](@cite) along the wind
and [Schwammle2005](@cite) across it: with `ĥ` the Fourier transform of the
height above the reference surface `zref` (m) and `kx`, `ky` the wave numbers
along and across the wind (per m),

    along  = 1 + F⁻¹{A kx (kx + i B |kx|) / |k| ĥ}
    across =     F⁻¹{A kx ky / |k| ĥ}

where `A = p.speedup` and `B = p.upwind_shift`. The wind speeds up over crests
and windward slopes, most just upwind of a crest, slows at the foot of a
slope, and turns towards the normal of a ridge it crosses at an angle: into a
gap through a ridge, out of it beyond. In the lee the theory is applied to the
envelope of the topography and the shadow lines (15° by default, `p.shadow`),
as in [Kroy2002](@cite). The domain is mirrored alongshore and, cross-shore,
extended by a smooth return from the land edge to the sea edge, so the
transform sees no jumps. Returns `along` and `across` (rows by columns; rows
run alongshore, increasing `across` points to higher rows).
"""
function dune_wind(u, p)
    m, n = size(u.topo)
    s = tand(p.shadow) / p.h
    hgt = zeros(2m, 2n)
    for r in 1:m
        M = -Inf
        for x in 1:n
            M = max(M, Float64(u.topo[r, x]))
            hgt[r, x] = (M - u.zref[r, x]) * p.h
            M -= s
        end
        a, b = hgt[r, n], hgt[r, 1]
        for x in (n + 1):(2n)
            hgt[r, x] = a + (b - a) * (1 - cospi((x - n) / (n + 1))) / 2
        end
    end
    for r in 1:m
        hgt[2m + 1 - r, :] .= view(hgt, r, :)
    end
    τx, τy = _dune_shear(hgt, p.speedup, p.upwind_shift)
    return (; along=1 .+ τx[1:m, 1:n], across=τy[1:m, 1:n])
end

# The shear stress perturbations along (second index) and across (first index) the wind over a
# periodic height field on 1 m cells.
function _dune_shear(hgt, A, B)
    F = FFTW.fft(hgt)
    kx = 2π .* FFTW.fftfreq(size(hgt, 2))
    ky = 2π .* FFTW.fftfreq(size(hgt, 1))
    Fx = similar(F)
    Fy = similar(F)
    for j in axes(F, 2), i in axes(F, 1)
        k = hypot(kx[j], ky[i])
        if k == 0
            Fx[i, j] = Fy[i, j] = 0
        else
            Fx[i, j] = A * kx[j] * (kx[j] + im * B * abs(kx[j])) / k * F[i, j]
            Fy[i, j] = A * kx[j] * ky[i] / k * F[i, j]
        end
    end
    return real.(FFTW.ifft(Fx)), real.(FFTW.ifft(Fy))
end

# One aeolian iteration of slab transport; returns the slabs crossing column `toe[r]` in each row.
# With a `wind` (the shear stress along and across the wind, relative to a flat bed), the erosion
# probability scales with the shear along the wind (at most 1), and at each cell a flying slab
# moves to the next row with probability |across| / along (at most 1), towards the side `across`
# points to; it stays in its row at the domain's edges.
function _dune_transport!(u, p, rng, shade, toe, wind=nothing)
    topo = u.topo
    m, n = size(topo)
    pe = zeros(m, n)
    pdend = zeros(m, n)
    pdmid = zeros(m, n)
    for I in CartesianIndices(topo)
        ρ = _dune_rho(u, I)
        dry = topo[I] > 0 && topo[I] > p.G * u.zref[I]
        loose = topo[I] > _dune_hard(u, I)
        pe[I] = shade[I] || !dry || !loose ? 0.0 : (1 - ρ) * p.pe0
        wind === nothing || (pe[I] = min(1.0, pe[I] * max(wind.along[I], 0.0)))
        pd = ρ * (1 - p.pd0) + p.pd0
        pdend[I] = shade[I] ? 1.0 : pd
        pdmid[I] = ρ >= p.rho_star ? pd : 0.0
    end
    change = zeros(Int32, m, n)
    crossed = zeros(Int, m)
    for x in 1:n, r0 in 1:m
        rand(rng) < pe[r0, x] || continue
        change[r0, x] -= 1
        r = r0
        k = 0
        while true
            k += 1
            y = x + k
            if wind !== nothing
                a, c = wind.along[r, y - 1], wind.across[r, y - 1]
                s = c > 0 ? 1 : -1
                a > 0 && rand(rng) < abs(c) / a && 1 <= r + s <= m && (r += s)
            end
            if y > n
                toe[r] == n + 1 && (crossed[r] += 1)
                break
            end
            toe[r] == y && (crossed[r] += 1)
            q = k % p.L == 0 ? pdend[r, y] : pdmid[r, y]
            if rand(rng) < q
                change[r, y] += 1
                break
            end
        end
    end
    topo .+= change
    u.balance .+= change
    return crossed
end

# Avalanches: until no slope exceeds the angle of repose, each cell steeper than it towards its
# steepest lower neighbour (of 8) passes it one slab, all at once.
function _dune_avalanche!(u, p, rng)
    topo = u.topo
    m, n = size(topo)
    lim(I, diag) = floor(Int,
                         (diag ? sqrt(2) : 1.0) *
                         tand(_dune_rho(u, I) >= 0.3 ? p.repose_veg : p.repose_bare) / p.h)
    moves = zeros(Int32, m, n)
    for _ in 1:10_000
        fill!(moves, 0)
        moved = false
        for I in CartesianIndices(topo)
            topo[I] > _dune_hard(u, I) || continue
            best, steep, ties = I, 0, 0
            for di in -1:1, dj in -1:1
                (di, dj) == (0, 0) && continue
                J = I + CartesianIndex(di, dj)
                checkbounds(Bool, topo, J) || continue
                excess = topo[I] - topo[J] - lim(I, di != 0 && dj != 0)
                excess > 0 || continue
                if excess > steep
                    best, steep, ties = J, excess, 1
                elseif excess == steep
                    ties += 1
                    rand(rng) < 1 / ties && (best = J)
                end
            end
            if steep > 0
                moves[I] -= 1
                moves[best] += 1
                moved = true
            end
        end
        moved || break
        topo .+= moves
        u.balance .+= moves
    end
    return u
end

# The sea: cells below the total water level and seaward of the first dry cell are reset to the
# reference surface with probability (1 - Rv ρ)(wave + Pinun); the wave strength falls with the
# shallow water crossed.
function _dune_hydro!(u, p, rng, η)
    topo = u.topo
    m, n = size(topo)
    level = (η + dune_runup(η, p.slope)) / p.h
    for r in 1:m
        D = 0.0
        for x in 1:n
            topo[r, x] >= level && break
            depth = max((level - topo[r, x]) * p.h, p.depth_min)
            wave = clamp(p.WE - p.Wdiss * D, p.wave_min, 1.0)
            P = min(1.0, (1 - p.Rv * _dune_rho(u, CartesianIndex(r, x))) * (wave + p.Pinun))
            u.sea[r, x] = 1 - (1 - u.sea[r, x]) * (1 - P)
            target = max(u.zref[r, x], _dune_hard(u, CartesianIndex(r, x)))
            if rand(rng) < P && topo[r, x] != target
                d = target - topo[r, x]
                topo[r, x] += d
                u.balance[r, x] += d
                u.pioneer[r, x] = 0.0
                u.stabilizer[r, x] = 0.0
            end
            D += 1 / depth
        end
    end
    return u
end

# A year of plant growth, from each cell's burial (m) over the year.
function _dune_vegetation!(u, p, rng)
    Δ = u.balance .* p.h
    veg = (u.pioneer .+ u.stabilizer) .> 0
    for I in CartesianIndices(u.topo)
        if u.topo[I] <= _dune_hard(u, I)
            u.pioneer[I] = 0.0
            u.stabilizer[I] = 0.0
            continue
        end
        vegetated = veg[I]
        if !vegetated
            u.topo[I] > 0 || continue
            near = any(J -> checkbounds(Bool, veg, J) && veg[J],
                       (I + CartesianIndex(di, dj)
                        for di in -1:1, dj in -1:1
                        if (di, dj) != (0, 0)))
            ((near && rand(rng) < p.p_lateral) || rand(rng) < p.p_establish) || continue
        end
        g1 = _dune_growth(Δ[I], p.pioneer_vertices, p.pioneer_peak)
        g2 = _dune_growth(Δ[I], p.stabilizer_vertices, p.stabilizer_peak)
        if vegetated
            u.pioneer[I] = clamp(u.pioneer[I] + g1, 0.0, 1.0)
            u.stabilizer[I] = clamp(u.stabilizer[I] + g2, 0.0, 1.0)
        else
            u.pioneer[I] = clamp(g1, 0.0, 1.0)
            u.stabilizer[I] = Δ[I] == 0 ? clamp(g2, 0.0, 1.0) : 0.0
        end
    end
    u.pioneer .*= 1 .- u.sea
    u.stabilizer .*= 1 .- u.sea
    u.sea .= 0
    u.balance .= 0
    return u
end

"""
    dune_toe(u, r; level=3, h=0.02) -> Int

The source's dune toe of row `r`: the most seaward column whose elevation is
above `level` m, with slabs of height `h`, or one past the last column if none
is.
"""
function dune_toe(u, r::Integer; level::Real=3, h::Real=0.02)
    x = findfirst(>(level / h), view(u.topo, r, :))
    return x === nothing ? size(u.topo, 2) + 1 : x
end

"""
    dune_step!(u, p, t, rng) -> Int

One aeolian iteration `t` of [Teixeira2023](@cite), applied to the state `u` in
place: the shadow zones and the probabilities of erosion and deposition are
fixed from the state at its start; every cell above the groundwater surface
(`G` times the reference surface) and above sea level, outside the shadow,
gives up a slab with probability `(1 - ρ) pe0`; each slab moves landward a cell
at a time, landing at the end of every hop of `L` cells with probability
`ρ (1 - pd0) + pd0` (1 in the shadow) and, part-way through a hop, with the
same probability where `ρ ≥ rho_star` (never otherwise); slabs leaving the
land side are lost. Avalanches then restore the angles of repose. Every
`hydro_every` iterations the sea acts with the next water level of
`water_levels` plus wave run-up ([`dune_runup`](@ref)), followed by
avalanches; every `iterations` iterations the plants grow, and their cover is
then multiplied by the chance that the sea did not reset the cell during the
year, as in DuBeVeg, and, if `references` is not empty, the reference surface
becomes that of the year just ended (a prograding beach). Returns the slabs
that crossed the dune toe (the first cell above 3 m) in all rows.

With `speedup > 0` (beyond the source), the shear stress of
[`dune_wind`](@ref) is computed from the state at the start of the iteration:
the erosion probability is multiplied by the shear along the wind (and
capped at 1), and at each cell a flying slab moves one row alongshore, towards
the side the shear across the wind points to, with probability
`|across| / along` (at most 1; never off the domain).
"""
function dune_step!(u, p, t::Integer, rng)
    m, n = size(u.topo)
    shade = falses(m, n)
    _dune_shadow!(shade, u.topo, p)
    toe = [dune_toe(u, r; h=p.h) for r in 1:m]
    wind = p.speedup == 0 ? nothing : dune_wind(u, p)
    crossed = _dune_transport!(u, p, rng, shade, toe, wind)
    _dune_avalanche!(u, p, rng)
    if t % p.hydro_every == 0
        η = p.water_levels[mod1(t ÷ p.hydro_every, length(p.water_levels))]
        _dune_hydro!(u, p, rng, η)
        _dune_avalanche!(u, p, rng)
    end
    if t % p.iterations == 0
        _dune_vegetation!(u, p, rng)
        if !isempty(p.references)
            row = p.references[min(t ÷ p.iterations + 1, length(p.references))]
            for r in 1:m
                u.zref[r, :] .= row
            end
        end
    end
    return sum(crossed)
end

"""
    dune_callback() -> DiscreteCallback

A callback that runs [`dune_step!`](@ref) after each step, with a random
number generator seeded from the problem's seed and the iteration.
"""
function dune_callback()
    affect! = function (integ)
        t = Int(integ.t)
        rng = Random.Xoshiro(hash((integ.prob.seed, t, :dune)))
        dune_step!(integ.u, integ.p, t, rng)
        return nothing
    end
    return DiscreteCallback((u, t, integ) -> true, affect!; save_positions=(false, false))
end

"""
    dune_landform(; form=:foredune, beach=500, offshore=60, toe=2.5, berm=0, berm_level=2,
                  embryo=40, embryo_crest=7, slack=160, wall=7.5, breach=5) -> ProfileSegment

A cross-shore profile as joined pieces ([`ProfileSegment`](@ref)), from the sea
to the land, `offshore + beach + 120` m long. With `form = :foredune`, the
source's: a Dean shoreface over `offshore` m ([`shoreface_segment`](@ref)), a
beach rising linearly from 0 at the shoreline to `toe` m at the foot of the
dune over `beach` m (2.5 m, as the source's figure shows; its text says 3 m),
a foredune ridge 100 m wide rising to 15 m (a smooth bump, its shape read from
the source's figure; [`ridge_segment`](@ref)), and a 5 m plateau behind it.
With `berm > 0` the beach is cut where it reaches `berm_level` m and a level
berm of `berm` m is inserted there: the beach of the source's prograding
scenario after it has built out by `berm` m. Berms add: the landform with
`berm = a + b` is the one with `berm = a` with `b` m more of berm joined in.

Other forms, which the source does not run, are composed of the same pieces
(and have the same length):

- `:double_ridge`: an embryo dune `embryo` m wide rising to `embryo_crest` m
  stands on the beach, its land side `slack` m seaward of the foredune, with
  a level dune slack (kind `:slack`) at `toe` m between them.
- `:seawall`: the beach ends at a wall (kind `:seawall`, 2 m, rising from
  `toe` to `wall` m) with a level promenade (kind `:promenade`) behind it, in
  place of the dune and plateau.
- `:breach`: the foredune is lowered to a trough rising from `toe` to `breach`
  m before falling to the plateau: a blowout through the ridge, for the rows
  it cuts.
"""
function dune_landform(; form::Symbol=:foredune, beach::Integer=500, offshore::Integer=60,
                       toe::Real=2.5, berm::Integer=0, berm_level::Real=2,
                       embryo::Integer=40, embryo_crest::Real=7, slack::Integer=160,
                       wall::Real=7.5, breach::Real=5)
    beach >= 10 || throw(ArgumentError("the beach must be at least 10 m wide, got $beach"))
    berm >= 0 || throw(ArgumentError("berm must be nonnegative, got $berm"))
    form in (:foredune, :double_ridge, :seawall, :breach) ||
        throw(ArgumentError("form must be :foredune, :double_ridge, :seawall or :breach, got :$form"))
    form === :foredune || berm == 0 ||
        throw(ArgumentError("a berm is inserted only in the :foredune form"))
    shore = shoreface_segment(offshore)
    plateau = flat_segment(20, 5; kind=:hinterland)
    if form === :double_ridge
        open = beach - embryo - slack
        open >= 10 ||
            throw(ArgumentError("the beach of $beach m leaves $open m seaward of an embryo dune of $embryo m with a slack of $slack m; at least 10 m are needed"))
        return shore * slope_segment(open, 0, toe) *
               ridge_segment(embryo, toe, embryo_crest, toe) *
               flat_segment(slack, toe; kind=:slack) * ridge_segment(100, toe, 15, 5) *
               plateau
    elseif form === :seawall
        wall > toe ||
            throw(ArgumentError("the wall top $wall m must be above the toe $toe m"))
        return shore * slope_segment(beach, 0, toe) *
               slope_segment(2, toe, wall; kind=:seawall) *
               flat_segment(118, wall; kind=:promenade)
    elseif form === :breach
        return shore * slope_segment(beach, 0, toe) * ridge_segment(100, toe, breach, 5) *
               plateau
    end
    dune = ridge_segment(100, toe, 15, 5) * plateau
    berm == 0 && return shore * slope_segment(beach, 0, toe) * dune
    0 < berm_level < toe ||
        throw(ArgumentError("berm_level must lie between 0 and the dune foot $toe m, got $berm_level"))
    lower = round(Int, beach * berm_level / toe)
    return shore * slope_segment(lower, 0, berm_level) * flat_segment(berm, berm_level) *
           slope_segment(beach - lower, berm_level, toe) * dune
end

"""
    dune_profile(; beach=500, offshore=60, toe=2.5, h=0.02) -> NamedTuple

The source's initial cross-shore profile ([`dune_landform`](@ref)), in slabs
of height `h` on 1 m cells from the sea to the land: its `elevation` and its
`reference` surface, the beach without the dune, level under it
([`reference_profile`](@ref)).
"""
function dune_profile(; beach::Integer=500, offshore::Integer=60, toe::Real=2.5,
                      h::Real=0.02)
    p = dune_landform(; beach, offshore, toe)
    n = length(p)
    return (elevation=profile_heights(p, n; h),
            reference=profile_heights(reference_profile(p), n; h))
end

"""
    dune_water_levels(rng; years=50, level=1.43, surge=0.24) -> Vector{Float64}

A synthetic series of the highest water level (m above mean sea level) in
each spring-neap cycle, 25 a year, standing in for the source's observed
levels at Scheveningen: `level` m plus a Gumbel-distributed surge of scale
`surge` m. The defaults are fitted to the series the source plots (its
Fig. 5d, read from the figure): a median of about 1.5 m, a tenth of the
cycles above 2 m and a highest level of about 3.2 m in 50 years.
"""
function dune_water_levels(rng; years::Integer=50, level::Real=1.43, surge::Real=0.24)
    surge >= 0 || throw(ArgumentError("surge must be nonnegative, got $surge"))
    return [level - surge * log(-log(rand(rng)))
            for _ in 1:(25years)]
end

"""
    problem_dune_formation(; scenario=:wide, representation=:saltation, alongshore=100,
                           years=50, rate=4, landform=nothing, seed=0, pmap=(;),
                           kwargs...) -> CAProblem

The setup of [Teixeira2023](@cite) for its `:narrow` (100 m) or `:wide`
(500 m) beach (see [`dune_profile`](@ref)), or its `:prograding` beach: the
narrow beach on a domain extending `rate × years` m further out to sea, whose
reference surface builds out by `rate` m a year (a berm at 2 m growing by
`rate` m, the profile below it moving seaward, [`dune_landform`](@ref)),
anchored at the land; the sea moves the beach towards it. With `alongshore` rows of 1 m,
for `years` years of 50 aeolian iterations. `representation = :saltation` is
the source's change (slabs of 0.02 m hopping 5 cells, trapped part-way by
vegetation of at least 0.3); `:default` is DuBeVeg as it was (slabs of
0.1 m hopping 1 cell). Plants cover the foredune fully above 6 m and with a
random cover between 4 and 6 m, pioneers seaward of each row's highest cell and
stabilizers landward of it (the source does not say which). Water levels
come from [`dune_water_levels`](@ref), drawn from `seed`.

Three more scenarios run landforms the source does not (see
[`dune_landform`](@ref)): `:double_ridge`, the wide beach with an embryo dune
seaward of the foredune; `:seawall`, the narrow beach ending at a seawall and
promenade, which neither wind, avalanches nor the sea can erode (the `hard`
layer); and `:blowout`, the narrow beach with the foredune breached in the
middle fifth of the rows.

`landform`, a [`ProfileSegment`](@ref) or a vector of them (one for each row,
all of the same length), replaces the scenario's profile: its reference
surface is [`reference_profile`](@ref)'s, pieces of kind `:seawall` or
`:promenade` are non-erodible, and the beach gradient for wave run-up is that
of its `:beach` pieces (rise over length) unless `pmap` gives `slope`.
"""
function problem_dune_formation(; scenario::Symbol=:wide, representation::Symbol=:saltation,
                                alongshore::Integer=100, years::Integer=50, rate::Integer=4,
                                landform=nothing, seed::Integer=0, pmap=(;), kwargs...)
    scenario in (:narrow, :wide, :prograding, :double_ridge, :seawall, :blowout) ||
        throw(ArgumentError("scenario must be :narrow, :wide, :prograding, :double_ridge, :seawall or :blowout, got :$scenario"))
    rate >= 0 || throw(ArgumentError("rate must be nonnegative, got $rate"))
    representation in (:saltation, :default) ||
        throw(ArgumentError("representation must be :saltation or :default, got :$representation"))
    rep = representation === :saltation ? (; h=0.02, L=5, rho_star=0.3) :
          (; h=0.1, L=1, rho_star=Inf)
    rng = Random.Xoshiro(hash((seed, 0, :dune)))
    beach = scenario in (:wide, :double_ridge) ? 500 : 100
    offshore = scenario === :prograding ? 60 + rate * years : 60
    slope = 3 / beach
    rows = if landform !== nothing
        forms = landform isa ProfileSegment ? fill(landform, alongshore) : collect(landform)
        length(forms) == alongshore ||
            throw(ArgumentError("landform gives $(length(forms)) rows, alongshore is $alongshore"))
        allequal(length.(forms)) ||
            throw(ArgumentError("the landforms of the rows differ in length: $(unique(length.(forms)))"))
        haskey(pmap, :slope) || (slope = _beach_gradient(first(forms)))
        forms
    elseif scenario === :blowout
        gap = (alongshore * 2 ÷ 5 + 1):(alongshore * 3 ÷ 5)
        [r in gap ? dune_landform(; form=:breach, beach, offshore) :
         dune_landform(; beach, offshore) for r in 1:alongshore]
    else
        form = scenario === :double_ridge ? :double_ridge :
               scenario === :seawall ? :seawall : :foredune
        fill(dune_landform(; form, beach, offshore), alongshore)
    end
    n = length(first(rows))
    references = scenario === :prograding && landform === nothing ?
                 [profile_heights(reference_profile(dune_landform(; beach, offshore,
                                                                  berm=rate * y)), n;
                                  h=rep.h) for y in 0:years] : Vector{Int32}[]
    given = merge(rep, (; water_levels=dune_water_levels(rng; years), slope, references),
                  (; (Symbol(k) => v for (k, v) in pairs(pmap))...))
    topo = permutedims(reduce(hcat, [profile_heights(f, n; h=rep.h) for f in rows]))
    zref = permutedims(reduce(hcat,
                              [profile_heights(reference_profile(f), n; h=rep.h)
                               for f in rows]))
    hard = fill(typemin(Int32), alongshore, n)
    for (r, f) in enumerate(rows), (x, k) in enumerate(profile_kinds(f, n))
        k in (:seawall, :promenade) && (hard[r, x] = topo[r, x])
    end
    pioneer = zeros(alongshore, n)
    stabilizer = zeros(alongshore, n)
    crest = [argmax(view(topo, r, :)) for r in 1:alongshore]
    for I in CartesianIndices(topo)
        topo[I] > hard[I] || continue
        z = topo[I] * rep.h
        cover = z > 6 ? 1.0 : z >= 4 ? rand(rng) : 0.0
        cover > 0 || continue
        I[2] <= crest[I[1]] ? (pioneer[I] = cover) : (stabilizer[I] = cover)
    end
    model = build_dune_formation()
    return CAProblem(model, (; topo, pioneer, stabilizer, zref, hard), 50years, given;
                     space=Grid(alongshore, n), seed, callback=dune_callback(), kwargs...)
end

# The gradient of a profile's beach: the rise of its `:beach` pieces over their length.
function _beach_gradient(p::ProfileSegment)
    beach = filter(s -> s.kind === :beach, _pieces(p))
    run = sum(length, beach; init=0)
    run > 0 ||
        throw(ArgumentError("the landform has no :beach piece to set the run-up slope; give slope in pmap"))
    return sum(s -> s.z1 - s.z0, beach) / run
end

"""
    dune_measures(u; h=0.02) -> NamedTuple

The source's measures of the foredune, averaged over the rows: its `volume`
(m³ per m alongshore, above 3 m), its `width` (m, from the dune toe to the
land side) and the position of its `toe` (the column of the first cell above
3 m), with slabs of height `h`.
"""
function dune_measures(u; h::Real=0.02)
    m, n = size(u.topo)
    toes = [dune_toe(u, r; h) for r in 1:m]
    volume = sum(max.(u.topo .* h .- 3, 0.0)) / m
    return (; volume, width=n + 1 - sum(toes) / m, toe=sum(toes) / m)
end

_register!(:dune_formation; build=build_dune_formation, problem=problem_dune_formation)
