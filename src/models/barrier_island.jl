# Reeves et al. (2025): MEEB v1.0, the Mesoscale Explicit Ecogeomorphic Barrier model, the evolution of an
# undeveloped barrier island over years to decades on 1 m cells and weekly steps. Each week wind-blown sand
# moves as slabs (DuBeVeg with the saltation of Teixeira et al. 2023) in one of four directions drawn from a
# wind rose; every second week a high-water event may strike, eroding the beach and dune face seaward of the
# foredune crest (Larson et al. 2004; Duran Vinent and Moore 2015) and routing overwash and its sand across
# the barrier landward of it (Barrier3D), killing the plants it floods; the ocean shoreline then moves with
# the shoreface (Lorenzo-Trueba and Ashton 2014) and alongshore diffusion by the wave climate (Ashton and
# Murray 2006), slopes are relaxed to the angle of repose, and once a year two plant types grow or die back
# with their burial. Where the paper is brief, the implementation follows the authors' code (MEEB 1.0.1).
#
# Slabs travel any distance and storms route water across whole rows, so the dynamics are a callback
# stored with the problem; the model's rule leaves cells as they are. The cells are 1 m, the paper's.

# ---------------------------------------------------------------------------------------------------
# Filters, as in scipy.ndimage and scipy.signal (the reflect mode mirrors about the edge, d c b a | a b c d)

_meeb_reflect(j, n) = (k=mod(j - 1, 2n); k < n ? k + 1 : 2n - k)

# A Gaussian filter of standard deviation σ cells, truncated at 4σ, along both axes; `mode` is
# `:reflect` or `:constant` (zero outside).
function _meeb_gaussian(A::AbstractMatrix, σ::Real; mode::Symbol=:reflect)
    r = floor(Int, 4σ + 0.5)
    w = [exp(-0.5 * (k / σ)^2) for k in (-r):r]
    w ./= sum(w)
    m, n = size(A)
    B = zeros(m, n)
    for j in 1:n, i in 1:m
        s = 0.0
        for k in (-r):r
            ii = i + k
            if 1 <= ii <= m
                s += w[k + r + 1] * A[ii, j]
            elseif mode === :reflect
                s += w[k + r + 1] * A[_meeb_reflect(ii, m), j]
            end
        end
        B[i, j] = s
    end
    C = zeros(m, n)
    for j in 1:n, i in 1:m
        s = 0.0
        for k in (-r):r
            jj = j + k
            if 1 <= jj <= n
                s += w[k + r + 1] * B[i, jj]
            elseif mode === :reflect
                s += w[k + r + 1] * B[i, _meeb_reflect(jj, n)]
            end
        end
        C[i, j] = s
    end
    return C
end

# The mean over a window of `w` rows (alongshore), reflected at the edges, as scipy's uniform_filter1d:
# the window of row i runs from i - w ÷ 2 to i + w - w ÷ 2 - 1.
function _meeb_alongshore_mean(A::AbstractMatrix, w::Integer)
    m, n = size(A)
    lo = w ÷ 2
    B = zeros(m, n)
    for j in 1:n, i in 1:m
        s = 0.0
        for k in (i - lo):(i + w - lo - 1)
            s += A[_meeb_reflect(k, m), j]
        end
        B[i, j] = s / w
    end
    return B
end

# A Savitzky-Golay filter of polynomial order 1 over an odd window of w points (scipy's `interp`
# mode): the mean of the window in the interior, and at each end the least-squares line through the
# first or last w points. A window longer than the series is shortened to it, and an even one by one.
function _meeb_savgol1(x::AbstractVector, w::Integer)
    n = length(x)
    w = min(w, n)
    iseven(w) && (w -= 1)
    y = float.(x)
    w >= 3 || return y
    h = w ÷ 2
    for i in (h + 1):(n - h)
        y[i] = sum(view(x, (i - h):(i + h))) / w
    end
    ts = 0:(w - 1)
    tm = sum(ts) / w
    line(seg, t) = begin
        sm = sum(seg) / w
        b = sum((ts .- tm) .* (seg .- sm)) / sum((ts .- tm) .^ 2)
        sm + b * (t - tm)
    end
    head = float.(x[1:w])
    tail = float.(x[(n - w + 1):n])
    for i in 1:h
        y[i] = line(head, i - 1)
        y[n - h + i] = line(tail, w - h + i - 1)
    end
    return y
end

# ---------------------------------------------------------------------------------------------------
# Shorelines and the foredune crest

# The ocean shoreline of row r: its first cell at or above MHW (the first cell if none is).
function _meeb_ocean_shoreline(topo, r, MHW)
    x = findfirst(>=(MHW), view(topo, r, :))
    return x === nothing ? 1 : x
end

# The back-barrier shoreline of row r: its last cell at or above MHW (the last cell if none is).
function _meeb_backbarrier(topo, r, MHW)
    x = findlast(>=(MHW), view(topo, r, :))
    return x === nothing ? size(topo, 2) : x
end

# The columns from the most seaward ocean shoreline to the cell before the most landward
# back-barrier shoreline, where the wind and avalanches act.
function _meeb_window(topo, MHW)
    m = size(topo, 1)
    lo = minimum(_meeb_ocean_shoreline(topo, r, MHW) for r in 1:m)
    hi = maximum(_meeb_backbarrier(topo, r, MHW) for r in 1:m) - 1
    return lo:max(lo, hi)
end

# Peaks of a profile: points higher than the point before and the first different point after.
function _meeb_local_maxima(p)
    pks = Int[]
    n = length(p)
    r = 2
    while r < n
        back, cur = p[r - 1], p[r]
        f = r + 1
        while f <= n && p[f] == cur
            f += 1
        end
        f > n && break
        back < cur && p[f] < cur && push!(pks, r)
        r = f
    end
    return pks
end

# The foredune peak of a profile (AutoMorph, Itzkin et al. 2021): the first peak above MHW with a drop
# of `threshold` m behind it before the ground rises above it, moved to the most seaward peak higher
# than (1 - crest_pct) times its elevation. Returns the index (or `nothing` without peaks above MHW)
# and whether a peak with that drop was found; without one, the highest point.
function _meeb_find_crest(p, MHW; threshold=0.6, crest_pct=0.1)
    pks = filter(i -> p[i] > MHW, _meeb_local_maxima(p))
    isempty(pks) && return nothing, false
    for i in pks
        drop = 0.0
        c = i
        while c + 1 <= length(p) && p[c + 1] <= p[i]
            drop = p[i] - p[c + 1]
            drop >= threshold && break
            c += 1
        end
        if drop >= threshold
            lo = p[i] * (1 - crest_pct)
            near = filter(k -> p[k] > lo, pks)
            return (isempty(near) ? i : first(near)), true
        end
    end
    return argmax(p), false
end

"""
    barrier_crest(topo, MHW) -> (crest, not_gap)

The cross-shore column of the foredune crest in each row of `topo` (m), by
MEEB's procedure [Reeves2025](@cite): the highest point of the topography averaged
over 150 m alongshore, within the first 100 m of a narrow barrier (under 300 m)
or its seaward 30% otherwise, smoothed alongshore (a Savitzky-Golay filter over
75 m); then in each row the most seaward peak within 25 m of that line with a
drop of 0.6 m behind it (AutoMorph's test); and that line smoothed over 11 m.
`not_gap` is false in rows with no such peak: gaps in the foredune.
"""
function barrier_crest(topo::AbstractMatrix, MHW::Real)
    m, n = size(topo)
    xs = [_meeb_ocean_shoreline(topo, r, MHW) for r in 1:m]
    xb = [_meeb_backbarrier(topo, r, MHW) for r in 1:m]
    lo, hi = minimum(xs), maximum(xb) - 1
    hi < lo && return fill(lo, m), zeros(Bool, m)
    crop = view(topo, :, lo:hi)
    W = size(crop, 2)
    avg = _meeb_alongshore_mean(crop, min(150, m))
    line = zeros(m)
    for r in 1:m
        s, b = xs[r] - lo, xb[r] - lo           # 0-based within the crop
        stop = b - s < 300 ? max(s + 100, b) : s + floor(Int, (b - s) * 0.3)
        stop = min(stop, W)
        line[r] = stop > s ? s + argmax(view(avg, r, (s + 1):stop)) - 1 : s
    end
    line = round.(Int, _meeb_savgol1(line, min(75, m)))
    crest = copy(line)
    not_gap = ones(Bool, m)
    for r in 1:m
        if line[r] > 0
            mini = max(0, line[r] - 25)
            maxi = min(W - 1, line[r] + 25)
            loc, found = maxi > mini ?
                         _meeb_find_crest(view(crop, r, (mini + 1):maxi), MHW) :
                         (nothing, false)
            found || (not_gap[r] = false)
            crest[r] = loc === nothing ? crest[r == 1 ? m : r - 1] : mini + loc - 1
        else
            crest[r] = argmax(view(crop, r, :)) - 1
        end
    end
    crest = round.(Int, _meeb_savgol1(crest, 11))
    return clamp.(crest .+ lo, 1, n - 1), not_gap
end

# ---------------------------------------------------------------------------------------------------
# Aeolian transport

# Wind directions: 1 landward (increasing column), 2 alongshore down (increasing row), 3 seaward,
# 4 alongshore up.
const _MEEB_STEPS = ((0, 1), (1, 0), (0, -1), (-1, 0))

# Cells in the lee of topography upwind: some cell i cells upwind, within MEEB's search range (the
# highest elevation over tan(shadow)), stands more than i tan(shadow) above it.
function _meeb_shadow(sub, dir, shadow)
    s = tand(shadow)
    m, n = size(sub)
    R = ceil(Int, maximum(sub) / s)
    di, dj = _MEEB_STEPS[dir]
    S = falses(m, n)
    for j in 1:n, i in 1:m
        for k in 1:R
            ii, jj = i - k * di, j - k * dj
            (1 <= ii <= m && 1 <= jj <= n) || break
            if sub[ii, jj] - sub[i, j] > s * k
                S[i, j] = true
                break
            end
        end
    end
    return S
end

# One week of slab transport by wind towards `dir`; returns the change in slabs.
function _meeb_aeolian(u, p, rng, MHW, dir)
    topo = u.topo
    cols = _meeb_window(topo, MHW)
    gw = max.(MHW .+ p.groundwater_depth .* (_meeb_gaussian(topo, 12.0) .- MHW), MHW)
    sub = view(topo, :, cols)
    shade = _meeb_shadow(sub, dir, p.shadow)
    veff = view(u.veg_eff, :, cols)
    gws = view(gw, :, cols)
    m, n = size(sub)
    Pe = zeros(m, n)
    Pd = zeros(m, n)
    for I in CartesianIndices(sub)
        z = sub[I]
        loose = !shade[I] && z > MHW && z > gws[I] && z - MHW > p.h
        Pe[I] = loose ? max(0.0, p.pe0 - p.pe0 / p.rho_entrain * veff[I]) : 0.0
        Pd[I] = shade[I] || z <= gws[I] ? 1.0 : min(1.0, p.pd0 + (p.pd1 - p.pd0) * veff[I])
    end
    hop = rand(rng, (p.L - p.L_dev):(p.L + p.L_dev))
    di, dj = _MEEB_STEPS[dir]
    change = zeros(m, n)
    for I in CartesianIndices(sub)
        rand(rng) < Pe[I] || continue
        change[I] -= 1
        i, j = Tuple(I)
        k = 0
        while true
            k += 1
            i += di
            j += dj
            (1 <= i <= m && 1 <= j <= n) || break
            q = k % hop == 0 || veff[i, j] >= p.rho_saltation ? Pd[i, j] : 0.0
            if rand(rng) < q
                change[i, j] += 1
                break
            end
        end
    end
    full = zeros(size(topo))
    full[:, cols] .= change
    return full
end

# Avalanches (DuBeVeg's rule): until no slab moves, each cell above MHW whose steepest drop to its 8
# neighbours exceeds the angle of repose (in whole slabs, times √2 along a diagonal; vegetated at a
# cover of at least repose_threshold) passes one slab down it, all cells at once, ties at random.
function _meeb_avalanche!(u, p, rng, MHW)
    topo = u.topo
    cols = _meeb_window(topo, MHW)
    T = topo[:, cols] ./ p.h
    vg = view(u.veg, :, cols)
    m, n = size(T)
    lim(I, diag) = floor((diag ? sqrt(2) : 1.0) *
                         tand(vg[I] >= p.repose_threshold ? p.repose_veg : p.repose_bare) /
                         p.h)
    moves = zeros(m, n)
    for _ in 1:100_000
        fill!(moves, 0)
        moved = false
        for I in CartesianIndices(T)
            T[I] * p.h > MHW || continue
            low = 0.0
            for di in -1:1, dj in -1:1
                J = I + CartesianIndex(di, dj)
                checkbounds(Bool, T, J) && (low = min(low, T[J] - T[I]))
            end
            low < 0 || continue
            best, ties = I, 0
            for di in -1:1, dj in -1:1
                (di, dj) == (0, 0) && continue
                J = I + CartesianIndex(di, dj)
                checkbounds(Bool, T, J) || continue
                d = T[J] - T[I]
                if d == low && -d > lim(I, di != 0 && dj != 0)
                    ties += 1
                    rand(rng) < 1 / ties && (best = J)
                end
            end
            if ties > 0
                moves[I] -= 1
                moves[best] += 1
                moved = true
            end
        end
        moved || break
        T .+= moves
    end
    topo[:, cols] .= T .* p.h
    return u
end

# ---------------------------------------------------------------------------------------------------
# High-water events

const _MEEB_PSTORM = [0.333, 0.333, 0.167, 0.310, 0.381, 0.310, 0.310, 0.310, 0.286, 0,
                      0.119,
                      0.024, 0.048, 0.048, 0.048, 0.071, 0.333, 0.286, 0.214, 0.190, 0.190,
                      0.262,
                      0.214, 0.262, 0.238]

"""
    barrier_storms() -> Matrix{Float64}
    barrier_storm_series() -> Matrix{Float64}

The high-water events of MEEB for North Core Banks, NC [Reeves2025](@cite) (see
`models/barrier_island/data/NOTICE.md`). `barrier_storms` gives the 9999
synthetic events of its stochastic forecasts, from a copula model of the
1979-2020 hindcast (columns: significant wave height Hs (m), duration (h),
still water level SL (m NAVD88), peak period Tp (s)). `barrier_storm_series`
gives the 222 events of that hindcast (columns: the week, counted in fiftieths
of a year from 1 January 1979, the total water level Rhigh and the lowest
runup Rlow (m NAVD88), and the duration (h)).
"""
barrier_storms() = _meeb_csv("synthetic_storms.csv")

@doc (@doc barrier_storms) barrier_storm_series() = _meeb_csv("storm_series.csv")

"""
    barrier_florence_site() -> NamedTuple

The test site of the paper's Fig. 5a, Table 1 [Reeves2025](@cite), from MEEB's
input files (see `models/barrier_island/data/NOTICE.md`): 500 m of North Core
Banks from 19250 m alongshore, the 701 m cross-shore from column 150 of MEEB's
domain, on 1 m cells (rows alongshore, the sea at the first column). Returns
the elevation (m NAVD88) before Hurricane Florence (2017 lidar) and after it
(October 2018) as `before` and `after`, the initial cover of the burial-tolerant
and burial-intolerant types (`spec1`, `spec2`), the overwash extent digitised
from imagery (`overwash`), and the event: `MHW` 0.39 m, a total water level
`Rhigh` of 3.32 m for `duration` 83 hours. Elevations are rounded to the
millimetre and covers to 1/250.
"""
function barrier_florence_site()
    path = joinpath(MODELS_DIR, "barrier_island", "data", "florence_site.bin")
    open(path) do io
        String(read(io, 8)) == "MEEBSITE" || throw(ArgumentError("not a MEEB site: $path"))
        m, n = Int(ltoh(read(io, Int32))), Int(ltoh(read(io, Int32)))
        layer(T) = ltoh.(read!(io, Matrix{T}(undef, m, n)))
        before = layer(Int16) ./ 1000
        after = layer(Int16) ./ 1000
        spec1 = layer(UInt8) ./ 250
        spec2 = layer(UInt8) ./ 250
        overwash = layer(UInt8) .!= 0
        return (; before, after, spec1, spec2, overwash, MHW=0.39, Rhigh=3.32, duration=83)
    end
end

function _meeb_csv(name)
    path = joinpath(MODELS_DIR, "barrier_island", "data", name)
    rows = [parse.(Float64, split(l, ',')) for l in Iterators.drop(eachline(path), 1)]
    return permutedims(reduce(hcat, rows))
end

# The beach slope of each row, as MEEB's code computes it: the elevation of the dune toe (the cell
# after the most landward cell seaward of the crest at most `toe_height` above MHW) over its distance
# from the shoreline; a row where that is not finite takes the previous row's (0.03 in the first).
function _meeb_beach_slopes(topo, xs, crest, MHW, toe_height)
    m = size(topo, 1)
    slopes = zeros(m)
    for r in 1:m
        seg = view(topo, r, 1:(crest[r] - 1))
        k = findfirst(<=(MHW + toe_height), reverse(seg))
        toe = crest[r] - (k === nothing ? 0 : k - 1)
        s = topo[r, toe] / (toe - xs[r])
        slopes[r] = isfinite(s) ? s : (r > 1 ? slopes[r - 1] : 0.03)
    end
    return slopes
end

# The total water level of event `e` (Hs, duration, SL, Tp) on each row's beach slope: the still water
# level plus the 2% exceedance runup of Stockdon et al. (2006).
function _meeb_twl(e, slopes)
    Hs, SL, Tp = e[1], e[3], e[4]
    L0 = 9.8 * Tp^2 / 2π
    return [SL +
            1.1 * (0.35β * sqrt(Hs * L0) +
                   sqrt((0.75β * sqrt(Hs * L0))^2 + (0.06 * sqrt(Hs * L0))^2) / 2)
            for β in slopes]
end

# A centred gradient with one-sided differences at the ends (numpy's), on 1 m cells.
function _meeb_gradient(y)
    n = length(y)
    g = zeros(n)
    n >= 2 || return g
    g[1] = y[2] - y[1]
    g[n] = y[n] - y[n - 1]
    for i in 2:(n - 1)
        g[i] = (y[i + 1] - y[i - 1]) / 2
    end
    return g
end

# Beach and dune-face change seaward of the crest over one overwash substep (Eq. 2-3): from the
# shoreline landward to the first cell below the water level whose landward neighbour is above it, or
# to the crest, the flux (βeq - β)(Ztwl - Z)², β the local slope, and the elevation falls by its
# divergence times Te / beach_substeps, beach_substeps times. Returns the change, the volume change of
# each row (m³) and the cells wetted.
function _meeb_beach!(E, crest, xs, MHW, Rhigh, p)
    m, n = size(E)
    pre = copy(E)
    wet = falses(m, n)
    Q = p.Te / p.beach_substeps
    for _ in 1:(p.beach_substeps), y in 1:m
        x0, xD, Rh = max(floor(Int, xs[y]), 1), crest[y], Rhigh[y]
        xD >= x0 || continue
        flux = zeros(xD - x0 + 1)
        wet[y, 1:(x0 - 1)] .= true
        for x in x0:xD
            zi = E[y, x]
            hprev = x == x0 ? MHW : E[y, x - 1]
            znext = E[y, min(x + 1, n)]
            stop = zi <= Rh < znext
            hnext = stop ? Rh : znext
            flux[x - x0 + 1] = (p.beta_eq - 0.5 * (hnext - hprev)) * (Rh - zi)^2
            stop && break
            wet[y, x] = true
        end
        E[y, x0:xD] .-= _meeb_gradient(flux) .* Q
    end
    Δ = E .- pre
    return Δ, vec(sum(Δ; dims=2)), wet
end

# Share discharge Q0 among the three landward neighbours by slope (Eq. 6-7, Barrier3D); S are the
# slopes down to the neighbours (NaN off the domain). Returns the three discharges.
function _meeb_share(Q0, S1, S2, S3, p, subaerial)
    nn = p.nn
    nz(q) = isnan(q) ? 0.0 : q
    if S1 > 0 || S2 > 0 || S3 > 0
        e1, e2, e3 = (S < 0 || isnan(S) ? 0.0 : S for S in (S1, S2, S3))
        t = e1^nn + e2^nn + e3^nn
        return nz(Q0 * e1^nn / t), nz(Q0 * e2^nn / t), nz(Q0 * e3^nn / t)
    elseif S1 == 0 || S2 == 0 || S3 == 0
        Qx = Q0 / ((S1 == 0) + (S2 == 0) + (S3 == 0))
        return (S1 == 0 ? Qx : 0.0), (S2 == 0 ? Qx : 0.0), (S3 == 0 ? Qx : 0.0)
    end
    if isnan(S1)
        a2, a3 = abs(S2)^(-nn), abs(S3)^(-nn)
        Q1, Q2, Q3 = 0.0, Q0 * a2 / (a2 + a3), Q0 * a3 / (a2 + a3)
    elseif isnan(S3)
        a1, a2 = abs(S1)^(-nn), abs(S2)^(-nn)
        Q1, Q2, Q3 = Q0 * a1 / (a1 + a2), Q0 * a2 / (a1 + a2), 0.0
    else
        a1, a2, a3 = abs(S1)^(-nn), abs(S2)^(-nn), abs(S3)^(-nn)
        t = a1 + a2 + a3
        Q1, Q2, Q3 = Q0 * a1 / t, Q0 * a2 / t, Q0 * a3 / t
    end
    Q1, Q2, Q3 = nz(Q1), nz(Q2), nz(Q3)
    if subaerial
        up(Q, S) = abs(S) > p.max_upslope ? 0.0 : Q * (1 - abs(S) / p.max_upslope)
        Q1, Q2, Q3 = up(Q1, S1), up(Q2, S2), up(Q3, S3)
    end
    return Q1, Q2, Q3
end

# One substep of overwash landward of the crest (Eq. 4-9, Barrier3D): water enters at each
# overtopped crest cell (Q = U (Ztwl - ZD) per metre, U = √(2 g (Ztwl - ZD)), in m³/h), loses Rin at
# each cross-shore step, is shared among the three landward neighbours by slope and reduced by their
# plants (Eq. 23, by the burial-tolerant type's coefficient where it grows), and carries sand
# Kow (Q (S + Cs))^m; in the bay the sand arriving is passed on, Cbb of it, in proportion to the
# discharge. Returns the cells with discharge and the sand in minus the sand out (m³/h).
function _meeb_overwash(E, crest, MHW, Rhigh, s1, s2, p)
    m, n = size(E)
    Dis = zeros(m, n)
    In = zeros(m, n)
    Out = zeros(m, n)
    zD = [E[r, crest[r]] for r in 1:m]
    any(Rhigh .> zD) || return Dis .> 0, In .- Out
    for r in 1:m
        R = max(Rhigh[r] - zD[r], 0.0)
        Dis[r, crest[r]] = sqrt(2 * 9.8 * R) * R * 3600
    end
    qmin = p.overwash_min_discharge / p.overwash_substeps
    qbbmin = p.overwash_min_bay_discharge / p.overwash_substeps
    slow(Q, a, b) = a > 0 ? Q * (1 - p.flow_reduction1 * a) :
                    Q * (1 - p.flow_reduction2 * b)
    # A negative stream power carries no sand (NaN, then zero, in MEEB's compiled code).
    sed(Q, S) = Q > qmin && Q * (S + p.Cs) >= 0 ? p.Kow * (Q * (S + p.Cs))^p.mm : 0.0
    for d in minimum(crest):(n - 1)
        d > maximum(crest) + 2 && sum(view(Dis, :, d)) <= 0 && break
        for i in 1:m
            d > 1 && Dis[i, d] > 0 && (Dis[i, d] -= p.Rin)
            Dis[i, d] < 0 && (Dis[i, d] = 0.0)
        end
        for i in 1:m
            Dis[i, d] > 0 || continue
            S1 = i > 1 ? (E[i, d] - E[i - 1, d + 1]) / sqrt(2) : NaN
            S2 = E[i, d] - E[i, d + 1]
            S3 = i < m ? (E[i, d] - E[i + 1, d + 1]) / sqrt(2) : NaN
            Q1, Q2, Q3 = _meeb_share(Dis[i, d], S1, S2, S3, p, E[i, d] > MHW)
            if i > 1
                Q1 = slow(Q1, s1[i - 1, d], s2[i - 1, d])
                Dis[i - 1, d + 1] += Q1
            end
            Q2 = slow(Q2, s1[i, d], s2[i, d])
            Dis[i, d + 1] += Q2
            if i < m
                Q3 = slow(Q3, s1[i + 1, d], s2[i + 1, d])
                Dis[i + 1, d + 1] += Q3
            end
            if E[i, d] > MHW || any(>(MHW), view(E, i, (d + 1):min(d + 9, n)))
                Qs1, Qs2, Qs3 = sed(Q1, S1), sed(Q2, S2), sed(Q3, S3)
                Qs1, Qs2, Qs3 = (isnan(q) ? 0.0 : q for q in (Qs1, Qs2, Qs3))
            else
                Qs0 = In[i, d] * p.Cbb
                t = Q1 + Q2 + Q3
                Qs1, Qs2, Qs3 = t != 0 ? (Qs0 * Q1 / t, Qs0 * Q2 / t, Qs0 * Q3 / t) :
                                (0.0, 0.0, 0.0)
                Qs1, Qs2, Qs3 = (isnan(q) || q < qbbmin ? 0.0 : q for q in (Qs1, Qs2, Qs3))
            end
            i > 1 && (In[i - 1, d + 1] += Qs1)
            In[i, d + 1] += Qs2
            i < m && (In[i + 1, d + 1] += Qs3)
            Out[i, d] = Qs1 + Qs2 + Qs3
        end
    end
    return Dis .> 0, In .- Out
end

# A high-water event of `dur` hours at total water level `Rhigh` (each row), in hourly steps of
# `overwash_substeps` substeps: the crest is found each hour, overwash moves sand landward of it
# (limited to `flux_limit` m an hour; the crest cells themselves change only from the seaward side)
# and the beach and dune face change seaward of it. Returns the elevation change, the overwash volume
# landward of the crest and the beach-dune volume change of each row (m³), and the cells flooded.
function _meeb_event!(u, p, MHW, Rhigh, dur, xs)
    E = copy(u.topo)
    m, n = size(E)
    ow = zeros(m)
    bd = zeros(m)
    flooded = falses(m, n)
    sub = p.overwash_substeps
    lim = p.flux_limit / sub
    crest = barrier_crest(E, MHW)[1]
    for TS in 0:(floor(Int, dur) * sub - 1)
        TS % sub == 0 && (crest = barrier_crest(E, MHW)[1])
        wetL, net = _meeb_overwash(E, crest, MHW, Rhigh, u.spec1, u.spec2, p)
        ΔL = clamp.(net ./ sub, -lim, lim)
        for r in 1:m
            ΔL[r, crest[r]] = 0.0
        end
        ow .+= vec(sum(ΔL; dims=2))
        flooded .|= wetL
        ΔS, dV, wetS = _meeb_beach!(copy(E), crest, xs, MHW, Rhigh, p)
        ΔS ./= sub
        bd .+= dV ./ sub
        flooded .|= wetS
        E .+= ΔS .+ ΔL
    end
    Δ = E .- u.topo
    u.topo .= E
    return Δ, ow, bd, flooded
end

"""
    barrier_event!(u, p, Rhigh, duration; MHW=p.MHW, shoreline=nothing, rng=Random.Xoshiro(0))
        -> NamedTuple

One high-water event of MEEB [Reeves2025](@cite) on the state `u` (with at
least `topo`, `spec1`, `spec2` and `veg`), in place: for `duration` hours at the
total water level `Rhigh` (m, one value or one for each row), in hourly steps
of `p.overwash_substeps`, the beach and dune face change seaward of the
foredune crest (Eq. 2-3) and overwash is routed landward of it (Eq. 4-9), as in
[`barrier_step!`](@ref); then slopes steeper than the angle of repose
avalanche, as in MEEB's hindcast of a single storm. `shoreline` gives the
cross-shore column of each row's ocean shoreline (by default the first cell at
or above `MHW`). Plants are left as they are. Returns the elevation `change`,
the overwash volume landward of the crest and the beach-dune volume change of
each row (`overwash`, `beach`, m³) and the cells the water reached (`flooded`).
"""
function barrier_event!(u, p, Rhigh, duration::Real; MHW::Real=p.MHW, shoreline=nothing,
                        rng=Random.Xoshiro(0))
    m = size(u.topo, 1)
    xs = shoreline === nothing ?
         Float64[_meeb_ocean_shoreline(u.topo, r, MHW) for r in 1:m] : Float64.(shoreline)
    R = Rhigh isa Real ? fill(Float64(Rhigh), m) : Float64.(Rhigh)
    pre = copy(u.topo)
    _, ow, bd, flooded = _meeb_event!(u, p, MHW, R, duration, xs)
    _meeb_avalanche!(u, p, rng, MHW)
    return (; change=u.topo .- pre, overwash=ow, beach=bd, flooded)
end

# ---------------------------------------------------------------------------------------------------
# The shoreline

"""
    barrier_shoreface(; Hs=0.98, T=6.6, D50=2e-4, R=1.65, es=0.01, cs=0.01) -> NamedTuple

The shoreface of MEEB, estimated from the wave climate and the sediment as in
Nienhuis and Lorenzo-Trueba (2019) (Eq. 15-18 of [Reeves2025](@cite)): the depth
`D` (m, Hallermeier), the equilibrium slope `S_eq`, the response rate `k`
(m³/m/yr) and the settling velocity `ws` (m/s, Ferguson and Church), with the
breaking depth `2 Hs / 0.78`, as in MEEB's code.
"""
function barrier_shoreface(; Hs::Real=0.98, T::Real=6.6, D50::Real=2e-4, R::Real=1.65,
                           es::Real=0.01, cs::Real=0.01)
    g = 9.81
    ws = R * g * D50^2 / (18e-6 + sqrt(0.75 * R * g * D50^3))
    z0 = 2Hs / 0.78
    D = 0.018 * Hs * T * sqrt(g / (R * D50))
    S_eq = 3ws / 4 / sqrt(D * g) * (5 + 3T^2 * g / 4 / π^2 / D)
    k = 3600 * 24 * 365 * (es * cs * g^(11 / 4) * Hs^5 * T^(5 / 2)) /
        (960 * R * π^(7 / 2) * ws^2) *
        ((1 / (11 / 4 * z0^(11 / 4)) - 1 / (11 / 4 * D^(11 / 4))) / (D - z0))
    return (; D, S_eq, k, ws)
end

# The alongshore transport constant of Nienhuis et al. (2015), as MEEB 1.0's
# calc_alongshore_transport_k (about 0.063 m^(3/5) s^(-6/5)).
function _meeb_transport_k()
    g = 9.80665
    return 5.3e-6 * 1050 * g^1.5 * 0.5^1.2 * (sqrt(g * 0.78) / 2π)^0.2
end

"""
    barrier_diffusivity(; a=0.6, h=0.39, Hs=0.98, T=6.6, D=barrier_shoreface().D) -> Vector{Float64}

The wave-climate-averaged shoreline diffusivity (m²/yr) of MEEB at the 181
wave angles from -90° to 90° in steps of 1° (Eq. 20-22 of [Reeves2025](@cite),
from Ashton and Murray 2006 and Nienhuis and Lorenzo-Trueba 2019): the
diffusivity at each wave angle convolved with the distribution of wave energy,
a fraction `a` of waves from the left and `h` at high angles (beyond 45°), as
in MEEB's code.
"""
function barrier_diffusivity(; a::Real=0.6, h::Real=0.39, Hs::Real=0.98, T::Real=6.6,
                             D::Real=barrier_shoreface().D)
    nb = 181
    start, stop = -π / 2, π / 2
    step = (stop - start) / (nb - 1)
    θ = [k == nb - 1 ? stop : k * step + start for k in 0:(nb - 1)]
    xb = deg2rad.([-90.0, -45.0, 0.0, 45.0, 90.0])
    f = [0.0, a * h, a * (1 - h), (1 - a) * (1 - h), (1 - a) * h] .* 4 ./ π
    shifted = nextfloat.(xb)
    pdf = [t < xb[1] || t > xb[end] ? 0.0 : f[min(count(<=(t), shifted) + 1, 5)] * step
           for t in θ]
    k = _meeb_transport_k()
    diff = [-(k / D * Hs^2.4 * T^0.2) * 365 * 24 * 3600 * cos(t)^0.2 *
            (1.2 * sin(t)^2 - cos(t)^2) for t in θ]
    full = zeros(2nb - 1)
    for i in 1:nb, j in 1:nb
        full[i + j - 1] += pdf[i] * diff[j]
    end
    first = (nb - 1) - (nb - 1) ÷ 2
    return full[(first + 1):(first + nb)]
end

# Cross-shore shoreline change (Eq. 13-14, LTA14) over one shoreline step of 1/25 year: the shoreface
# flux k (S_eq - S), S = D / (xs - xt), moves the toe by 4 qsf / D + 2 RSLR / S and the shoreline by
# 2 (qow + qbd) / D - 4 qsf / D, each averaged over sections of `section` rows. Returns S.
function _meeb_cross_shore!(xs, xt, ow, bd, p)
    m = length(xs)
    rslr = p.RSLR / 25
    k = p.k_sf / 25
    S = p.D_sf ./ (xs .- xt)
    qsf = k .* (p.S_sf_eq .- S)
    dt = 4 .* qsf ./ p.D_sf .+ 2 .* rslr ./ S
    ds = 2 .* (max.(ow, 0.0) .+ bd) ./ p.D_sf .- 4 .* qsf ./ p.D_sf
    for lo in 1:(p.section):m
        rows = lo:min(lo + p.section - 1, m)
        xt[rows] .+= sum(view(dt, rows)) / length(rows)
        xs[rows] .+= sum(view(ds, rows)) / length(rows)
    end
    return S
end

# Solve the tridiagonal system with sub-, main and super-diagonals a, b, c (Thomas).
function _meeb_tridiagonal(a, b, c, d)
    n = length(b)
    cp = zeros(n)
    dp = zeros(n)
    cp[1] = n > 1 ? c[1] / b[1] : 0.0
    dp[1] = d[1] / b[1]
    for i in 2:n
        den = b[i] - a[i - 1] * cp[i - 1]
        i < n && (cp[i] = c[i] / den)
        dp[i] = (d[i] - a[i - 1] * dp[i - 1]) / den
    end
    x = zeros(n)
    x[n] = dp[n]
    for i in (n - 1):-1:1
        x[i] = dp[i] - cp[i] * x[i + 1]
    end
    return x
end

# Alongshore shoreline diffusion (Eq. 19) over `dt` years by Crank-Nicolson on the shoreline of the
# first row of each section, the end sections held fixed.
function _meeb_alongshore!(xs, diffusivity, p, dt)
    m = length(xs)
    y = xs[1:(p.section):m]
    ny = length(y)
    ny >= 3 || return xs
    nxt = [y[2:end]; y[1]]
    prv = [y[end]; y[1:(end - 1)]]
    angle = 180 .* atan.(nxt .- y, p.section) ./ π
    r = [max(0.0,
             diffusivity[clamp(round(Int, 90 - a), 1, 180) + 1] * dt / 2 / p.section^2)
         for a in angle]
    r[1] = r[end] = 0.0
    rhs = y .+ r .* (nxt .- 2 .* y .+ prv)
    new = _meeb_tridiagonal(-r[2:end], 1 .+ 2 .* r, -r[1:(end - 1)], rhs)
    for (j, lo) in enumerate(1:(p.section):m)
        xs[lo:min(lo + p.section - 1, m)] .= new[j]
    end
    return xs
end

# Regrade the beach and shoreface where the shoreline has moved a whole cell: on erosion the old
# shoreline cell is lowered below MHW and a linear shoreface of slope S drawn seaward of it; on
# accretion the new beach cells take the mean height of the five most seaward old beach cells.
function _meeb_regrade!(topo, xs, prev, MHW, S, p)
    m, n = size(topo)
    rslr = p.RSLR / 25
    for r in 1:m
        target, old = floor(Int, xs[r]), floor(Int, prev[r])
        if target > old
            if target <= n
                topo[r, target] = min(MHW - (MHW - topo[r, old]) / 2, MHW - rslr)
                for c in 1:(target - 1)
                    topo[r, c] = (c - target) * S[r] + topo[r, target]
                end
            else
                for c in 1:n
                    topo[r, c] = (c - 1 - n) * S[r] + MHW
                end
            end
        elseif target < old && 1 < target <= n + 1
            beach = old:min(old + 4, n)
            topo[r, target:(old - 1)] .= sum(view(topo, r, beach)) / length(beach) + rslr
            for c in 1:(target - 1)
                topo[r, c] = (c - target) * S[r] + MHW
            end
        end
    end
    return topo
end

# ---------------------------------------------------------------------------------------------------
# Vegetation

# Piecewise-linear growth of a plant type with vertices `v` (m/yr of burial) and peak growth `peak`:
# -1 at and beyond the outer vertices, 0 at the second and fourth, `peak` at the third.
function _meeb_growth(Δ, v, peak)
    x1, x2, x3, x4, x5 = v
    Δ < x1 && return -1.0
    Δ < x2 && return -1 + (Δ - x1) / (x2 - x1)
    Δ < x3 && return peak * (Δ - x2) / (x3 - x2)
    Δ < x4 && return peak - peak * (Δ - x3) / (x4 - x3)
    Δ < x5 && return -(Δ - x4) / (x5 - x4)
    return -1.0
end

# A year of plant growth from the burial in `balance`: each type keeps what it loses and gains only
# where present, spread from a neighbour of its own type or newly established, at most the cover
# left free.
function _meeb_vegetation!(u, p, rng, MHW)
    s1, s2 = copy(u.spec1), copy(u.spec2)
    near(A, I) = any(J -> checkbounds(Bool, A, J) && A[J] > 0,
                     (I + CartesianIndex(di, dj) for di in -1:1, dj in -1:1))
    for I in CartesianIndices(u.topo)
        z = u.topo[I]
        g1 = clamp(s1[I] + _meeb_growth(u.balance[I], p.vertices1, p.peak1), 0.0, 1.0) -
             s1[I]
        g2 = clamp(s2[I] + _meeb_growth(u.balance[I], p.vertices2, p.peak2), 0.0, 1.0) -
             s2[I]
        dry = z > MHW
        settles = 0 <= z - MHW < 20
        lat1 = dry && near(s1, I) && rand(rng) < p.p_lateral1
        lat2 = dry && near(s2, I) && rand(rng) < p.p_lateral2
        new1 = dry && settles && s1[I] <= 0 && rand(rng) < p.p_pioneer1
        new2 = dry && settles && s2[I] <= 0 && u.change[I] == 0 && rand(rng) < p.p_pioneer2
        free = 1 - u.veg[I]
        a1 = (s1[I] > 0 || lat1 || new1) ? min(free, max(g1, 0.0)) : 0.0
        a2 = (s2[I] > 0 || lat2 || new2) ? min(free, max(g2, 0.0)) : 0.0
        u.spec1[I] = z <= p.elev_min1 + MHW ? 0.0 : s1[I] + a1 + min(g1, 0.0)
        u.spec2[I] = z <= p.elev_min2 + MHW ? 0.0 : s2[I] + a2 + min(g2, 0.0)
    end
    u.veg .= clamp.(clamp.(u.spec1, 0.0, 1.0) .+ clamp.(u.spec2, 0.0, 1.0), 0.0, p.maxveg)
    u.veg_eff .= _meeb_gaussian(u.veg, p.veg_sigma; mode=:constant)
    return u
end

# ---------------------------------------------------------------------------------------------------
# The model

"""
    build_barrier_island() -> CAModel

MEEB v1.0, the Mesoscale Explicit Ecogeomorphic Barrier model of
[Reeves2025](@cite). Each cell of 1 m holds its elevation (`topo`, m against
the datum of `MHW`, NAVD88 at North Core Banks), the cover of a burial-tolerant
(`spec1`, dune grasses) and a burial-intolerant (`spec2`, woody shrubs) plant
type, their combined cover (`veg`, at most `maxveg`) and that cover smoothed by
a Gaussian of `veg_sigma` m (`veg_eff`, updated once a year, which the wind
sees), the elevation gained over the current balance period (`balance`, m) and
its total absolute change (`change`, m), the cross-shore positions of the ocean
shoreline and the shoreface toe of its row (`shoreline`, `toe`, in columns, the
same along each row), and the number of high-water events that flooded it
(`inundated`). Rows run alongshore and columns cross-shore, from the sea at the
first column to the back-barrier bay. The dynamics are
[`barrier_step!`](@ref), a callback stored with the problem.
"""
function build_barrier_island()
    rule = LocalRule(Offsets((1, 0))) do c, nbrs, ctx
        return c
    end
    states = [StateVariable(:topo, 0.0; description="elevation (m)"),
              StateVariable(:spec1, 0.0; description="cover of the burial-tolerant type"),
              StateVariable(:spec2, 0.0; description="cover of the burial-intolerant type"),
              StateVariable(:veg, 0.0; description="combined plant cover"),
              StateVariable(:veg_eff, 0.0;
                            description="combined cover smoothed by a Gaussian, as the wind sees it"),
              StateVariable(:balance, 0.0;
                            description="elevation gained over the balance period (m)"),
              StateVariable(:change, 0.0;
                            description="absolute elevation change over the balance period (m)"),
              StateVariable(:shoreline, 0.0;
                            description="cross-shore position of the row's ocean shoreline (columns)"),
              StateVariable(:toe, 0.0;
                            description="cross-shore position of the row's shoreface toe (columns)"),
              StateVariable(:inundated, 0;
                            description="high-water events that flooded the cell")]
    sf = barrier_shoreface()
    parameters = [Parameter(:MHW, 0.39; description="initial mean high water (m)"),
                  Parameter(:RSLR, 0.0; description="relative sea-level rise (m/yr)"),
                  Parameter(:h, 0.02; description="aeolian slab height (m)"),
                  Parameter(:L, 5; description="mean saltation length (cells)"),
                  Parameter(:L_dev, 2;
                            description="half-range of the uniform saltation length (cells)"),
                  Parameter(:groundwater_depth, 0.4;
                            description="fraction of the smoothed height above MHW at the water table"),
                  Parameter(:wind_rose, (0.81, 0.04, 0.06, 0.09);
                            description="probabilities of wind towards the land, alongshore down, the sea, alongshore up"),
                  Parameter(:pe0, 0.10; description="erosion probability without plants"),
                  Parameter(:pd0, 0.22;
                            description="deposition probability without plants"),
                  Parameter(:pd1, 0.54; description="deposition probability at full cover"),
                  Parameter(:rho_entrain, 0.10;
                            description="cover at which entrainment stops (ρq0)"),
                  Parameter(:rho_saltation, 0.35;
                            description="cover that traps slabs part-way through a hop (ρv)"),
                  Parameter(:shadow, 12.0; description="shadow angle (degrees)"),
                  Parameter(:repose_bare, 20.0; description="angle of repose of bare sand"),
                  Parameter(:repose_veg, 30.0;
                            description="angle of repose of vegetated sand"),
                  Parameter(:repose_threshold, 0.3;
                            description="cover from which the vegetated angle applies"),
                  Parameter(:D_bb, 1.5;
                            description="equilibrium depth of the back-barrier bay (m)"),
                  Parameter(:beta_eq, 0.022; description="equilibrium beach slope"),
                  Parameter(:Te, 1.48;
                            description="erosive timescale coefficient seaward of the crest"),
                  Parameter(:beach_substeps, 25;
                            description="substeps of the beach equation in each overwash substep"),
                  Parameter(:Rin, 249.0;
                            description="infiltration and drag of overwash flow (m³/h)"),
                  Parameter(:max_upslope, 1.5; description="steepest slope water flows up"),
                  Parameter(:Kow, 0.0001684;
                            description="overwash sediment transport coefficient"),
                  Parameter(:Cs, 0.0283; description="overwash flow momentum constant"),
                  Parameter(:mm, 1.04;
                            description="exponent of the overwash sediment flux"),
                  Parameter(:nn, 0.5; description="flow routing exponent"),
                  Parameter(:Cbb, 0.7;
                            description="fraction of the sand passed on in the back-barrier bay"),
                  Parameter(:flux_limit, 1.0;
                            description="largest overwash elevation change an hour (m)"),
                  Parameter(:overwash_min_discharge, 1.0;
                            description="least discharge that carries sand (m³/h)"),
                  Parameter(:overwash_min_bay_discharge, 1.0;
                            description="least sand flux passed on in the bay (m³/h)"),
                  Parameter(:overwash_substeps, 25;
                            description="substeps of each storm hour"),
                  Parameter(:toe_height, 1.67;
                            description="mean dune toe height above MHW (m), for beach slopes"),
                  Parameter(:pstorm, copy(_MEEB_PSTORM);
                            description="probability of an event in each fortnight of the year"),
                  Parameter(:storms, Matrix{Float64}(undef, 0, 4);
                            description="synthetic events (Hs, duration, SL, Tp), drawn at random"),
                  Parameter(:storm_series, Matrix{Float64}(undef, 0, 4);
                            description="hindcast events (week, Rhigh, Rlow, duration), used in place of the synthetic ones if not empty"),
                  Parameter(:start_week, 0;
                            description="week of the hindcast series at which the run starts"),
                  Parameter(:D_sf, sf.D; description="shoreface depth (m)"),
                  Parameter(:S_sf_eq, sf.S_eq; description="equilibrium shoreface slope"),
                  Parameter(:k_sf, sf.k; description="shoreface flux rate (m³/m/yr)"),
                  Parameter(:section, 25;
                            description="alongshore length of the shoreline sections (rows)"),
                  Parameter(:diffusivity, barrier_diffusivity();
                            description="shoreline diffusivity by wave angle (m²/yr)"),
                  Parameter(:vertices1, (-1.5, -0.05, 0.5, 1.5, 2.2);
                            description="burial rates (m/yr) shaping burial-tolerant growth"),
                  Parameter(:peak1, 0.2; description="peak yearly growth, burial-tolerant"),
                  Parameter(:vertices2, (-1.6, -0.7, 0.0, 0.2, 2.1);
                            description="burial rates (m/yr) shaping burial-intolerant growth"),
                  Parameter(:peak2, 0.05;
                            description="peak yearly growth, burial-intolerant"),
                  Parameter(:p_lateral1, 0.2;
                            description="yearly probability of lateral expansion, burial-tolerant"),
                  Parameter(:p_lateral2, 0.2;
                            description="yearly probability of lateral expansion, burial-intolerant"),
                  Parameter(:p_pioneer1, 0.05;
                            description="yearly probability of establishment, burial-tolerant"),
                  Parameter(:p_pioneer2, 0.05;
                            description="yearly probability of establishment, burial-intolerant"),
                  Parameter(:elev_min1, 0.25;
                            description="least elevation above MHW, burial-tolerant (m)"),
                  Parameter(:elev_min2, 0.25;
                            description="least elevation above MHW, burial-intolerant (m)"),
                  Parameter(:maxveg, 1.0; description="largest combined cover"),
                  Parameter(:flow_reduction1, 0.02;
                            description="overwash flow reduction at full cover, burial-tolerant"),
                  Parameter(:flow_reduction2, 0.05;
                            description="overwash flow reduction at full cover, burial-intolerant"),
                  Parameter(:veg_sigma, 3.0;
                            description="standard deviation of the smoothing of cover (m)"),
                  Parameter(:balance_period, :year;
                            description="period whose burial the plants respond to: :year (the paper) or :week (MEEB's code)")]
    return CAModel(rule; name=:barrier_island, states, parameters,
                   description="Reeves et al. (2025) MEEB barrier island evolution")
end

"""
    barrier_step!(u, p, t, rng, storm_rng) -> NamedTuple

Week `t` of MEEB [Reeves2025](@cite) (50 a year), applied to the state `u` in
place, with `rng` for the wind, avalanches and plants and `storm_rng` for the
high-water events. MHW stands at `MHW + RSLR t / 50`. The wind blows towards a
direction drawn from `wind_rose`. Between the most seaward ocean and most
landward back-barrier shorelines, each cell outside the shadow (`shadow`
degrees), above the water table (MHW plus `groundwater_depth` times the height
above MHW of the topography smoothed by a Gaussian of 12 m, at least MHW) and
more than a slab above MHW gives up a slab of `h` m with probability
`pe0 (1 - ρ / rho_entrain)` (at least 0), ρ being the smoothed cover `veg_eff`.
Each slab hops a length drawn each week uniformly from `L ± L_dev` cells,
landing at the end of each hop with probability `pd0 + (pd1 - pd0) ρ` (1 in the
shadow or at or below the water table) and, part-way, with that probability
where `ρ ≥ rho_saltation`; slabs leaving the domain are lost.

Every second week (`t` odd) a high-water event occurs with the probability of
its fortnight (`pstorm`), drawn from `storms` (its runup computed on each row's
beach, Stockdon et al. 2006, then raised by the sea-level rise so far) or taken
from `storm_series`. For its duration it changes the beach and dune face seaward
of the foredune crest ([`barrier_crest`](@ref)) and routes overwash landward of
it, hour by hour; plants die where it floods (the flooded cells smoothed by a
Gaussian of 2 m and rounded). The ocean shoreline then moves with the shoreface
and alongshore diffusion, the beach and shoreface are regraded where it moved a
whole cell, slopes steeper than the angle of repose avalanche, and the bay is
kept at least `D_bb` deep.

In the first week of each year after the first, the plants grow with their
burial over the year (`balance_period = :year`) or over that week only
(`balance_period = :week`, as MEEB's code does): each type only where present,
spread from a neighbour of its own type (`p_lateral1`, `p_lateral2`) or newly
established (`p_pioneer1`, `p_pioneer2`, 0 to 20 m above MHW; for `spec2` only
where the elevation did not change), its gain at most the cover left free; they
die below `elev_min1`, `elev_min2` above MHW. Returns the wind direction,
whether an event occurred, and its highest total water level (NaN without one).
"""
function barrier_step!(u, p, t::Integer, rng, storm_rng)
    it = t - 1
    MHW = p.MHW + p.RSLR * t / 50
    pre = copy(u.topo)
    rose = cumsum(collect(p.wind_rose))
    draw = rand(rng) * rose[end]
    dir = something(findfirst(>=(draw), rose), 4)
    u.topo .+= _meeb_aeolian(u, p, rng, MHW, dir) .* p.h
    u.balance .+= u.topo .- pre
    u.change .+= abs.(u.topo .- pre)
    event, twl = false, NaN
    if it % 2 == 0
        m = size(u.topo, 1)
        xs = u.shoreline[:, 1]
        xt = u.toe[:, 1]
        Rhigh, dur = Float64[], 0.0
        if !isempty(p.storm_series)
            k = findfirst(==(it + p.start_week), view(p.storm_series, :, 1))
            if k !== nothing && p.storm_series[k, 2] > MHW
                event = true
                Rhigh, dur = fill(p.storm_series[k, 2], m), p.storm_series[k, 4]
            end
        elseif !isempty(p.storms) && rand(storm_rng) < p.pstorm[(it % 50) ÷ 2 + 1]
            e = view(p.storms, rand(storm_rng, axes(p.storms, 1)), :)
            crest = barrier_crest(u.topo, MHW)[1]
            Rhigh = _meeb_twl(e, _meeb_beach_slopes(u.topo, xs, crest, MHW, p.toe_height))
            dur = e[2]
            event = any(>(MHW), Rhigh)
            Rhigh .+= MHW - p.MHW
        end
        ow, bd = zeros(m), zeros(m)
        if event
            twl = maximum(Rhigh)
            Δ, ow, bd, flooded = _meeb_event!(u, p, MHW, Rhigh, dur, xs)
            dead = round.(_meeb_gaussian(Float64.(flooded), 2.0)) .> 0
            u.spec1[dead] .= 0.0
            u.spec2[dead] .= 0.0
            u.veg .= u.spec1 .+ u.spec2
            u.inundated .+= dead
            u.balance .+= Δ
            u.change .+= abs.(Δ)
        end
        prev = copy(xs)
        S = _meeb_cross_shore!(xs, xt, ow, bd, p)
        _meeb_alongshore!(xs, p.diffusivity, p, 2 / 50)
        before = copy(u.topo)
        _meeb_regrade!(u.topo, xs, prev, MHW, S, p)
        _meeb_avalanche!(u, p, rng, MHW)
        u.balance .+= u.topo .- before
        u.change .+= abs.(u.topo .- before)
        for r in 1:m
            u.shoreline[r, :] .= xs[r]
            u.toe[r, :] .= xt[r]
            for c in _meeb_backbarrier(u.topo, r, MHW):size(u.topo, 2)
                u.topo[r, c] < MHW - p.D_bb && (u.topo[r, c] = MHW - p.D_bb)
            end
        end
    end
    it % 50 == 0 && it > 0 && _meeb_vegetation!(u, p, rng, MHW)
    if p.balance_period === :week || it % 50 == 0
        u.balance .= 0
        u.change .= 0
    end
    return (; wind=dir, event, twl)
end

"""
    barrier_callback() -> DiscreteCallback

A callback that runs [`barrier_step!`](@ref) after each step, with random number
generators seeded from the problem's seed and the week: one for the wind,
avalanches and plants, one for the high-water events.
"""
function barrier_callback()
    affect! = function (integ)
        t = Int(integ.t)
        rng = Random.Xoshiro(hash((integ.prob.seed, t, :barrier)))
        storm_rng = Random.Xoshiro(hash((integ.prob.seed, t, :barrier_storm)))
        barrier_step!(integ.u, integ.p, t, rng, storm_rng)
        return nothing
    end
    return DiscreteCallback((u, t, integ) -> true, affect!; save_positions=(false, false))
end

"""
    barrier_profile(; MHW=0.39, offshore=60, beach=76, dune=4.0, dune_width=30,
                    interior=200, interior_height=1.2, bay=60, D_bb=1.5,
                    S_sf=barrier_shoreface().S_eq) -> Vector{Float64}

An idealised cross-shore profile of an undeveloped barrier (m against the datum
of `MHW`), from the sea to the back-barrier bay, standing in for the lidar of
North Core Banks that [Reeves2025](@cite) start from: a shoreface of slope
`S_sf` over `offshore` m below MHW (MEEB's initial elevations add one), a beach
rising linearly over `beach` m to the dune toe 1.67 m above MHW, a foredune
whose crest stands `dune` m above MHW (a sine bump `dune_width` m wide on the
line from the toe to the interior; none if `dune` is 0), a level interior
`interior` m wide at `interior_height` m, a back-barrier slope down to `D_bb` m
below MHW over 30 m, and a level bay `bay` m wide.
"""
function barrier_profile(; MHW::Real=0.39, offshore::Integer=60, beach::Integer=76,
                         dune::Real=4.0, dune_width::Integer=30, interior::Integer=200,
                         interior_height::Real=1.2, bay::Integer=60, D_bb::Real=1.5,
                         S_sf::Real=barrier_shoreface().S_eq)
    dune >= 0 || throw(ArgumentError("dune must be nonnegative, got $dune"))
    toe = MHW + 1.67
    z = [MHW - k * S_sf for k in offshore:-1:1]
    append!(z, [MHW + (toe - MHW) * k / beach for k in 0:(beach - 1)])
    rise = dune > 0 ? max(0.0, MHW + dune - (toe + interior_height) / 2) : 0.0
    for k in 0:(dune_width - 1)
        f = k / (dune_width - 1)
        push!(z, toe + (interior_height - toe) * f + rise * sinpi(f))
    end
    append!(z, fill(float(interior_height), interior))
    append!(z, [interior_height + (MHW - D_bb - interior_height) * k / 30 for k in 1:30])
    append!(z, fill(MHW - D_bb, bay))
    return z
end

"""
    problem_barrier_island(; scenario=:dune, alongshore=100, years=10, storms=:stochastic,
                           start_week=0, seed=0, pmap=(;), kwargs...) -> CAProblem

A barrier island for MEEB [Reeves2025](@cite): `alongshore` rows of 1 m, each the
profile of [`barrier_profile`](@ref), for `years` years of 50 weeks. With
`scenario = :dune` the foredune stands 4 m above MHW everywhere; with `:gap` it
is absent in the middle fifth of the rows (an overwash gap); with `:low` it
stands only 2 m above MHW. Burial-tolerant grasses cover the foredune and the
10 m of beach before it (cover 0.6 to 1, at random), and burial-intolerant
shrubs a random half of the interior behind it (0.4 to 0.8), as the paper's
classes of initial density; neither grows within 0.25 m of MHW. High-water
events are drawn from MEEB's synthetic events for North Core Banks
([`barrier_storms`](@ref)) with `storms = :stochastic`, or taken from its
1979-2020 hindcast ([`barrier_storm_series`](@ref)) from week `start_week` with
`storms = :hindcast`. The shoreline starts at each section's mean first cell at
or above MHW and the shoreface toe `D_sf / S_sf_eq` m seaward of it. `pmap`
overrides parameters.
"""
function problem_barrier_island(; scenario::Symbol=:dune, alongshore::Integer=100,
                                years::Integer=10, storms::Symbol=:stochastic,
                                start_week::Integer=0, seed::Integer=0, pmap=(;),
                                kwargs...)
    scenario in (:dune, :gap, :low) ||
        throw(ArgumentError("scenario must be :dune, :gap or :low, got :$scenario"))
    storms in (:stochastic, :hindcast) ||
        throw(ArgumentError("storms must be :stochastic or :hindcast, got :$storms"))
    alongshore >= 1 || throw(ArgumentError("alongshore must be positive, got $alongshore"))
    model = build_barrier_island()
    defaults = Dict(pr.name => pr.default for pr in model.parameters)
    none = Matrix{Float64}(undef, 0, 4)
    given = merge((; storms=storms === :stochastic ? barrier_storms() : none,
                   storm_series=storms === :hindcast ? barrier_storm_series() : none,
                   start_week),
                  (; (Symbol(k) => v for (k, v) in pairs(pmap))...))
    par(k) = haskey(given, k) ? given[k] : defaults[k]
    MHW = par(:MHW)
    rng = Random.Xoshiro(hash((seed, 0, :barrier)))
    gap = scenario === :gap ? ((alongshore * 2 ÷ 5 + 1):(alongshore * 3 ÷ 5)) : (1:0)
    dune(r) = scenario === :low ? 2.0 : r in gap ? 0.0 : 4.0
    topo = permutedims(reduce(hcat,
                              [barrier_profile(; MHW, dune=dune(r))
                               for r in 1:alongshore]))
    m, n = size(topo)
    spec1 = zeros(m, n)
    spec2 = zeros(m, n)
    toe = 60 + 76
    for I in CartesianIndices(topo)
        c = I[2]
        topo[I] > MHW + 0.25 || continue
        if toe - 10 <= c < toe + 30
            spec1[I] = 0.6 + 0.4rand(rng)
        elseif c >= toe + 30 && rand(rng) < 0.5
            spec2[I] = 0.4 + 0.4rand(rng)
        end
    end
    veg = clamp.(spec1 .+ spec2, 0.0, par(:maxveg))
    veg_eff = _meeb_gaussian(veg, par(:veg_sigma); mode=:constant)
    section = par(:section)
    raw = [Float64(_meeb_ocean_shoreline(topo, r, MHW)) for r in 1:m]
    xs = similar(raw)
    for lo in 1:section:m
        rr = lo:min(lo + section - 1, m)
        xs[rr] .= sum(raw[rr]) / length(rr)
    end
    xt = xs .- floor(Int, par(:D_sf) / par(:S_sf_eq))
    return CAProblem(model,
                     (; topo, spec1, spec2, veg, veg_eff, shoreline=repeat(xs, 1, n),
                      toe=repeat(xt, 1, n)), 50years, given; space=Grid(m, n), seed,
                     callback=barrier_callback(), kwargs...)
end

"""
    barrier_measures(u, MHW) -> NamedTuple

Measures of a barrier state `u` at mean high water `MHW`, averaged over its
rows: the height of the foredune crest above MHW (`crest`, m, at the columns
of [`barrier_crest`](@ref)), the fraction of rows with a gap in the foredune
(`gaps`), the subaerial volume (`volume`, m³ per m alongshore, above MHW), the
mean plant cover of the cells above MHW (`cover`) and the ocean shoreline
(`shoreline`, columns).
"""
function barrier_measures(u, MHW::Real)
    m = size(u.topo, 1)
    crest, ok = barrier_crest(u.topo, MHW)
    dry = u.topo .> MHW
    return (; crest=sum(u.topo[r, crest[r]] - MHW for r in 1:m) / m,
            gaps=count(!, ok) / m,
            volume=sum(max.(u.topo .- MHW, 0.0)) / m,
            cover=count(dry) > 0 ? sum(u.veg[dry]) / count(dry) : 0.0,
            shoreline=sum(view(u.shoreline, :, 1)) / m)
end

"""
    brier_skill_score(simulated, observed, baseline, mask) -> Float64

The Brier skill score of MEEB's hindcasts (Eq. 24 of [Reeves2025](@cite)): one
minus the mean squared error of `simulated` against `observed` over the cells
of the boolean `mask`, divided by that of `baseline` (no change: the initial
state).
"""
function brier_skill_score(simulated, observed, baseline, mask)
    mse = sum(abs2, simulated[mask] .- observed[mask]) / count(mask)
    ref = sum(abs2, baseline[mask] .- observed[mask]) / count(mask)
    return 1 - mse / ref
end

_register!(:barrier_island; build=build_barrier_island, problem=problem_barrier_island)
