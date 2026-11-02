# Liu & Estep (2025): a multiscale stochastic cellular automaton of the mountain pine beetle on a grid of
# hexagons. Within each cell, the beetles' numbers follow a birth-death process run several times and
# averaged (the microscale, the model's rule); between cells, part of the density of the more infested
# cells of each interaction neighbourhood is pooled and shared among cells drawn into the dispersal (the
# macroscale, a callback stored with the problem, because the pooling reaches across a whole cluster of
# infested cells). Cells are hexagons in the source's layout: m rows and n columns, odd columns half a
# cell lower, stored on a square lattice.

# The six neighbours of the cell in row i, column j: odd columns sit half a cell lower than even ones.
function _mpb_neighbours(i, j)
    return isodd(j) ?
           ((i - 1, j - 1), (i - 1, j), (i - 1, j + 1), (i, j - 1), (i, j + 1),
            (i + 1, j)) :
           ((i - 1, j), (i, j - 1), (i, j + 1), (i + 1, j - 1), (i + 1, j + 1), (i + 1, j))
end

# The number of steps between two cells of the hexagonal layout.
function _mpb_distance((i1, j1), (i2, j2))
    dq = j2 - j1
    dr = (i2 - cld(j2, 2)) - (i1 - cld(j1, 2))
    return (abs(dq) + abs(dr) + abs(dq + dr)) ÷ 2
end

# The birth and death rates per unit time of macroscale step `t` (the step from t - 1 to t).
function _mpb_rates(p, t)
    if p.schedule === :constant
        return p.birth, p.death
    elseif p.schedule === :periodic
        return p.birth * (1 + p.amplitude * sinpi(2t / p.period)), p.death
    end
    s, c = abs(sinpi(t / p.half_period)), abs(cospi(t / p.half_period))
    b1, d1 = t <= p.switch ? (p.b1_endemic, p.d1_endemic) : (p.b1_epidemic, p.d1_epidemic)
    b = p.b0 + b1 * (s < p.f1 ? s^p.shape : s^(1 / p.shape))
    d = p.d0 + d1 * (c < p.f2 ? c^p.shape : c^(1 / p.shape))
    return b, d
end

"""
    build_mpb_multiscale() -> CAModel

The microscale of the mountain pine beetle model of [Liu2025](@cite): each
cell's normalised beetle density (`density`, from 0 to 1) becomes the mean,
over `replicates` independent runs, of a birth-death process started from
`density × capacity` beetles and run for `micro_steps` steps of length
`1/micro_steps`, divided by `capacity`. In each step a population of `x`
gains a beetle with probability `b x (1 - x/capacity) / micro_steps`, loses
one with probability `d x / micro_steps`, and is set to 0 when below 1 at
the start of a step, as in the source's code. The rates `b` and `d` are
`birth` and `death` with `schedule = :constant`; with `:periodic`, the birth
rate is multiplied by `1 + amplitude sin(2π t / period)`; with `:outbreak`,
seasonal rates switch from endemic to epidemic amplitudes after step
`switch` (see the source's code). With `integer_population`, the starting
number is rounded, as the source's article defines it. Blocked cells
(`blocked`) hold no beetles.

The macroscale dispersal is [`mpb_disperse!`](@ref), a callback stored with
the problem.
"""
function build_mpb_multiscale()
    rule = LocalRule(Offsets((1, 0))) do c, nbrs, ctx
        p = ctx.p
        (c.density > 0 && !c.blocked) || return (density=0.0, blocked=c.blocked)
        b, d = _mpb_rates(p, ctx.t + 1)
        C, N = p.capacity, p.micro_steps
        x0 = c.density * C
        p.integer_population && (x0 = round(x0))
        total = 0.0
        for rep in 1:(p.replicates)
            x = x0
            for s in 1:N
                x < 1 && (x = 0.0)
                x == 0 && break
                pb = max(b * x * (1 - x / C), 0.0) / N
                pd = d * x / N
                u = uniform(ctx, (rep - 1) * N + s)
                x += u < pd ? -1 : u >= 1 - pb ? 1 : 0
            end
            total += x
        end
        return (density=total / (p.replicates * C), blocked=c.blocked)
    end
    states = [StateVariable(:density, 0.0;
                            description="beetle density, as a share of the cell's capacity"),
              StateVariable(:blocked, false;
                            description="a barrier cell, which holds no beetles and does not neighbour others")]
    parameters = [Parameter(:birth, 1.0; description="birth rate b"),
                  Parameter(:death, 0.2; description="death rate d"),
                  Parameter(:capacity, 100.0;
                            description="beetles a cell holds at most, C"),
                  Parameter(:micro_steps, 100; description="microscale steps a step, N"),
                  Parameter(:replicates, 20; description="microscale runs averaged, M"),
                  Parameter(:integer_population, false;
                            description="whether the starting number of beetles is rounded"),
                  Parameter(:schedule, :constant;
                            description=":constant, :periodic or :outbreak birth and death rates"),
                  Parameter(:amplitude, 2 / 3;
                            description="relative amplitude of periodic birth rates"),
                  Parameter(:period, 20;
                            description="steps in a year of periodic birth rates"),
                  Parameter(:b0, 0.8;
                            description="base birth rate of the outbreak schedule"),
                  Parameter(:d0, 0.5;
                            description="base death rate of the outbreak schedule"),
                  Parameter(:b1_endemic, 0.2;
                            description="seasonal birth amplitude before the switch"),
                  Parameter(:d1_endemic, 0.4;
                            description="seasonal death amplitude before the switch"),
                  Parameter(:b1_epidemic, 0.5;
                            description="seasonal birth amplitude after the switch"),
                  Parameter(:d1_epidemic, 0.2;
                            description="seasonal death amplitude after the switch"),
                  Parameter(:f1, 0.5;
                            description="level of |sin| switching the birth season"),
                  Parameter(:f2, 0.9;
                            description="level of |cos| switching the death season"),
                  Parameter(:shape, 10.0; description="exponent shaping the seasons"),
                  Parameter(:half_period, 20; description="steps between seasonal peaks"),
                  Parameter(:switch, 60;
                            description="last endemic step of the outbreak schedule"),
                  Parameter(:alpha, 1.0; description="exponent of the inclusion weight"),
                  Parameter(:beta, 1.0;
                            description="exponent of the inclusion probability"),
                  Parameter(:gamma1, 0.8;
                            description="share of the excess density of contributing cells dispersed"),
                  Parameter(:gamma2, 0.8;
                            description="share of the dispersed density reaching a receiving cell"),
                  Parameter(:gamma2_contributors, false;
                            description="whether contributing cells' share is also multiplied by gamma2 (the article's equations)"),
                  Parameter(:phi, 1 / 3;
                            description="base probability a receiving cell keeps its beetles"),
                  Parameter(:rho, 1 / 3;
                            description="increase of that probability with density"),
                  Parameter(:xi, 1.0;
                            description="scale of density in the survival probability"),
                  Parameter(:linkage, 0;
                            description="steps (2 or more) within which affected cells form one interaction neighbourhood; 0 for one neighbourhood of all (the source's code)"),
                  Parameter(:wind, :none;
                            description=":none, or the direction (:east, :west, :north, :south, :northeast, ...) beetles drift towards"),
                  Parameter(:wind_factor, 0.05;
                            description="factor on the inclusion probability of low-density cells upwind")]
    return CAModel(rule; name=:mpb_multiscale, states, parameters,
                   description="Liu & Estep (2025) multiscale stochastic CA of mountain pine beetle dispersal")
end

# Groups of affected cells (indices into `cells`) linked by chains of cells at most `linkage` apart.
function _mpb_groups(cells, linkage)
    linkage <= 0 && return [collect(eachindex(cells))]
    parent = collect(eachindex(cells))
    function root(k)
        while parent[k] != k
            parent[k] = parent[parent[k]]
            k = parent[k]
        end
        return k
    end
    for a in eachindex(cells), b in (a + 1):lastindex(cells)
        if _mpb_distance(cells[a], cells[b]) <= linkage
            parent[root(a)] = root(b)
        end
    end
    groups = Dict{Int,Vector{Int}}()
    for k in eachindex(cells)
        push!(get!(groups, root(k), Int[]), k)
    end
    return sort!(collect(values(groups)); by=first)
end

# Whether a cell at (row, col) lies upwind of `wind`, relative to the means of its neighbourhood.
function _mpb_upwind(wind, row, col, mrow, mcol)
    wind === :east && return col < mcol
    wind === :west && return col > mcol
    wind === :north && return row < mrow
    wind === :south && return row > mrow
    wind === :northeast && return col < mcol && row < mrow
    wind === :northwest && return col > mcol && row < mrow
    wind === :southeast && return col < mcol && row > mrow
    wind === :southwest && return col > mcol && row > mrow
    return false
end

"""
    mpb_disperse!(u, p, rng) -> u

The macroscale dispersal of [Liu2025](@cite), applied to the state `u` in
place after the microscale step. Affected cells (density above 0) form
interaction neighbourhoods: with `linkage = 0` all of them form one, as in the
source's code; otherwise those linked by chains of affected cells at most
`linkage` steps apart. Each neighbourhood, with all the immediate neighbours
of its cells (including auxiliary cells beyond the lattice, which can receive
but not keep beetles), disperses independently:

1. each cell is drawn into the dispersal with probability
   `(w^alpha d / (1 + w^alpha))^beta`, where `w` counts its affected
   neighbours (plus 1 if it is affected itself) and `d` is the largest
   absolute difference between its density and its neighbours'; with `wind`,
   the probability of cells below the neighbourhood's mean density and upwind
   of its centre is multiplied by `wind_factor`;
2. among the `k` cells drawn, those at or above their mean density `π̄`
   contribute `gamma1 (π - π̄)`, pooled into `π̃`;
3. the others receive `gamma2 π̃ / k` and keep their beetles with probability
   `phi + rho atan(xi π)` at their new density, or lose them all;
4. contributing cells become `π - gamma1 (π - π̄) + π̃ / k`, as in the source's
   code, or `π - gamma1 (π - π̄) + gamma2 π̃ / k` with `gamma2_contributors`,
   as its article writes it.

Blocked cells take no part and are no one's neighbours.
"""
function mpb_disperse!(u, p, rng)
    m, n = size(u.density)
    P = zeros(m + 2, n + 2)
    blocked(i, j) = 1 <= i <= m && 1 <= j <= n && u.blocked[i, j]
    for j in 1:n, i in 1:m
        P[i + 1, j + 1] = u.blocked[i, j] ? 0.0 : u.density[i, j]
    end
    dens(i, j) = 0 <= i <= m + 1 && 0 <= j <= n + 1 ? P[i + 1, j + 1] : 0.0
    affected = [(i, j) for j in 1:n for i in 1:m if P[i + 1, j + 1] > 0]
    isempty(affected) && return u
    nbrs(i, j) = Iterators.filter(c -> !blocked(c...), _mpb_neighbours(i, j))
    for group in _mpb_groups(affected, p.linkage)
        cells = Set{Tuple{Int,Int}}()
        for k in group
            push!(cells, affected[k])
            for c in nbrs(affected[k]...)
                0 <= c[1] <= m + 1 && 0 <= c[2] <= n + 1 && push!(cells, c)
            end
        end
        cells = sort!(collect(cells); by=c -> (c[2], c[1]))
        mrow = sum(first, cells) / length(cells)
        mcol = sum(last, cells) / length(cells)
        mdens = sum(c -> dens(c...), cells) / length(cells)
        drawn = Tuple{Int,Int}[]
        for (i, j) in cells
            π0 = dens(i, j)
            w = (π0 > 0) + count(c -> dens(c...) > 0, nbrs(i, j))
            d = maximum((abs(π0 - dens(c...)) for c in nbrs(i, j)); init=0.0)
            pI = (w^p.alpha * d / (1 + w^p.alpha))^p.beta
            p.wind !== :none && π0 < mdens && _mpb_upwind(p.wind, i, j, mrow, mcol) &&
                (pI *= p.wind_factor)
            rand(rng) < pI && push!(drawn, (i, j))
        end
        isempty(drawn) && continue
        k = length(drawn)
        avg = sum(c -> dens(c...), drawn) / k
        pool = p.gamma1 * sum(c -> max(dens(c...) - avg, 0.0), drawn)
        share = p.gamma2 * pool / k
        own = p.gamma2_contributors ? share : pool / k
        for (i, j) in drawn
            π0 = dens(i, j)
            if π0 >= avg
                P[i + 1, j + 1] = π0 - p.gamma1 * (π0 - avg) + own
            else
                π1 = π0 + share
                keep = rand(rng) < p.phi + p.rho * atan(p.xi * π1)
                P[i + 1, j + 1] = keep ? π1 : 0.0
            end
        end
    end
    for j in 1:n, i in 1:m
        u.density[i, j] = u.blocked[i, j] ? 0.0 : P[i + 1, j + 1]
    end
    return u
end

"""
    mpb_callback() -> DiscreteCallback

A callback that runs [`mpb_disperse!`](@ref) after each step's microscale,
with a random number generator seeded from the problem's seed and the time.
"""
function mpb_callback()
    affect! = function (integ)
        rng = Random.Xoshiro(hash((integ.prob.seed, Int(integ.t), :mpb)))
        mpb_disperse!(integ.u, integ.p, rng)
        return nothing
    end
    return DiscreteCallback((u, t, integ) -> true, affect!; save_positions=(false, false))
end

# The source's starting patterns on a 50 × 50 grid, as (column, first row, densities).
const MPB_PATTERNS = Dict(:behaviour => [(22, 23, [0.6, 0.5, 0.4, 0.7]),
                                         (23, 23, [0.5, 0.5, 0.4, 0.2, 0.8, 0.3]),
                                         (24, 22, [0.2, 0.3, 0.4, 0.7, 0.5, 0.4, 0.2]),
                                         (25, 21, [0.2, 0.5, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4]),
                                         (26, 23, [0.4, 0.7, 0.4, 0.3, 0.3]),
                                         (27, 23, [0.6, 0.3, 0.6, 0.6, 0.1]),
                                         (28, 23, [0.3, 0.5, 0.2, 0.1, 0.1]),
                                         (29, 23, [0.2, 0.6, 0.6, 0.4, 0.1])],
                          :endings => [(22, 23, [0.6, 0.5, 0.4, 0.7]),
                                       (23, 23, [0.5, 0.5, 0.4, 0.2, 0.8, 0.3]),
                                       (24, 22, [0.2, 0.3, 0.4, 0.7, 0.5, 0.4, 0.2]),
                                       (25, 21, [0.2, 0.5, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4]),
                                       (26, 23, [0.4, 0.7, 0.4, 0.3, 0.3]),
                                       (27, 23, [0.6, 0.3, 0.6, 0.6, 0.1]),
                                       (28, 23, [0.3, 0.5, 0.2, 0.1, 0.1]),
                                       (29, 23, [0.2, 0.6, 0.6, 0.4, 0.1]),
                                       (30, 23, [0.6, 0.5, 0.4, 0.3]),
                                       (31, 23, [0.6, 0.5, 0.4, 0.2])],
                          :sweep => [(22, 23, [0.6, 0.5, 0.4, 0.7]),
                                     (23, 23, [0.8, 0.7, 0.6, 0.1, 0.8, 0.1]),
                                     (24, 22, [0.7, 0.0, 0.8, 0.7, 0.6, 0.4, 0.2]),
                                     (25, 21, [0.2, 0.7, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8]),
                                     (26, 23, [0.8, 0.7, 0.6, 0.1, 0.8]),
                                     (27, 23, [0.6, 0.5, 0.4, 0.7]),
                                     (28, 23, [0.6, 0.5, 0.4, 0.7])])

"""
    mpb_pattern(name=:behaviour; n=50) -> Matrix{Float64}

One of the starting patterns of the source's code, on an `n × n` grid (rows
by columns): `:behaviour` (its Fig. 7), `:endings` (its Fig. 10) or `:sweep`
(its parameter sweeps), clusters of about 40 cells near the centre of a 50 × 50 grid, moved by an
even number of columns to stay near the centre of another size.
"""
function mpb_pattern(name::Symbol=:behaviour; n::Integer=50)
    haskey(MPB_PATTERNS, name) ||
        throw(ArgumentError("pattern must be :behaviour, :endings or :sweep, got :$name"))
    n >= 12 || throw(ArgumentError("the grid must be at least 12 cells wide, got $n"))
    shift = 2 * ((n - 50) ÷ 4)
    density = zeros(n, n)
    for (col, row, values) in MPB_PATTERNS[name], (k, v) in enumerate(values)
        density[row + k - 1 + shift, col + shift] = v
    end
    return density
end

"""
    problem_mpb_multiscale(; pattern=:behaviour, n=50, steps=150, blocked=nothing,
                           seed=0, pmap=(;), kwargs...) -> CAProblem

The setup of [Liu2025](@cite): an `n × n` grid of hexagons with absorbing
edges, starting from `pattern` (a name for [`mpb_pattern`](@ref), or a matrix
of densities), for `steps` steps. `blocked` is an optional `n × n` Boolean
matrix of barrier cells. The defaults are the article's Table 2: `birth = 1`,
`death = 0.2`, `gamma1 = gamma2 = 0.8`, `phi = rho = 1/3`, `capacity = 100`,
`micro_steps = 100`, `replicates = 20`.
"""
function problem_mpb_multiscale(; pattern=:behaviour, n::Integer=50, steps::Integer=150,
                                blocked=nothing, seed::Integer=0, pmap=(;), kwargs...)
    density = pattern isa Symbol ? mpb_pattern(pattern; n) : Float64.(pattern)
    size(density) == (n, n) ||
        throw(ArgumentError("the pattern must be $n × $n, got $(size(density))"))
    all(x -> 0 <= x <= 1, density) || throw(ArgumentError("densities must lie in [0, 1]"))
    mask = blocked === nothing ? fill(false, n, n) : collect(Bool, blocked)
    size(mask) == (n, n) || throw(ArgumentError("blocked must be $n × $n"))
    model = build_mpb_multiscale()
    given = (; (Symbol(k) => v for (k, v) in pairs(pmap))...)
    p = merge(defaults(model), given)
    p.schedule in (:constant, :periodic, :outbreak) ||
        throw(ArgumentError("schedule must be :constant, :periodic or :outbreak, got :$(p.schedule)"))
    (p.linkage == 0 || p.linkage >= 2) ||
        throw(ArgumentError("linkage must be 0 or at least 2, got $(p.linkage)"))
    u0 = (; density=density .* .!mask, blocked=mask)
    return CAProblem(model, u0, steps, given; space=Grid(n, n), seed,
                     callback=mpb_callback(), kwargs...)
end

"""
    mpb_complexity(u) -> Float64

The source's complexity of the pattern of affected cells (density above 0):
the edges between affected and unaffected cells (or the edge of the grid)
per affected cell, `(6 N - Σ n_a) / N` for `N` affected cells with `n_a`
affected neighbours each. `NaN` when no cell is affected.
"""
function mpb_complexity(u)
    D = u.density
    m, n = size(D)
    N, links = 0, 0
    for j in 1:n, i in 1:m
        D[i, j] > 0 || continue
        N += 1
        links += count(c -> 1 <= c[1] <= m && 1 <= c[2] <= n && D[c...] > 0,
                       _mpb_neighbours(i, j))
    end
    return N == 0 ? NaN : (6N - links) / N
end

"""
    mpb_radius(u) -> Float64

The radius of the circle with the area of the affected cells, as the source's
code computes it: `sqrt(6√3 N / π)` for `N` affected cells (hexagons of side
2). Its change from step to step is the source's rate of dispersion.
"""
mpb_radius(u) = sqrt(6 * sqrt(3) * count(>(0), u.density) / π)

_register!(:mpb_multiscale; build=build_mpb_multiscale, problem=problem_mpb_multiscale)
