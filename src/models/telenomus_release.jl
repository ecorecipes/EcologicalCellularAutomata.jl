# Weber, Garcia, Bueno, Bueno & Godoy (2022): releases of the egg parasitoid Telenomus podisi against
# the neotropical brown stink bug Euschistus heros in a soybean field of 1 m cells, in daily steps.
# Each cell holds female stink bugs (nymphs and adults, by age), up to 60 egg masses, and female
# wasps (by age). Ageing, death, hatching and emergence happen within each cell: that is the model's
# rule. Movement, egg laying, parasitism, colonisation and releases reach other cells and depend on
# what earlier individuals did that day: that is a step callback stored with the problem.

const TEL_NYMPH_AGES = 30    # nymphs become adults at age 30
const TEL_ADULT_AGES = 36    # adult stink bugs die at age 37
const TEL_WASP_AGES = 26     # wasps die at age 26
const TEL_SLOTS = 60         # egg masses a cell can hold
const TEL_HATCH_AGE = 5      # egg masses hatch at age 5
const TEL_EMERGE_AGE = 14    # wasps emerge from parasitised eggs at age 14
const TEL_K_NYMPH = 500
const TEL_K_ADULT = 300
const TEL_K_WASP = 550

"""
    StinkBugEggMass(eggs, age, para, para_age)

One of the up to 60 egg masses of a cell in the `telenomus_release` model:
`eggs` live unparasitised stink bug eggs laid `age` days ago (0 when none
remain), and `para` parasitoid immatures of mean age `para_age` (0 when none).
A place holds a new egg mass only when both ages are 0.
"""
struct StinkBugEggMass
    eggs::Int8
    age::Int8
    para::Int8
    para_age::Int8
end

const TEL_NO_MASS = StinkBugEggMass(Int8(0), Int8(0), Int8(0), Int8(0))

# Binomial draw by inversion from one counter-based uniform; p is small and n at most a few hundred.
function _tel_binomial(ctx, n::Integer, p::Real, k::Integer)
    (n <= 0 || p <= 0) && return 0
    p >= 1 && return Int(n)
    u = uniform(ctx, k)
    q = 1 - p
    f = q^n
    F = f
    x = 0
    while u > F && x < n
        f *= (n - x) / (x + 1) * (p / q)
        x += 1
        F += f
    end
    return x
end

"""
    telenomus_dispersal_steps(mean, k) -> NTuple{4,Float64}

The cumulative probabilities that a stink bug moves 0, 1, 2 or 3 cells along
each axis in a day: a negative binomial distribution with the given mean and
aggregation `k`, truncated at 3. The source's tables are
`telenomus_dispersal_steps(0.5, 2.8433)` for nymphs and
`telenomus_dispersal_steps(0.5, 3.1728)` for adults (to four decimals).
"""
function telenomus_dispersal_steps(mean::Real, k::Real)
    mean > 0 && k > 0 ||
        throw(ArgumentError("mean and k must be positive, got $mean and $k"))
    r = mean / (k + mean)
    p0 = (k / (k + mean))^k
    p1 = p0 * k * r
    p2 = p1 * (k + 1) / 2 * r
    return (p0, p0 + p1, p0 + p1 + p2, 1.0)
end

"""
    build_telenomus_release() -> CAModel

The within-cell dynamics of the soybean field of [Weber2022](@cite), in daily
steps on 1 m cells. Each cell holds female stink bugs (*Euschistus heros*) as
counts by age, `nymphs` (ages 1 to 30) and `adults` (ages 1 to 36), up to 60
egg masses ([`StinkBugEggMass`](@ref)) in `masses`, and female wasps
(*Telenomus podisi*) as counts by age in `wasps` (ages 1 to 26). Each day, in
every cell:

- adult stink bugs age a day (not on the first day) and die with probability
  `mort_adult` or at age 37; wasps die with probability `mort_wasp` or at
  age 26;
- nymphs of age 30 become adults of age 1 while the cell has room (300
  adults); the others age a day and die with probability `mort_nymph`;
- egg masses of age 5 hatch `⌊eggs × sex_ratio × egg_survival⌋` female nymphs
  (3 for an intact mass of 11) while the cell has room (500 nymphs), and the
  others age a day; parasitoid immatures of age 14 are removed or, with
  `offspring = true`, give `⌊para × sex_ratio_wasp × immature_survival⌋`
  female wasps while the cell has room (550 wasps), and the others age a day.

Movement, egg laying, parasitism, colonisation and releases are
[`telenomus_callback`](@ref).
"""
function build_telenomus_release()
    rule = LocalRule(Offsets((1, 0))) do c, nbrs, ctx
        p = ctx.p
        shift = ctx.t != 0
        adults0 = ntuple(Val(TEL_ADULT_AGES)) do a
            src = shift ? a - 1 : a
            src < 1 && return Int16(0)
            n = c.adults[src]
            return Int16(n - _tel_binomial(ctx, n, p.mort_adult, src))
        end
        wasps0 = ntuple(Val(TEL_WASP_AGES)) do a
            a == TEL_WASP_AGES && return Int16(0)
            n = c.wasps[a]
            return Int16(n - _tel_binomial(ctx, n, p.mort_wasp, TEL_ADULT_AGES + a))
        end
        nymphs0 = ntuple(Val(TEL_NYMPH_AGES)) do a
            a == 1 && return Int16(0)
            n = c.nymphs[a - 1]
            return Int16(n - _tel_binomial(ctx, n, p.mort_nymph, 70 + a))
        end
        emerging = Int(c.nymphs[TEL_NYMPH_AGES])
        new_adults = min(emerging, TEL_K_ADULT - sum(Int, adults0))
        adults = Base.setindex(adults0, adults0[1] + Int16(new_adults), 1)
        hatching = 0
        wasp_brood = 0
        for m in c.masses
            m.age == TEL_HATCH_AGE &&
                (hatching += floor(Int, m.eggs * p.sex_ratio * p.egg_survival))
            m.para_age == TEL_EMERGE_AGE &&
                (wasp_brood += floor(Int, m.para * p.sex_ratio_wasp * p.immature_survival))
        end
        hatched = min(hatching, TEL_K_NYMPH - sum(Int, nymphs0))
        nymphs = Base.setindex(nymphs0, Int16(hatched), 1)
        emerged = p.offspring ? min(wasp_brood, TEL_K_WASP - sum(Int, wasps0)) : 0
        wasps = Base.setindex(wasps0, wasps0[1] + Int16(emerged), 1)
        masses = map(c.masses) do m
            eggs, age = m.age == TEL_HATCH_AGE ? (Int8(0), Int8(0)) :
                        m.age > 0 ? (m.eggs, m.age + Int8(1)) : (m.eggs, m.age)
            para, para_age = m.para_age == TEL_EMERGE_AGE ? (Int8(0), Int8(0)) :
                             m.para_age > 0 ? (m.para, m.para_age + Int8(1)) :
                             (m.para, m.para_age)
            return StinkBugEggMass(eggs, age, para, para_age)
        end
        return (; nymphs, adults, masses, wasps)
    end
    zeros_of(n) = ntuple(_ -> Int16(0), Val(n))
    states = [StateVariable(:nymphs, zeros_of(TEL_NYMPH_AGES);
                            description="female stink bug nymphs by age (days 1-30)"),
              StateVariable(:adults, zeros_of(TEL_ADULT_AGES);
                            description="female adult stink bugs by age (days 1-36)"),
              StateVariable(:masses, ntuple(_ -> TEL_NO_MASS, Val(TEL_SLOTS));
                            description="the cell's 60 places for stink bug egg masses"),
              StateVariable(:wasps, zeros_of(TEL_WASP_AGES);
                            description="female T. podisi by age (days 1-26)")]
    parameters = [Parameter(:mort_adult, 0.03;
                            description="daily death rate of adult stink bugs"),
                  Parameter(:mort_nymph, 0.014; description="daily death rate of nymphs"),
                  Parameter(:egg_survival, 1 - 0.37; description="stink bug egg viability"),
                  Parameter(:sex_ratio, 0.5; description="share of female stink bugs"),
                  Parameter(:p_oviposit, 0.7;
                            description="daily probability a mature female lays an egg mass"),
                  Parameter(:eggs_per_mass, 11; description="eggs in an egg mass"),
                  Parameter(:maturity, 13; description="age at which stink bugs lay eggs"),
                  Parameter(:dispersal_age, 5;
                            description="nymphs up to this age do not move"),
                  Parameter(:nymph_steps, (0.6309, 0.8992, 0.9763, 1.0);
                            description="cumulative probabilities of 0-3 cells per axis, nymphs"),
                  Parameter(:adult_steps, (0.6286, 0.9001, 0.9772, 1.0);
                            description="cumulative probabilities of 0-3 cells per axis, adults"),
                  Parameter(:mort_wasp, 0.0566; description="daily death rate of wasps"),
                  Parameter(:immature_survival, 1 - 0.1757;
                            description="survival of parasitoid immatures"),
                  Parameter(:sex_ratio_wasp, 0.7; description="share of female wasps"),
                  Parameter(:offspring, false;
                            description="whether parasitised eggs give new wasps"),
                  Parameter(:radius, 8; description="wasp search radius (m)"),
                  Parameter(:fecundity, 13.62737;
                            description="eggs a wasp lays on its second day"),
                  Parameter(:fecundity_decay, 0.10221;
                            description="daily decay of wasp fecundity"),
                  Parameter(:colonists, 8000;
                            description="adult stink bugs arriving in the first days"),
                  Parameter(:colonisation_days, 10;
                            description="days over which colonists arrive"),
                  Parameter(:colonised, 40;
                            description="half-width of the central square colonists land in"),
                  Parameter(:releases, 0; description="number of wasp releases"),
                  Parameter(:per_ha, 5000; description="female wasps per release"),
                  Parameter(:spacing, 25;
                            description="distance between release points or strips (m)"),
                  Parameter(:mode, :strips; description="release in :strips or at :points"),
                  Parameter(:first_release, 4; description="day of the first release"),
                  Parameter(:interval, 7; description="days between releases"),
                  Parameter(:release_area, 100;
                            description="side of the central square of releases (m)")]
    return CAModel(rule; name=:telenomus_release, states, parameters,
                   description="Weber et al. (2022) releases of Telenomus podisi against Euschistus heros")
end

_tel_inside(x, y, n1, n2) = 1 <= x <= n1 && 1 <= y <= n2

# One axis of a stink bug's daily move: a direction and a distance from the cumulative table.
function _tel_axis_step(rng, steps)
    d = 2 * rand(rng) - 1
    s = d < 0 ? -1 : 1
    d = abs(d)
    for k in 0:2
        d <= steps[k + 1] && return s * k
    end
    return 3s
end

function _tel_bug_destination(rng, i, j, steps, cur, new, cap)
    n1, n2 = size(cur)
    while true
        dx = _tel_axis_step(rng, steps)
        dy = _tel_axis_step(rng, steps)
        x, y = i + dx, j + dy
        _tel_inside(x, y, n1, n2) || continue
        dx == 0 && dy == 0 && return (i, j)
        cur[x, y] + new[x, y] < cap && return (x, y)
    end
end

const TEL_MOORE = ((-1, -1), (0, -1), (1, -1), (-1, 0), (0, 0), (1, 0), (-1, 1), (0, 1),
                   (1, 1))

function _tel_lay!(M, rng, order, x, y, eggs)
    _, n1, n2 = size(M)
    Random.shuffle!(rng, order)
    for o in order
        a, b = x + TEL_MOORE[o][1], y + TEL_MOORE[o][2]
        _tel_inside(a, b, n1, n2) || continue
        for s in axes(M, 1)
            m = M[s, a, b]
            if m.age == 0 && m.para_age == 0
                M[s, a, b] = StinkBugEggMass(Int8(eggs), Int8(1), m.para, m.para_age)
                return true
            end
        end
    end
    return false
end

function _tel_move_stink_bugs!(N, A, M, p, rng)
    _, n1, n2 = size(N)
    newN = zero(N)
    newA = zero(A)
    curN = Int[sum(Int, view(N, :, i, j)) for i in 1:n1, j in 1:n2]
    curA = Int[sum(Int, view(A, :, i, j)) for i in 1:n1, j in 1:n2]
    totN = zeros(Int, n1, n2)
    totA = zeros(Int, n1, n2)
    order = collect(1:9)
    for j in 1:n2, i in 1:n1
        for a in axes(N, 1), _ in 1:N[a, i, j]
            x, y = a <= p.dispersal_age ? (i, j) :
                   _tel_bug_destination(rng, i, j, p.nymph_steps, curN, totN, TEL_K_NYMPH)
            curN[i, j] -= 1
            newN[a, x, y] += 1
            totN[x, y] += 1
        end
        for a in axes(A, 1), _ in 1:A[a, i, j]
            x, y = _tel_bug_destination(rng, i, j, p.adult_steps, curA, totA, TEL_K_ADULT)
            curA[i, j] -= 1
            newA[a, x, y] += 1
            totA[x, y] += 1
            a >= p.maturity && rand(rng) <= p.p_oviposit &&
                _tel_lay!(M, rng, order, x, y, p.eggs_per_mass)
        end
    end
    N .= newN
    A .= newA
    return nothing
end

# Parasitise the egg masses of one cell in order; returns the eggs left to lay.
function _tel_parasitise_cell!(M, x, y, budget)
    for s in axes(M, 1)
        m = M[s, x, y]
        m.eggs > 0 || continue
        k = min(budget, Int(m.eggs))
        para_age = (Int(m.para_age) * m.para + k) ÷ (m.para + k)
        left = m.eggs - k
        M[s, x, y] = StinkBugEggMass(Int8(left), left == 0 ? Int8(0) : m.age,
                                     Int8(m.para + k), Int8(para_age))
        budget -= k
        budget == 0 && return 0
    end
    return budget
end

function _tel_parasitise!(M, rng, order, x, y, budget)
    _, n1, n2 = size(M)
    budget = _tel_parasitise_cell!(M, x, y, budget)
    budget == 0 && return nothing
    Random.shuffle!(rng, order)
    for o in order
        a, b = x + TEL_MOORE[o][1], y + TEL_MOORE[o][2]
        (a, b) == (x, y) && continue
        _tel_inside(a, b, n1, n2) || continue
        budget = _tel_parasitise_cell!(M, a, b, budget)
        budget == 0 && return nothing
    end
    return nothing
end

function _tel_wasp_destination(rng, i, j, cands, r, cur, new)
    n1, n2 = size(cur)
    left = length(cands)
    tries = left > 1 ? left : 1
    while true
        if left == 0
            dx, dy = rand(rng, (-r):r), rand(rng, (-r):r)
            x, y = i + dx, j + dy
            (isqrt(dx^2 + dy^2) > r || !_tel_inside(x, y, n1, n2)) && continue
        else
            x, y = cands[rand(rng, 1:length(cands))]
            tries == 1 && (left = 0)
        end
        (x, y) == (i, j) && return (i, j)
        cur[x, y] + new[x, y] < TEL_K_WASP && return (x, y)
        tries = max(tries - 1, 1)
    end
end

function _tel_move_wasps!(W, M, p, rng)
    _, n1, n2 = size(W)
    newW = zero(W)
    cur = Int[sum(Int, view(W, :, i, j)) for i in 1:n1, j in 1:n2]
    tot = zeros(Int, n1, n2)
    cells = [(i, j) for j in 1:n2 for i in 1:n1 if cur[i, j] > 0]
    Random.shuffle!(rng, cells)
    r = p.radius
    disc = [(dx, dy) for dx in (-r):r for dy in (-r):r if dx^2 + dy^2 <= r^2]
    cands = Tuple{Int,Int}[]
    order = collect(1:9)
    for (i, j) in cells
        empty!(cands)
        for (dx, dy) in disc
            x, y = i + dx, j + dy
            _tel_inside(x, y, n1, n2) && any(m -> m.eggs > 0, view(M, :, x, y)) &&
                push!(cands, (x, y))
        end
        for a in axes(W, 1), _ in 1:W[a, i, j]
            x, y = _tel_wasp_destination(rng, i, j, cands, r, cur, tot)
            cur[i, j] -= 1
            newW[a, x, y] += 1
            tot[x, y] += 1
            if 1 < a < TEL_WASP_AGES
                budget = floor(Int, p.fecundity * exp(-p.fecundity_decay * (a - 1)))
                budget > 0 && _tel_parasitise!(M, rng, order, x, y, budget)
            end
        end
    end
    # wasps age a day
    W[1, :, :] .= 0
    W[2:end, :, :] .= view(newW, 1:(size(W, 1) - 1), :, :)
    return nothing
end

function _tel_colonise!(A, rng, n, p)
    _, n1, n2 = size(A)
    ci, cj = n1 ÷ 2, n2 ÷ 2
    h = p.colonised
    while n > 0
        i, j = ci + rand(rng, (-h):h), cj + rand(rng, (-h):h)
        _tel_inside(i, j, n1, n2) || continue
        sum(Int, view(A, :, i, j)) < TEL_K_ADULT || continue
        A[rand(rng, 1:(p.maturity)), i, j] += 1
        n -= 1
    end
    return nothing
end

"""
    telenomus_release_cells(p, n1, n2) -> Vector{Tuple{Int,Int}}

The cells that receive wasps at a release of the strategy in the parameters
`p` (`mode`, `spacing`, `release_area`) on an `n1 × n2` field: a grid of points
`spacing` m apart, or strips `spacing` m apart, centred in the central square
of side `release_area`.
"""
function telenomus_release_cells(p, n1::Integer, n2::Integer)
    A, s = p.release_area, p.spacing
    1 <= s <= A || throw(ArgumentError("spacing must be between 1 and $A m, got $s"))
    p.mode in (:strips, :points) ||
        throw(ArgumentError("mode must be :strips or :points, got $(repr(p.mode))"))
    rows = A ÷ s
    cols, c0, cstep = p.mode === :strips ? (A, 0, 1) : (rows, s ÷ 2, s)
    r0 = s ÷ 2
    first_row = (n1 - A) ÷ 2 + 1 + r0 + (A - 2r0 - s * (rows - 1)) ÷ 2
    first_col = (n2 - A) ÷ 2 + 1 + c0 + (A - 2c0 - cstep * (cols - 1)) ÷ 2
    return [(first_row + s * x, first_col + cstep * y) for x in 0:(rows - 1)
            for y in 0:(cols - 1)]
end

function _tel_release!(W, p)
    _, n1, n2 = size(W)
    cells = telenomus_release_cells(p, n1, n2)
    each = p.per_ha / length(cells) * p.release_area^2 / 10_000
    whole = floor(Int, each)
    residue = each - whole
    carried = residue
    spill = p.mode === :strips ? (1, 0) : (3, 3)
    for (i, j) in cells
        k = whole + floor(Int, carried)
        carried += residue - floor(Int, carried)
        for a in i:min(i + spill[1], n1), b in j:min(j + spill[2], n2)
            k == 0 && break
            put = min(k, TEL_K_WASP - sum(Int, view(W, :, a, b)))
            put > 0 || continue
            W[1, a, b] += put
            k -= put
        end
    end
    return nothing
end

function _tel_release_day(p, day)
    p.releases >= 1 || return false
    d = day - p.first_release
    return day > 0 && d >= 0 && d % p.interval == 0 && d ÷ p.interval < p.releases
end

"""
    telenomus_step!(u, p, day, rng) -> u

The part of a day of [Weber2022](@cite) that reaches beyond a cell, done on
the state `u` after the day's within-cell step, so that `u` becomes the field
at the start of day `day`:

1. Stink bugs move, cell by cell: each nymph older than `dispersal_age` and
   each adult moves along each axis by a distance drawn from `nymph_steps` or
   `adult_steps`, staying put on a draw of 0 and 0 and otherwise taking a
   free place in the cell drawn (drawing again if it is full or outside the
   field). A mature adult then lays an egg mass with probability
   `p_oviposit` in the first free place among the 3 × 3 cells around it, in
   random order.
2. Wasps move, cell by cell in random order: each picks at random a cell
   within `radius` that holds unparasitised eggs, or any cell within
   `radius` when none does, and, from its second day, parasitises
   `⌊fecundity e^(-fecundity_decay (age - 1))⌋` eggs in its new cell and then
   in the 8 cells around it, in random order. Wasps then age a day.
3. Colonists arrive (in the first `colonisation_days` days) in the central
   square, aged 1 to 13 days, and wasps are released on release days at
   [`telenomus_release_cells`](@ref), filling a cell to 550 and spilling over
   into the 4 × 4 block (points) or the next row (strips).
"""
function telenomus_step!(u, p, day::Integer, rng)
    N = reinterpret(reshape, Int16, u.nymphs)
    A = reinterpret(reshape, Int16, u.adults)
    W = reinterpret(reshape, Int16, u.wasps)
    M = reinterpret(reshape, StinkBugEggMass, u.masses)
    _tel_move_stink_bugs!(N, A, M, p, rng)
    _tel_move_wasps!(W, M, p, rng)
    day < p.colonisation_days &&
        _tel_colonise!(A, rng, p.colonists ÷ p.colonisation_days, p)
    _tel_release_day(p, day) && _tel_release!(W, p)
    return u
end

"""
    telenomus_callback() -> DiscreteCallback

A callback that runs [`telenomus_step!`](@ref) after each day's within-cell
step, with a random number generator seeded from the problem's seed and the
day.
"""
function telenomus_callback()
    affect! = function (integ)
        day = Int(integ.t)
        rng = Random.Xoshiro(hash((integ.prob.seed, day, :telenomus)))
        telenomus_step!(integ.u, integ.p, day, rng)
        return nothing
    end
    return DiscreteCallback((u, t, integ) -> true, affect!; save_positions=(false, false))
end

"""
    problem_telenomus_release(; releases=0, per_ha=5000, spacing=25, mode=:strips,
                              offspring=false, n=200, days=80, seed=0, pmap=(;),
                              kwargs...) -> CAProblem

The setup of [Weber2022](@cite): an `n × n` soybean field of 1 m cells (4 ha),
colonised by 8000 adult female stink bugs over the first 10 days in its
central 81 × 81 m, for 80 days. `releases` releases of `per_ha` female wasps,
from day 4 and a week apart, go in `mode` (`:strips` or `:points`) `spacing`
m apart in the central hectare; `offspring = false` (the source's main
experiment) removes parasitised eggs at the wasps' emergence. The state at
time `t` is the field at the start of day `t`, after that day's arrivals.
Nothing lies beyond the field.
"""
function problem_telenomus_release(; releases::Integer=0, per_ha::Real=5000,
                                   spacing::Integer=25, mode::Symbol=:strips,
                                   offspring::Bool=false, n::Integer=200, days::Integer=80,
                                   seed::Integer=0, pmap=(;), kwargs...)
    releases >= 0 || throw(ArgumentError("releases must not be negative, got $releases"))
    per_ha >= 0 || throw(ArgumentError("per_ha must not be negative, got $per_ha"))
    model = build_telenomus_release()
    given = merge((; releases, per_ha, spacing, mode, offspring),
                  (; (Symbol(k) => v for (k, v) in pairs(pmap))...))
    p = merge(defaults(model), given)
    n >= p.release_area + 4 ||
        throw(ArgumentError("the field must be at least $(p.release_area + 4) m wide, got $n"))
    telenomus_release_cells(p, n, n)
    zeros_of(k) = ntuple(_ -> Int16(0), Val(k))
    u0 = (nymphs=fill(zeros_of(TEL_NYMPH_AGES), n, n),
          adults=fill(zeros_of(TEL_ADULT_AGES), n, n),
          masses=fill(ntuple(_ -> TEL_NO_MASS, Val(TEL_SLOTS)), n, n),
          wasps=fill(zeros_of(TEL_WASP_AGES), n, n))
    rng = Random.Xoshiro(hash((seed, 0, :telenomus)))
    _tel_colonise!(reinterpret(reshape, Int16, u0.adults), rng,
                   p.colonists ÷ p.colonisation_days, p)
    return CAProblem(model, u0, days, given; space=Grid(n, n; boundary=Reflecting()), seed,
                     callback=telenomus_callback(), kwargs...)
end

"""
    telenomus_et(u; area=100) -> Float64

The source's measure against the economic threshold of 2: stink bugs of at
least 0.5 cm (adults and nymphs older than 12 days) per m² in the central
square of side `area`. As in the source, the count covers one more row and
column than the square (101 × 101 cells for the default) and is divided by
`area²`.
"""
function telenomus_et(u; area::Integer=100)
    n1, n2 = size(u.adults)
    r0, c0 = (n1 - area) ÷ 2 + 1, (n2 - area) ÷ 2 + 1
    total = 0
    for j in c0:min(c0 + area, n2), i in r0:min(r0 + area, n1)
        total += sum(Int, u.adults[i, j]) + sum(Int, u.nymphs[i, j][13:end])
    end
    return total / area^2
end

"""
    telenomus_counts(u) -> NamedTuple

Whole-field totals: stink bug `eggs`, `nymphs` and `adults`, parasitoid
`immatures` and adult `wasps`.
"""
function telenomus_counts(u)
    return (eggs=sum(c -> sum(m -> Int(m.eggs), c), u.masses),
            nymphs=sum(c -> sum(Int, c), u.nymphs), adults=sum(c -> sum(Int, c), u.adults),
            immatures=sum(c -> sum(m -> Int(m.para), c), u.masses),
            wasps=sum(c -> sum(Int, c), u.wasps))
end

_register!(:telenomus_release; build=build_telenomus_release,
           problem=problem_telenomus_release)
