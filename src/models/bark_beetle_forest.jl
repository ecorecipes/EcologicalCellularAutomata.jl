# Königswieser & Jahn (2020): a spruce forest on a grid of 1 m cells, each holding at most one tree,
# which ages, grows, dies and spreads saplings, with contaminated ground, strip harvest, drought and bark
# beetles, in yearly steps. The rule is each tree's own year: removal of the trees felled or killed by
# beetles the year before, ageing, growth, and death by age or contaminated ground. Spreading, harvest,
# infestation and the beetles' spread are a callback stored with the problem, because each depends on
# what was done before it in the same year.

const BBF_STANDING = Int8(0)
const BBF_FELLED = Int8(1)
const BBF_KILLED = Int8(2)

_bbf_drought(p, year) = any(s -> s < year <= s + p.drought_length, p.droughts)

# The mass a tree of `age` (at the start of the year) adds in a year, in kg.
function _bbf_growth(age, p)
    return age < p.mature_age ? p.growth_young : age < p.final_age ? p.growth_old : 0.0
end

# The yearly death probability of a tree that has just reached `age`.
function _bbf_mortality(age, p)
    return age <= p.sapling_age ? p.mort_sapling :
           age <= p.old_age ? p.mort_adult : p.mort_old
end

# The probability that a tree with `n` neighbours spreads: by class of `n`.
function _bbf_spread_probability(n, p)
    for (limit, prob) in zip(p.spread_limits, p.spread_probabilities)
        n <= limit && return prob
    end
    return 0.0
end

"""
    build_bark_beetle_forest() -> CAModel

The spruce forest model of [Konigswieser2020](@cite), in yearly steps. Each
cell of 1 m² holds at most one tree (`tree`), with an age in years (`age`), a
dry wood mass in kg (`mass`), a number of bark beetles (`beetles`), a fate
for the coming year (`fate`: 0 standing, 1 felled, 2 killed by beetles) and
a fixed flag for contaminated ground (`contaminated`). Each year, in every
cell:

- a tree felled or killed by beetles the year before is removed, with its
  beetles;
- a tree ages a year and gains `growth_young` kg a year until `mature_age`,
  then `growth_old` kg until `final_age`;
- it dies with probability `mort_sapling` up to `sapling_age`, `mort_adult`
  up to `old_age`, and `mort_old` beyond, and at once on contaminated ground
  past `contaminated_age`.

In a drought year (`droughts`, each lasting `drought_length` years) growth
is multiplied by `drought_growth`, mortality by `drought_mortality` and the
chance of spreading by `drought_spread`.

Spreading, strip harvest, infestation and the beetles' spread are
[`bark_beetle_step!`](@ref), a callback stored with the problem.
"""
function build_bark_beetle_forest()
    rule = LocalRule(Offsets((1, 0))) do c, nbrs, ctx
        p = ctx.p
        gone = (tree=false, age=Int16(0), mass=0.0, beetles=0.0, fate=BBF_STANDING,
                contaminated=c.contaminated)
        (c.tree && c.fate == BBF_STANDING) || return gone
        dry = _bbf_drought(p, ctx.t + 1)
        age = c.age + Int16(1)
        pd = _bbf_mortality(age, p) * (dry ? p.drought_mortality : 1.0)
        dies = bernoulli(ctx, pd, 1) || (c.contaminated && age > p.contaminated_age)
        dies && return gone
        mass = c.mass + _bbf_growth(c.age, p) * (dry ? p.drought_growth : 1.0)
        return (tree=true, age, mass, beetles=c.beetles, fate=BBF_STANDING,
                contaminated=c.contaminated)
    end
    states = [StateVariable(:tree, false; description="a tree stands in the cell"),
              StateVariable(:age, Int16(0); description="age of the tree (years)"),
              StateVariable(:mass, 0.0; description="dry wood mass of the tree (kg)"),
              StateVariable(:beetles, 0.0; description="bark beetles on the tree"),
              StateVariable(:fate, BBF_STANDING; domain=Int8(0):Int8(2),
                            description="0 standing, 1 felled, 2 killed by beetles; removed next year"),
              StateVariable(:contaminated, false;
                            description="contaminated ground, where trees die past a young age")]
    parameters = [Parameter(:growth_young, 28.0;
                            description="mass gained a year until mature_age (kg)"),
                  Parameter(:growth_old, 14.0;
                            description="mass gained a year from mature_age to final_age (kg)"),
                  Parameter(:mature_age, 70; description="age of slower growth (years)"),
                  Parameter(:final_age, 140; description="age growth stops (years)"),
                  Parameter(:mort_sapling, 0.05;
                            description="yearly death probability up to sapling_age"),
                  Parameter(:mort_adult, 0.01;
                            description="yearly death probability up to old_age"),
                  Parameter(:mort_old, 0.2;
                            description="yearly death probability beyond old_age"),
                  Parameter(:sapling_age, 3;
                            description="last age of high sapling mortality"),
                  Parameter(:old_age, 150; description="last age of low mortality"),
                  Parameter(:contaminated_age, 3;
                            description="trees older than this die on contaminated ground"),
                  Parameter(:neighbour_radius, 10;
                            description="half-width of the square counting neighbours (m)"),
                  Parameter(:spread_limits, (0, 3, 9, 50);
                            description="upper neighbour counts of the spreading classes"),
                  Parameter(:spread_probabilities, (0.0, 0.5, 1.0, 0.5);
                            description="chance of spreading in each class; none beyond the last"),
                  Parameter(:spread_interval, 2; description="years between spreading"),
                  Parameter(:sapling_radius, 6;
                            description="half-width of the square a sapling lands in (m)"),
                  Parameter(:spacing, 2.0;
                            description="distance within which a sapling's site must be free of trees (m)"),
                  Parameter(:site_limit, 50;
                            description="most trees allowed around a sapling's site (negative for no limit)"),
                  Parameter(:droughts, ();
                            description="years after which a drought starts"),
                  Parameter(:drought_length, 10; description="years a drought lasts"),
                  Parameter(:drought_mortality, 1.2;
                            description="factor on mortality in a drought"),
                  Parameter(:drought_growth, 0.5;
                            description="factor on growth in a drought"),
                  Parameter(:drought_spread, 0.5;
                            description="factor on the chance of spreading in a drought"),
                  Parameter(:harvest, false; description="whether strips are harvested"),
                  Parameter(:harvest_start, 100;
                            description="years after which harvest starts"),
                  Parameter(:harvest_years, 150; description="years harvest lasts"),
                  Parameter(:harvest_age, 90;
                            description="trees older than this are felled"),
                  Parameter(:strips, 10;
                            description="strips harvested in turn, one a year"),
                  Parameter(:infestation, -1;
                            description="year beetles are brought in (negative for never)"),
                  Parameter(:infest_rows, (100, 120);
                            description="rows of the area beetles are brought into"),
                  Parameter(:infest_cols, (100, 120);
                            description="columns of the area beetles are brought into"),
                  Parameter(:infest_draws, 100;
                            description="trees drawn, with replacement, to receive beetles"),
                  Parameter(:infest_beetles, 20.0;
                            description="beetles added at each draw"),
                  Parameter(:beetle_growth, 10.0;
                            description="factor on an infested tree's beetles at each spread"),
                  Parameter(:beetle_spreads, 3;
                            description="spreads of the beetles a year"),
                  Parameter(:beetle_radius, 10;
                            description="half-width of the square beetles spread over (m)"),
                  Parameter(:whole_beetles, true;
                            description="whether beetles are whole, each tree's share rounded down"),
                  Parameter(:beetle_kill, 200.0;
                            description="beetles that kill a tree at the end of the year")]
    return CAModel(rule; name=:bark_beetle_forest, states, parameters,
                   description="Königswieser & Jahn (2020) spruce forest with harvest, drought and bark beetles")
end

# Summed-area table of `A` into `S` (one larger on each axis).
function _bbf_integral!(S, A)
    n1, n2 = size(A)
    @inbounds for j in 1:n2, i in 1:n1
        S[i + 1, j + 1] = A[i, j] + S[i, j + 1] + S[i + 1, j] - S[i, j]
    end
    return S
end

# The sum of the table's array over the square of half-width `r` around (i, j), clipped to the grid.
@inline function _bbf_box(S, i, j, r)
    n1, n2 = size(S) .- 1
    i0, i1, j0, j1 = max(i - r, 1), min(i + r, n1), max(j - r, 1), min(j + r, n2)
    return @inbounds S[i1 + 1, j1 + 1] - S[i0, j1 + 1] - S[i1 + 1, j0] + S[i0, j0]
end

# The trees within the square of half-width `r` around (x, y).
function _bbf_count(tree, x, y, r)
    n1, n2 = size(tree)
    n = 0
    @inbounds for j in max(y - r, 1):min(y + r, n2), i in max(x - r, 1):min(x + r, n1)
        n += tree[i, j]
    end
    return n
end

# Whether no tree stands within distance `d` of (x, y).
function _bbf_free(tree, x, y, d)
    n1, n2 = size(tree)
    h = floor(Int, d)
    @inbounds for j in max(y - h, 1):min(y + h, n2), i in max(x - h, 1):min(x + h, n1)
        tree[i, j] && (i - x)^2 + (j - y)^2 <= d^2 && return false
    end
    return true
end

function _bbf_spread!(u, p, year, rng)
    tree = u.tree
    n1, n2 = size(tree)
    S = _bbf_integral!(zeros(Int, n1 + 1, n2 + 1), tree)
    parents = findall(tree)
    counts = [_bbf_box(S, I[1], I[2], p.neighbour_radius) - 1 for I in parents]
    factor = _bbf_drought(p, year) ? p.drought_spread : 1.0
    r, d = p.sapling_radius, p.spacing
    for (I, n) in zip(parents, counts)
        rand(rng) < _bbf_spread_probability(n, p) * factor || continue
        x, y = I[1] + rand(rng, (-r):r), I[2] + rand(rng, (-r):r)
        (1 <= x <= n1 && 1 <= y <= n2 && _bbf_free(tree, x, y, d)) || continue
        (p.site_limit < 0 || _bbf_count(tree, x, y, p.neighbour_radius) <= p.site_limit) ||
            continue
        tree[x, y], u.age[x, y], u.mass[x, y] = true, Int16(0), 0.0
        u.beetles[x, y], u.fate[x, y] = 0.0, BBF_STANDING
    end
    return u
end

# The rows of the strip harvested in `year`, or an empty range.
function _bbf_strip(p, n1, year)
    (p.harvest && p.harvest_start < year <= p.harvest_start + p.harvest_years) || return 1:0
    s = mod(year - p.harvest_start - 1, p.strips)
    width = n1 ÷ p.strips
    return (s * width + 1):(s == p.strips - 1 ? n1 : (s + 1) * width)
end

function _bbf_harvest!(u, p, year)
    for j in axes(u.tree, 2), i in _bbf_strip(p, size(u.tree, 1), year)
        if u.tree[i, j] && u.fate[i, j] == BBF_STANDING && u.age[i, j] > p.harvest_age
            u.fate[i, j], u.beetles[i, j] = BBF_FELLED, 0.0
        end
    end
    return u
end

_bbf_host(u, I) = @inbounds u.tree[I] && u.fate[I] == BBF_STANDING

function _bbf_infest!(u, p, rng)
    rows = max(p.infest_rows[1], 1):min(p.infest_rows[2], size(u.tree, 1))
    cols = max(p.infest_cols[1], 1):min(p.infest_cols[2], size(u.tree, 2))
    hosts = [I for I in CartesianIndices((rows, cols)) if _bbf_host(u, I)]
    isempty(hosts) && return u
    for _ in 1:p.infest_draws
        u.beetles[rand(rng, hosts)] += p.infest_beetles
    end
    return u
end

function _bbf_beetles!(u, p)
    B = u.beetles
    any(>(0), B) || return u
    n1, n2 = size(B)
    host = [_bbf_host(u, I) for I in CartesianIndices(B)]
    S = _bbf_integral!(zeros(Int, n1 + 1, n2 + 1), host)
    r = p.beetle_radius
    added = zeros(n1, n2)
    for _ in 1:p.beetle_spreads
        fill!(added, 0.0)
        for I in findall(>(0), B)
            i, j = Tuple(I)
            B[I] *= p.beetle_growth
            k = _bbf_box(S, i, j, r) - host[I]
            k > 0 || continue
            share = p.whole_beetles ? floor(B[I] / k) : B[I] / k
            share > 0 || continue
            for y in max(j - r, 1):min(j + r, n2), x in max(i - r, 1):min(i + r, n1)
                (host[x, y] && (x, y) != (i, j)) && (added[x, y] += share)
            end
        end
        B .+= added
    end
    for I in eachindex(B)
        if B[I] >= p.beetle_kill
            u.fate[I], B[I] = BBF_KILLED, 0.0
        end
    end
    return u
end

"""
    bark_beetle_step!(u, p, year, rng) -> u

The part of the forest's `year` that depends on the order of events, applied
to the state `u` in place after the rule's step, in this order:

1. With `harvest`, in the years after `harvest_start` for `harvest_years`
   years, the trees older than `harvest_age` in one of `strips` strips of
   rows, in turn, are felled, with their beetles.
2. In year `infestation`, `infest_draws` trees, drawn with replacement from
   those standing (not felled) in the rows `infest_rows` and columns
   `infest_cols`, each gain `infest_beetles` beetles. If no tree stands there,
   none do.
3. `beetle_spreads` times, all at once, every tree with beetles multiplies
   them by `beetle_growth`, and the same number is shared equally among the
   other standing trees in the square of half-width `beetle_radius` around
   it; with `whole_beetles` each share is rounded down.
4. Every tree with `beetle_kill` beetles or more is killed, and its beetles
   with it.
5. In years divisible by `spread_interval`, the trees standing at that
   moment, those felled or killed this year included, are visited in column
   order. A tree with `n` other trees in the square of half-width
   `neighbour_radius` around it spreads with probability
   `spread_probabilities[k]`, for the first class `k` with `n ≤
   spread_limits[k]` (none beyond the last class). Its sapling, of age 0,
   lands in a random cell of the square of half-width `sapling_radius`, and
   takes root only if no tree stands within distance `spacing` of it and at
   most `site_limit` trees stand in the square of half-width
   `neighbour_radius` around it, saplings planted earlier in the same year
   included.

Trees felled or killed are removed at the start of the next year.
"""
function bark_beetle_step!(u, p, year::Integer, rng)
    _bbf_harvest!(u, p, year)
    year == p.infestation && _bbf_infest!(u, p, rng)
    _bbf_beetles!(u, p)
    year % p.spread_interval == 0 && _bbf_spread!(u, p, year, rng)
    return u
end

"""
    bark_beetle_callback() -> DiscreteCallback

A callback that runs [`bark_beetle_step!`](@ref) after each year's rule step,
with a random number generator seeded from the problem's seed and the year.
"""
function bark_beetle_callback()
    affect! = function (integ)
        year = Int(integ.t)
        rng = Random.Xoshiro(hash((integ.prob.seed, year, :bark_beetle)))
        bark_beetle_step!(integ.u, integ.p, year, rng)
        return nothing
    end
    return DiscreteCallback((u, t, integ) -> true, affect!; save_positions=(false, false))
end

"""
    problem_bark_beetle_forest(; n=300, years=300, contamination=false, harvest=false,
                               droughts=(), infestation=nothing, seed=0, pmap=(;),
                               kwargs...) -> CAProblem

The setup of [Konigswieser2020](@cite): an `n × n` forest of 1 m cells (9 ha
for the default), starting from 20 trees aged 10 (and 280 kg) in distinct
cells within 20 m of the centre, placed at random from `seed`, for `years`
years. Each scenario is off unless given:

- `contamination = true` contaminates the rows 90 to 120 and columns 120 to
  180;
- `harvest = true` fells the strips from year 101 to year 250;
- `droughts` lists the years after which a drought of 10 years starts, such
  as `(50, 250)`;
- `infestation` is the year 2000 beetles are brought into the rows and
  columns 100 to 120.

The state at time `t` is the forest at the end of year `t`. Nothing lies
beyond the forest.
"""
function problem_bark_beetle_forest(; n::Integer=300, years::Integer=300,
                                    contamination::Bool=false, harvest::Bool=false,
                                    droughts=(),
                                    infestation::Union{Nothing,Integer}=nothing,
                                    trees::Integer=20, planting_radius::Real=20,
                                    initial_age::Integer=10, seed::Integer=0, pmap=(;),
                                    kwargs...)
    n >= 1 || throw(ArgumentError("n must be positive, got $n"))
    all(d -> d isa Integer, droughts) ||
        throw(ArgumentError("droughts must list years, got $droughts"))
    model = build_bark_beetle_forest()
    given = merge((; harvest, droughts=Tuple(droughts),
                   infestation=something(infestation, -1)),
                  (; (Symbol(k) => v for (k, v) in pairs(pmap))...))
    p = merge(defaults(model), given)
    c = (n + 1) / 2
    near = [I
            for I in CartesianIndices((n, n))
            if (I[1] - c)^2 + (I[2] - c)^2 <= planting_radius^2]
    length(near) >= trees ||
        throw(ArgumentError("$trees trees do not fit within $planting_radius m of the centre"))
    rng = Random.Xoshiro(hash((seed, 0, :bark_beetle)))
    tree = falses(n, n)
    age = zeros(Int16, n, n)
    mass = zeros(n, n)
    for I in near[Random.randperm(rng, length(near))[1:trees]]
        tree[I], age[I] = true, initial_age
        mass[I] = sum(a -> _bbf_growth(a, p), 0:(initial_age - 1); init=0.0)
    end
    contaminated = falses(n, n)
    contamination && (contaminated[90:min(120, n), 120:min(180, n)] .= true)
    u0 = (; tree=collect(tree), age, mass, contaminated=collect(contaminated))
    return CAProblem(model, u0, years, given; space=Grid(n, n), seed,
                     callback=bark_beetle_callback(), kwargs...)
end

"""
    bark_beetle_counts(u) -> NamedTuple

The forest's totals at the end of a year: standing `trees` (not felled or
killed that year), their wood `mass` in tonnes, the `felled` mass in tonnes,
the trees `killed` by beetles that year, and the `infested` trees and the
`beetles` on them.
"""
function bark_beetle_counts(u)
    trees, mass, felled, killed, infested, beetles = 0, 0.0, 0.0, 0, 0, 0.0
    for I in eachindex(u.tree)
        u.tree[I] || continue
        f = u.fate[I]
        if f == BBF_STANDING
            trees += 1
            mass += u.mass[I]
            u.beetles[I] > 0 && (infested += 1; beetles += u.beetles[I])
        elseif f == BBF_FELLED
            felled += u.mass[I]
        else
            killed += 1
        end
    end
    return (; trees, mass=mass / 1000, felled=felled / 1000, killed, infested, beetles)
end

"""
    bark_beetle_growth_rate(mass; from=50, to=200) -> Float64

The source's growth rate of the forest: the mean yearly relative change of
`mass` (indexed by year from 0) over the years `from` to `to`, in percent.
"""
function bark_beetle_growth_rate(mass::AbstractVector; from::Integer=50, to::Integer=200)
    (1 <= from <= to && to + 1 <= length(mass)) ||
        throw(ArgumentError("the years $from to $to must lie within the series"))
    return 100 * sum(t -> mass[t + 1] / mass[t] - 1, from:to) / (to - from + 1)
end

_register!(:bark_beetle_forest; build=build_bark_beetle_forest,
           problem=problem_bark_beetle_forest)
