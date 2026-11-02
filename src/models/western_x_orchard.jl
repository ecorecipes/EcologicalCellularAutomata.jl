# Ghorbani, Atallah & Gallardo (2026): Western X disease, a phytoplasma spread by leafhoppers, in a cherry
# orchard of 20 x 30 trees, in monthly steps over 25 years. A tree is unvisited, visited (inoculated),
# infected, visible (symptomatic) or high (late, symptomatic), with an age, or the cell is empty after a
# tree is cut. The disease is the model's rule; removal, replanting and testing are a callback stored with
# the problem, which flags what it did each month so that yields, costs and net present values can be
# computed afterwards from the saved states. The source's decision model (a partially observable Markov
# decision process solved by rollout) is not implemented: only its fixed strategies are.

const WXD_EMPTY = Int8(0)
const WXD_UNVISITED = Int8(1)
const WXD_VISITED = Int8(2)
const WXD_INFECTED = Int8(3)
const WXD_VISIBLE = Int8(4)
const WXD_HIGH = Int8(5)
const WXD_STRATEGIES = (:none, :threshold_cut, :threshold_replant, :test_cut, :test_replant)

# The monthly probability that an unvisited tree with `k` symptomatic neighbours is inoculated.
function _wxd_infection(k, p)
    q = min(p.alpha + p.beta * k, 1.0)
    return p.rates === :annual ? -expm1(log1p(-q) / 12) : q
end

"""
    build_western_x_orchard() -> CAModel

The Western X disease model of [Ghorbani2026](@cite), in monthly steps. Each
cell is a cherry tree with a disease `state` (0 an empty cell where a tree was
cut, 1 unvisited, 2 visited by an infective leafhopper, 3 infected, 4 visible,
5 high), an `age` in months, a `clock` of months in its state, and flags set by
the management callback (`cut`, `replanted`, `tested`). Each month, in every
cell at once:

- an unvisited tree is inoculated with probability `alpha + beta k`, where `k`
  counts the visible and high trees among its eight neighbours; with
  `rates = :annual` this is a yearly probability, converted to a monthly one;
- a visited, infected or visible tree that has spent a whole number of years
  in its state moves to the next with probability `progression`;
- every tree ages a month.

Management is [`western_x_manage!`](@ref); the economics are
[`western_x_economics`](@ref).
"""
function build_western_x_orchard()
    rule = LocalRule(Moore(1)) do c, nbrs, ctx
        p = ctx.p
        s = c.state
        s == WXD_EMPTY &&
            return (state=s, age=Int16(0), clock=Int16(0), cut=false, replanted=false,
                    tested=false)
        clock = c.clock + Int16(1)
        if s == WXD_UNVISITED
            k = count(n -> n.state >= WXD_VISIBLE, nbrs)
            if bernoulli(ctx, _wxd_infection(k, p), 1)
                s, clock = WXD_VISITED, Int16(0)
            end
        elseif s < WXD_HIGH && clock % 12 == 0 && bernoulli(ctx, p.progression, 2)
            s, clock = s + Int8(1), Int16(0)
        end
        return (state=s, age=c.age + Int16(1), clock, cut=false, replanted=false,
                tested=false)
    end
    states = [StateVariable(:state, WXD_UNVISITED; domain=Int8(0):Int8(5),
                            description="0 empty, 1 unvisited, 2 visited, 3 infected, 4 visible, 5 high"),
              StateVariable(:age, Int16(0); description="age of the tree (months)"),
              StateVariable(:clock, Int16(0); description="months in the present state"),
              StateVariable(:cut, false; description="cut this month"),
              StateVariable(:replanted, false; description="replanted this month"),
              StateVariable(:tested, false; description="tested this month")]
    parameters = [Parameter(:alpha, 0.015;
                            description="baseline probability of inoculation (a year with rates = :annual)"),
                  Parameter(:beta, 0.2;
                            description="increase of that probability for each symptomatic neighbour"),
                  Parameter(:rates, :annual;
                            description=":annual (alpha + beta k a year) or :monthly"),
                  Parameter(:progression, 0.5;
                            description="probability of moving to the next state at each full year in a state"),
                  Parameter(:strategy, :none;
                            description="management: :none, :threshold_cut, :threshold_replant, :test_cut or :test_replant"),
                  Parameter(:threshold, 0.01;
                            description="share of trees symptomatic at which threshold strategies cut the symptomatic trees"),
                  Parameter(:interval, 12;
                            description="months between management decisions"),
                  Parameter(:first_decision, 12; description="month of the first decision"),
                  Parameter(:radius, 1;
                            description="Moore radius around symptomatic trees within which trees are tested"),
                  Parameter(:test_share, 0.1;
                            description="share of the trees within the radius that are tested"),
                  Parameter(:sensitivity, 0.85;
                            description="probability a test detects a visited or infected tree")]
    return CAModel(rule; name=:western_x_orchard, states, parameters,
                   description="Ghorbani, Atallah & Gallardo (2026) Western X disease in a cherry orchard")
end

"""
    western_x_manage!(u, p, month, rng) -> u

The fixed management strategies of [Ghorbani2026](@cite), applied to the
state `u` in place after the disease's step, in decision months (`month` at
or after `first_decision`, every `interval` months). Visual scouting finds
every visible and high tree.

- `:threshold_cut` and `:threshold_replant` cut every symptomatic tree once at
  least `threshold` of the cells hold one (and replant it, as an unvisited tree
  of age 0).
- `:test_cut` and `:test_replant` also test a share `test_share` of the trees
  within Moore radius `radius` of each symptomatic tree that are not
  symptomatic themselves; a test detects a visited or infected tree with
  probability `sensitivity` and never errs on an unvisited one. Symptomatic
  and detected trees are cut (and replanted).

The trees cut, replanted and tested are flagged for the month.
"""
function western_x_manage!(u, p, month::Integer, rng)
    p.strategy === :none && return u
    (month >= p.first_decision && (month - p.first_decision) % p.interval == 0) || return u
    S = u.state
    symptomatic = findall(>=(WXD_VISIBLE), S)
    isempty(symptomatic) && return u
    targets = Set(symptomatic)
    if p.strategy in (:threshold_cut, :threshold_replant)
        length(symptomatic) >= p.threshold * length(S) || return u
    else
        candidates = Set{CartesianIndex{2}}()
        r = p.radius
        for I in symptomatic, J in CartesianIndices(S)
            maximum(abs.(Tuple(J - I))) <= r || continue
            WXD_UNVISITED <= S[J] <= WXD_INFECTED && push!(candidates, J)
        end
        for J in sort!(collect(candidates))
            rand(rng) < p.test_share || continue
            u.tested[J] = true
            S[J] in (WXD_VISITED, WXD_INFECTED) && rand(rng) < p.sensitivity &&
                push!(targets, J)
        end
    end
    replant = p.strategy in (:threshold_replant, :test_replant)
    for I in targets
        u.cut[I] = true
        u.replanted[I] = replant
        S[I] = replant ? WXD_UNVISITED : WXD_EMPTY
        u.age[I], u.clock[I] = 0, 0
    end
    return u
end

"""
    western_x_callback() -> DiscreteCallback

A callback that runs [`western_x_manage!`](@ref) after each month's disease
step, with a random number generator seeded from the problem's seed and the
month.
"""
function western_x_callback()
    affect! = function (integ)
        month = Int(integ.t)
        rng = Random.Xoshiro(hash((integ.prob.seed, month, :western_x)))
        western_x_manage!(integ.u, integ.p, month, rng)
        return nothing
    end
    return DiscreteCallback((u, t, integ) -> true, affect!; save_positions=(false, false))
end

"""
    problem_western_x_orchard(; strategy=:none, rows=20, cols=30, initial=(0.02, 0.08),
                              disease=true, years=25, seed=0, pmap=(;), kwargs...)
        -> CAProblem

The setup of [Ghorbani2026](@cite): a newly planted orchard of `rows × cols`
trees (an acre), of which a share drawn uniformly from `initial` (a range, or
a number), chosen at random from `seed`, start visited, for `years` years in
monthly steps, under a management `strategy` (see
[`western_x_manage!`](@ref); `:none` for none). With `disease = false`, no
tree is ever inoculated (the source's disease-free scenario). Nothing lies
beyond the orchard. Management decides every other year, as in the source:
from month 24 for threshold strategies and from month 48 for testing
strategies, unless `pmap` says otherwise.
"""
function problem_western_x_orchard(; strategy::Symbol=:none, rows::Integer=20,
                                   cols::Integer=30, initial=(0.02, 0.08),
                                   disease::Bool=true, years::Integer=25,
                                   seed::Integer=0, pmap=(;), kwargs...)
    strategy in WXD_STRATEGIES ||
        throw(ArgumentError("strategy must be one of $WXD_STRATEGIES, got :$strategy"))
    lo, hi = initial isa Real ? (initial, initial) : initial
    0 <= lo <= hi <= 1 || throw(ArgumentError("initial must lie in [0, 1], got $initial"))
    testing = strategy in (:test_cut, :test_replant)
    schedule = (; interval=24, first_decision=testing ? 48 : 24)
    given = merge((; strategy), schedule,
                  disease ? (;) : (; alpha=0.0, beta=0.0),
                  (; (Symbol(k) => v for (k, v) in pairs(pmap))...))
    model = build_western_x_orchard()
    p = merge(defaults(model), given)
    p.rates in (:annual, :monthly) ||
        throw(ArgumentError("rates must be :annual or :monthly, got :$(p.rates)"))
    rng = Random.Xoshiro(hash((seed, 0, :western_x)))
    share = disease ? lo + (hi - lo) * rand(rng) : 0.0
    state = fill(WXD_UNVISITED, rows, cols)
    n = round(Int, share * rows * cols)
    state[Random.randperm(rng, rows * cols)[1:n]] .= WXD_VISITED
    return CAProblem(model, (; state), 12years, given; space=Grid(rows, cols), seed,
                     callback=western_x_callback(), kwargs...)
end

# The source's yield (lb a month) and cost ($ a month) of a healthy tree of `age` months.
_wxd_yield(age) = age < 24 ? 0.0 : age < 36 ? 0.9 : age < 48 ? 1.8 : age < 60 ? 2.2 : 2.7
function _wxd_cost(age)
    return age < 12 ? 2.0 :
           age < 24 ? 0.6 : age < 36 ? 1.9 : age < 48 ? 2.9 :
                                             age < 60 ? 3.4 : 3.9
end
const WXD_LOSS = (0.0, 0.1, 0.25, 0.5, 0.75)

"""
    western_x_economics(sol; price=2.35, discount=0.05, no_harvest=0.25, cut_cost=33,
                        plant_cost=14.3, test_cost=300, scale_costs=false) -> NamedTuple

The economics of a run saved every month, per orchard, as in
[Ghorbani2026](@cite): each tree yields `y(age)` lb a month and costs `c(age)`
dollars a month; the yield is multiplied by `1 - loss`, with losses 0.1, 0.25,
0.5 and 0.75 in the visited, infected, visible and high states, and so is the
cost with `scale_costs = true`, as the source's equation writes it; an empty
cell yields and costs nothing. A month in
which at least `no_harvest` of the cells hold symptomatic trees brings no
revenue (the source's growers do not harvest a badly diseased block). Each
tree cut costs `cut_cost`, each replanted `plant_cost` more, each tested
`test_cost`. Values are discounted at `discount` a year. Returns the yearly
`yield` (lb) and net revenue (`net`, undiscounted), and the net present value
`npv`.
"""
function western_x_economics(sol; price::Real=2.35, discount::Real=0.05,
                             no_harvest::Real=0.25, cut_cost::Real=33,
                             plant_cost::Real=14.3, test_cost::Real=300,
                             scale_costs::Bool=false)
    all(diff(sol.t) .== 1) || throw(ArgumentError("the run must be saved every month"))
    months = length(sol.u) - 1
    δ = (1 + discount)^(-1 / 12)
    yield = zeros(cld(months, 12))
    net = zeros(cld(months, 12))
    npv = 0.0
    for k in 1:months
        u = sol.u[k + 1]
        S = u.state
        harvest = count(>=(WXD_VISIBLE), S) < no_harvest * length(S)
        y, c = 0.0, 0.0
        for I in eachindex(S)
            S[I] == WXD_EMPTY && continue
            f = 1 - WXD_LOSS[S[I]]
            y += f * _wxd_yield(u.age[I])
            c += (scale_costs ? f : 1.0) * _wxd_cost(u.age[I])
        end
        y = harvest ? y : 0.0
        mgmt = cut_cost * count(u.cut) + plant_cost * count(u.replanted) +
               test_cost * count(u.tested)
        v = price * y - c - mgmt
        year = cld(k, 12)
        yield[year] += y
        net[year] += v
        npv += δ^k * v
    end
    return (; yield, net, npv)
end

_register!(:western_x_orchard; build=build_western_x_orchard,
           problem=problem_western_x_orchard)
