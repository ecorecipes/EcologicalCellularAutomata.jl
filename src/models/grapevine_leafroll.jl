# Atallah, Gómez, Conrad & Nyrop (2015) and Atallah (2014, Chapter 2): grapevine leafroll disease in a
# vineyard of 130 x 44 vines, in monthly steps, with roguing, replanting and virus testing of
# neighbours. A vine is healthy, exposed-undetectable, exposed-detectable, infective-moderate or
# infective-high, and has an age in months. The disease is the model's rule; a management strategy is a
# step callback stored with the problem, which flags the vines it rogued or tested, so that the net
# present value can be computed afterwards from the saved states.

const GLR_H = Int8(0)
const GLR_EU = Int8(1)
const GLR_ED = Int8(2)
const GLR_IM = Int8(3)
const GLR_IH = Int8(4)

# Triangular(min, mode, max) by inverting its distribution function at u.
function _triangular(u, lo, mode, hi)
    f = (mode - lo) / (hi - lo)
    return u < f ? lo + sqrt(u * (hi - lo) * (mode - lo)) :
           hi - sqrt((1 - u) * (hi - lo) * (hi - mode))
end
_undetectable_months(u, p) = Int16(ceil(_triangular(u, p.eu_min, p.eu_mode, p.eu_max)))
function _triangular_cdf(x, lo, mode, hi)
    x <= lo && return 0.0
    x >= hi && return 1.0
    return x <= mode ? (x - lo)^2 / ((hi - lo) * (mode - lo)) :
           1 - (hi - x)^2 / ((hi - lo) * (hi - mode))
end

"""
    build_grapevine_leafroll(; undetectable=:duration) -> CAModel

The grapevine leafroll model of [Atallah2015](@cite), in monthly steps. Each
cell is a vine with an infection state (`state`: 0 healthy, 1
exposed-undetectable, 2 exposed-detectable, 3 infective-moderate, 4
infective-high), an age in months (`age`), a clock (`clock`: the months left
undetectable, or the months spent in an infective state), and flags set by
the management callback (`rogued`, `tested`).
Columns of the lattice are trellis rows. Each month, synchronously:

- a healthy vine is infected with probability `1 - exp(-(nα α + nβ β))`, where
  `nα` and `nβ` count the infective (3 or 4) neighbours along its trellis row
  and in the neighbouring rows;
- with `undetectable = :duration` a new infection stays undetectable for a
  triangular number of months (`eu_min`, `eu_mode`, `eu_max`, rounded up);
  with `:cdf`, as the source writes it, an undetectable vine becomes detectable
  each month with probability `F(x)`, the triangular distribution function at
  the months `x` it has been undetectable;
- an exposed-detectable vine becomes infective-moderate with probability
  `1 - exp(-1/L)`, where `L` is `L_young`, `L_mature` or `L_old` by age (under
  60 months, under 240, older);
- an infective-moderate vine becomes infective-high with probability
  `1 - exp(-1/inf_period)`;
- every vine ages by a month.

Management is [`grapevine_strategy_callback`](@ref).
"""
function build_grapevine_leafroll(; undetectable::Symbol=:duration)
    undetectable in (:duration, :cdf) ||
        throw(ArgumentError("undetectable must be :duration or :cdf, got :$undetectable"))
    duration = undetectable === :duration
    # along the trellis row (first axis), then across rows
    nbhd = Offsets([CartesianIndex(-1, 0), CartesianIndex(1, 0), CartesianIndex(0, -1),
                    CartesianIndex(0, 1)])
    rule = LocalRule(nbhd) do cell, nbrs, ctx
        p = ctx.p
        s = cell.state
        age = cell.age + Int16(1)
        clock = cell.clock
        if s == GLR_H
            na = (nbrs[1].state >= GLR_IM) + (nbrs[2].state >= GLR_IM)
            nb = (nbrs[3].state >= GLR_IM) + (nbrs[4].state >= GLR_IM)
            if na + nb > 0 && bernoulli(ctx, -expm1(-(na * p.α + nb * p.β)), 1)
                s = GLR_EU
                clock = duration ? _undetectable_months(uniform(ctx, 2), p) : Int16(0)
            end
        elseif s == GLR_EU
            if duration
                clock -= Int16(1)
                clock <= 0 && ((s, clock) = (GLR_ED, Int16(0)))
            else
                F = _triangular_cdf(clock, p.eu_min, p.eu_mode, p.eu_max)
                clock += Int16(1)
                bernoulli(ctx, F, 2) && ((s, clock) = (GLR_ED, Int16(0)))
            end
        elseif s == GLR_ED
            L = cell.age < 60 ? p.L_young : cell.age < 240 ? p.L_mature : p.L_old
            bernoulli(ctx, -expm1(-1 / L), 3) && ((s, clock) = (GLR_IM, Int16(-1)))
        elseif s == GLR_IM
            bernoulli(ctx, -expm1(-1 / p.inf_period), 4) && (s = GLR_IH)
        end
        s >= GLR_IM && (clock += Int16(1))
        return (state=s, age, clock, rogued=false, tested=false)
    end
    parameters = [Parameter(:α, 4.2;
                            description="infection rate per infective neighbour along the trellis row (1/month)"),
                  Parameter(:β, 0.014;
                            description="infection rate per infective neighbour across rows (1/month)"),
                  Parameter(:eu_min, 4.0;
                            description="shortest undetectable period (months)"),
                  Parameter(:eu_mode, 12.0;
                            description="commonest undetectable period (months)"),
                  Parameter(:eu_max, 18.0;
                            description="longest undetectable period (months)"),
                  Parameter(:L_young, 24.0;
                            description="mean latency of young vines (months)"),
                  Parameter(:L_mature, 48.0;
                            description="mean latency of mature vines (months)"),
                  Parameter(:L_old, 72.0; description="mean latency of old vines (months)"),
                  Parameter(:inf_period, 36.0;
                            description="mean time infective-moderate (months)")]
    return CAModel(rule; name=:grapevine_leafroll,
                   states=[StateVariable(:state, GLR_H; domain=Int8(0):Int8(4),
                                         description="0 H, 1 exposed-undetectable, 2 exposed-detectable, 3 infective-moderate, 4 infective-high"),
                           StateVariable(:age, Int16(0); description="age (months)"),
                           StateVariable(:clock, Int16(0);
                                         description="months left undetectable, or months infective"),
                           StateVariable(:rogued, false;
                                         description="rogued and replanted this month"),
                           StateVariable(:tested, false;
                                         description="virus-tested this month")],
                   parameters,
                   description="Atallah et al. (2015) grapevine leafroll, undetectable period as $undetectable")
end

const GLR_STRATEGIES = (:none, :I, :Im, :Ih, :ImY, :ImM, :ImO, :IhM, :IhO, :ImNS, :ImNSEW,
                        :ImNS2EW, :ImNS2EW2, :ImY_NS, :ImM_NS, :ImO_NS, :ImY_NSEW,
                        :ImM_NSEW, :ImO_NSEW)

_age_class(a) = a < 60 ? :Y : a < 240 ? :M : :O

# The trigger of a strategy (which vines are rogued) and the neighbours its spatial form tests.
function _strategy(strategy::Symbol)
    strategy in GLR_STRATEGIES ||
        throw(ArgumentError("strategy must be one of $GLR_STRATEGIES, got :$strategy"))
    name = String(strategy)
    test = if endswith(name, "NS2EW2")
        ((-1, 0), (1, 0), (-2, 0), (2, 0), (0, -1), (0, 1), (0, -2), (0, 2))
    elseif endswith(name, "NS2EW")
        ((-1, 0), (1, 0), (-2, 0), (2, 0), (0, -1), (0, 1))
    elseif endswith(name, "NSEW")
        ((-1, 0), (1, 0), (0, -1), (0, 1))
    elseif endswith(name, "NS")
        ((-1, 0), (1, 0))
    else
        ()
    end
    base = replace(split(name, "_")[1], r"NS.*$" => "")
    states = base == "none" ? () :
             startswith(base, "Ih") ? (GLR_IH,) :
             base == "I" ? (GLR_IM, GLR_IH) : (GLR_IM,)
    ages = length(base) > 2 && base[end] in ('Y', 'M', 'O') ? (Symbol(base[end]),) :
           (:Y, :M, :O)
    return states, ages, test
end

"""
    grapevine_manage!(u, strategy; delay=0, maturity=36, protect_replants=true)
        -> (rogued, tested)

Apply a management strategy of [Atallah2015](@cite) to the vineyard state `u`
in place, and flag the vines rogued and tested. Every vine meeting the trigger
(`:I`, infective; `:Im`, infective-moderate; `:Ih`, infective-high; with an
age class `Y`, `M` or `O` appended, only vines of that class) is rogued and
replaced by a healthy vine of age 0. A spatial strategy (`:ImNS`, `:ImNSEW`,
`:ImNS2EW`, `:ImNS2EW2`, and their age-class forms such as `:ImM_NSEW`) also
tests the listed neighbours of each triggering vine that are not infective,
once each, and rogues those found exposed-detectable. With `protect_replants`,
vines younger than `maturity` months are never rogued. A vine triggers only
once it has been infective for `delay` months, so that with `delay = 1` it
infects its neighbours for a month before it is found. Returns the numbers of
vines rogued and tested.
"""
function grapevine_manage!(u, strategy::Symbol; delay::Integer=0, maturity::Integer=36,
                           protect_replants::Bool=true)
    states, ages, test = _strategy(strategy)
    isempty(states) && return (0, 0)
    s, age = u.state, u.age
    can_rogue(I) = !protect_replants || age[I] >= maturity
    triggers = [I
                for I in CartesianIndices(s)
                if s[I] in states && u.clock[I] >= delay && _age_class(age[I]) in ages &&
                       can_rogue(I)]
    totest = Set{CartesianIndex{2}}()
    for I in triggers, (di, dj) in test
        J = I + CartesianIndex(di, dj)
        checkbounds(Bool, s, J) && s[J] <= GLR_ED && push!(totest, J)
    end
    positive = [J for J in totest if s[J] == GLR_ED && can_rogue(J)]
    for J in totest
        u.tested[J] = true
    end
    for I in Iterators.flatten((triggers, positive))
        s[I], age[I], u.clock[I], u.rogued[I] = GLR_H, Int16(0), Int16(0), true
    end
    return (length(triggers) + length(positive), length(totest))
end

"""
    grapevine_strategy_callback(strategy; delay=0, maturity=36, protect_replants=true)
        -> DiscreteCallback

A callback that applies [`grapevine_manage!`](@ref) after every monthly step.
"""
function grapevine_strategy_callback(strategy::Symbol; delay::Integer=0,
                                     maturity::Integer=36, protect_replants::Bool=true)
    _strategy(strategy)
    affect! = integ -> (grapevine_manage!(integ.u, strategy; delay, maturity,
                                          protect_replants);
                        nothing)
    return DiscreteCallback((u, t, integ) -> true, affect!; save_positions=(false, false))
end

"""
    problem_grapevine_leafroll(; strategy=:none, undetectable=:duration, delay=1,
                               initial=:Eu, rows=130, cols=44, infected=0.02, months=600,
                               protect_replants=true, seed=0, pmap=(;), kwargs...)
        -> CAProblem

The setup of [Atallah2015](@cite): a vineyard of `cols` trellis rows of `rows`
vines, all healthy and newly planted, with a share `infected` of them, chosen
at random from `seed`, infected at the start, for 600 months under a
management `strategy` (see [`grapevine_manage!`](@ref); `:none` for no
control) whose detection takes `delay` months. The initial infections are
exposed-undetectable (`initial = :Eu`, as the source says "exposed"),
exposed-detectable (`:Ed`) or infective-moderate (`:Im`). Nothing lies beyond
the vineyard.
"""
function problem_grapevine_leafroll(; strategy::Symbol=:none,
                                    undetectable::Symbol=:duration,
                                    delay::Integer=1, initial::Symbol=:Eu,
                                    rows::Integer=130, cols::Integer=44,
                                    infected::Real=0.02, months::Integer=600,
                                    protect_replants::Bool=true, seed::Integer=0,
                                    pmap=(;), kwargs...)
    0 <= infected <= 1 || throw(ArgumentError("infected must lie in [0, 1], got $infected"))
    start = (Eu=GLR_EU, Ed=GLR_ED, Im=GLR_IM)
    haskey(start, initial) ||
        throw(ArgumentError("initial must be :Eu, :Ed or :Im, got :$initial"))
    model = build_grapevine_leafroll(; undetectable)
    p = merge(defaults(model), (; (Symbol(k) => v for (k, v) in pairs(pmap))...))
    rng = Random.Xoshiro(seed)
    state = fill(GLR_H, rows, cols)
    clock = zeros(Int16, rows, cols)
    for I in Random.randperm(rng, rows * cols)[1:round(Int, infected * rows * cols)]
        state[I] = start[initial]
        draw = rand(rng)
        clock[I] = initial === :Eu && undetectable === :duration ?
                   _undetectable_months(draw, p) : 0
    end
    callback = grapevine_strategy_callback(strategy; delay, protect_replants)
    return CAProblem(model, (; state, clock), months, pmap; space=Grid(rows, cols), seed,
                     callback, kwargs...)
end

"""
    grapevine_revenue(u; price=0.4266, maturity=36, yield_loss=(0.3, 0.5, 0.75),
                      quality_loss=0.1) -> Float64

The vineyard's revenue for a month, in dollars: each vine at least `maturity`
months old earns `price` when healthy, reduced by the yield loss of its
state (exposed, infective-moderate, infective-high) and by `quality_loss` when
infected.
"""
function grapevine_revenue(u; price::Real=0.4266, maturity::Integer=36,
                           yield_loss=(0.3, 0.5, 0.75), quality_loss::Real=0.1)
    r = 0.0
    for I in eachindex(u.state)
        u.age[I] >= maturity || continue
        s = u.state[I]
        r += s == GLR_H ? price :
             price * (1 - yield_loss[s <= GLR_ED ? 1 : s - 1]) * (1 - quality_loss)
    end
    return r
end

"""
    grapevine_npv(sol; discount=0.9959, rogue_cost=7.25, test_cost=2.61, kwargs...) -> Float64

The net present value of a run saved every month: the revenue of each month's
state ([`grapevine_revenue`](@ref), with `kwargs`) less the cost of the vines
rogued and tested at the end of that month, discounted by `discount` a month.
"""
function grapevine_npv(sol; discount::Real=0.9959, rogue_cost::Real=7.25,
                       test_cost::Real=2.61, kwargs...)
    all(diff(sol.t) .== 1) || throw(ArgumentError("the run must be saved every month"))
    v = 0.0
    for k in 1:(length(sol.u) - 1)
        next = sol.u[k + 1]
        cost = rogue_cost * count(next.rogued) + test_cost * count(next.tested)
        v += discount^sol.t[k] * (grapevine_revenue(sol.u[k]; kwargs...) - cost)
    end
    return v
end

"""
    grapevine_half_life(sol) -> Union{Int,Nothing}

The first saved month at which half the vines or more are infected (not
healthy), or `nothing` if that never happens in the run.
"""
function grapevine_half_life(sol)
    for (t, u) in zip(sol.t, sol.u)
        2count(!=(GLR_H), u.state) >= length(u.state) && return t
    end
    return nothing
end

_register!(:grapevine_leafroll; build=build_grapevine_leafroll,
           problem=problem_grapevine_leafroll)
