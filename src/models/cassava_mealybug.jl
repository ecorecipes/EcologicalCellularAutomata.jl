# Promrak & Rattanakul (2015, 2017): the cassava mealybug and its biological control by green lacewing
# larvae in a cassava field, one plant per cell, in daily steps. Within each plant, the life stages of
# mealybugs (eggs, instars, adults) and lacewings (eggs, larvae, pupae, adults) follow difference
# equations, with lacewing larvae eating mealybugs; that is the model's rule. Between plants, a sweep
# in random order moves instars from infested plants to their neighbours, wind brings instars from
# outside, plants are surveyed and removed, and lacewings are released: that is a step callback stored
# with the problem, as instars move between cells and the release decision needs the whole field.

const CMB_S = Int8(0)
const CMB_I = Int8(1)
const CMB_E = Int8(2)

"""
    build_cassava_mealybug() -> CAModel

The within-plant dynamics of the cassava mealybug model of
[Promrak2015](@cite), with the algorithm of [Promrak2017](@cite). Each cell is
a cassava plant, susceptible (0), infested (1) or removed (2), with mealybug
eggs, instars and adults (`Pe`, `Pi`, `Pm`), lacewing eggs, larvae, pupae and
adults (`Me`, `Mi`, `Md`, `Mm`), the day it was infested (`since`, -1 if never)
and the worst yield loss it has suffered (`loss`: 0 none, 1 a tenth, 2 three
tenths, 3 total). Each day, in every plant:

    Pi' = Pi + r1 α1 Pe - α2 Pi - min(20 Mi, Pi)
    Pm' = Pm + r2 α2 Pi - α3 Pm - min(20 Mi, Pm)
    Pe' = Pe + r3 α4 ν1 Pm - α1 Pe - min(20 Mi, Pe)
    Mi' = Mi + s1 γ1 Me - γ2 Mi
    Md' = Md + s2 γ2 δ1 Mi - γ3 Md,  δ1 = 0.0521 min(1, (Pi + Pm + Pe) / (60 Mi))
    Mm' = Mm + s3 γ3 Md - δ2 Mm
    Me' = Me + s4 ν2 Mm - γ1 Me

(a larva eats at most 20 of each mealybug stage a day). A plant infested for
two weeks or more suffers the loss of its period: total in days 1 to 120,
three tenths in days 121 to 210, a tenth afterwards. Movement between plants,
surveys, removal and releases are [`cassava_management_callback`](@ref).
"""
function build_cassava_mealybug()
    rule = LocalRule(Moore(1)) do c, nbrs, ctx
        p = ctx.p
        c.state == CMB_E && return c
        eat(P) = min(20 * c.Mi, max(P, 0.0))
        # as in the source, counts are not clamped: predation can take a stage below zero,
        # and a plant whose total is not positive is cleared by the sweep
        Pi = c.Pi + p.r1 * p.α1 * c.Pe - p.α2 * c.Pi - eat(c.Pi)
        Pm = c.Pm + p.r2 * p.α2 * c.Pi - p.α3 * c.Pm - eat(c.Pm)
        Pe = c.Pe + p.r3 * p.α4 * p.ν1 * c.Pm - p.α1 * c.Pe - eat(c.Pe)
        prey = c.Pi + c.Pm + c.Pe
        δ1 = c.Mi > 0 ? 0.0521 * min(1.0, prey / (60 * c.Mi)) : 0.0
        Mi = c.Mi + p.s1 * p.γ1 * c.Me - p.γ2 * c.Mi
        Md = c.Md + p.s2 * p.γ2 * δ1 * c.Mi - p.γ3 * c.Md
        Mm = c.Mm + p.s3 * p.γ3 * c.Md - p.δ2 * c.Mm
        Me = c.Me + p.s4 * p.ν2 * c.Mm - p.γ1 * c.Me
        day = ctx.t + 1
        loss = c.loss
        if c.state == CMB_I && c.since >= 0 && day - c.since >= 14
            loss = max(loss, day <= 120 ? Int8(3) : day <= 210 ? Int8(2) : Int8(1))
        end
        return (state=c.state, Pe, Pi, Pm, Me, Mi, Md, Mm, since=c.since, loss)
    end
    states = [StateVariable(:state, CMB_S; domain=Int8(0):Int8(2),
                            description="0 susceptible, 1 infested, 2 removed"),
              StateVariable(:Pe, 0.0; description="mealybug eggs"),
              StateVariable(:Pi, 0.0; description="mealybug instars"),
              StateVariable(:Pm, 0.0; description="adult mealybugs"),
              StateVariable(:Me, 0.0; description="lacewing eggs"),
              StateVariable(:Mi, 0.0; description="lacewing larvae"),
              StateVariable(:Md, 0.0; description="lacewing pupae"),
              StateVariable(:Mm, 0.0; description="adult lacewings"),
              StateVariable(:since, Int16(-1); description="day infested (-1 never)"),
              StateVariable(:loss, Int8(0);
                            description="worst yield loss: 0 none, 1 10%, 2 30%, 3 100%")]
    parameters = [Parameter(:α1, 0.1493; description="mealybug eggs hatching a day"),
                  Parameter(:α2, 0.0388; description="instars maturing a day"),
                  Parameter(:α3, 0.0596; description="adult mealybug death rate"),
                  Parameter(:α4, 0.4468; description="adult females reproducing"),
                  Parameter(:r1, 0.912; description="egg survival"),
                  Parameter(:r2, 0.6589; description="instar survival"),
                  Parameter(:r3, 0.85; description="share of females"),
                  Parameter(:ν1, 2.9126; description="eggs per female mealybug a day"),
                  Parameter(:γ1, 0.2703; description="lacewing eggs hatching a day"),
                  Parameter(:γ2, 0.0625; description="larvae pupating a day"),
                  Parameter(:γ3, 0.1053; description="pupae emerging a day"),
                  Parameter(:δ2, 0.0206; description="adult lacewing death rate"),
                  Parameter(:s1, 0.817; description="lacewing egg survival"),
                  Parameter(:s2, 0.7938; description="larva survival"),
                  Parameter(:s3, 0.7402; description="pupa survival"),
                  Parameter(:s4, 0.485; description="share of female lacewings"),
                  Parameter(:ν2, 2.3876; description="eggs per female lacewing a day"),
                  Parameter(:n1, 0.05;
                            description="instar transfer from a neighbour, ring 1"),
                  Parameter(:n2, 0.005; description="instar transfer, ring 2"),
                  Parameter(:n3, 0.0005; description="instar transfer, ring 3"),
                  Parameter(:w1, 0.01;
                            description="wind infestation on the two border rows"),
                  Parameter(:w2, 0.001; description="wind infestation inside"),
                  Parameter(:m1, Inf;
                            description="mealybugs above which a plant surveyed in days 121-210 is removed"),
                  Parameter(:survey_plants, 60;
                            description="plants surveyed every two weeks"),
                  Parameter(:interval, 60;
                            description="days between releases while mealybugs are found"),
                  Parameter(:G1, 1000.0;
                            description="larvae released when over half the surveyed plants are infested"),
                  Parameter(:G2, 800.0; description="larvae released otherwise"),
                  # the field's management memory, updated by the callback
                  Parameter(:last_release, 0;
                            description="day of the last release (0: none pending)"),
                  Parameter(:found, 0;
                            description="infested plants found at the last survey"),
                  Parameter(:sampled, 0; description="plants sampled at the last survey"),
                  Parameter(:released, 0.0; description="larvae released so far")]
    return CAModel(rule; name=:cassava_mealybug, states, parameters,
                   description="Promrak & Rattanakul (2015) cassava mealybugs and green lacewings")
end

# The three rings of the 7 x 7 neighbourhood, with their ring numbers.
const CMB_RINGS = [(CartesianIndex(i, j), max(abs(i), abs(j))) for j in -3:3
                   for i in -3:3
                   if (i, j) != (0, 0)]

"""
    cassava_sweep!(u, p, day, rng) -> (found, sampled)

One day of movement, survey and removal in the cassava mealybug model of
[Promrak2015](@cite), in place on the state `u`, following the flowcharts of
[Promrak2017](@cite). Plants are visited once each, in random order:

- a susceptible plant is infested by wind (probability `w1` on the two border
  rows, `w2` inside) and, for each infested neighbour not yet visited that has
  an instar to give, gains it with probability `n1`, `n2` or `n3` by ring,
  becoming infested;
- an infested plant whose mealybug total is not positive becomes
  susceptible, its counts set to zero; on a survey day
  (every 14 days from day 31) it is sampled with probability `survey_plants`
  over the plants not removed, and a sampled plant is removed, unless it is
  sampled in days 121 to 210 with at most `m1` mealybugs; otherwise it gains
  instars by wind and from neighbours as above.

Returns the infested plants found and the plants sampled on a survey day, and
`(0, 0)` on other days.
"""
function cassava_sweep!(u, p, day::Integer, rng)
    L1, L2 = size(u.state)
    border(I) = I[1] <= 2 || I[2] <= 2 || I[1] >= L1 - 1 || I[2] >= L2 - 1
    survey = day >= 31 && (day - 31) % 14 == 0
    alive = count(!=(CMB_E), u.state)
    f = survey && alive > 0 ? min(1.0, p.survey_plants / alive) : 0.0
    visited = falses(L1, L2)
    found = 0
    sampled = survey ? min(alive, Int(p.survey_plants)) : 0
    ringp = (p.n1, p.n2, p.n3)
    function gather!(I)
        gained = false
        for (o, ring) in CMB_RINGS
            J = I + o
            checkbounds(Bool, u.state, J) && u.state[J] == CMB_I && !visited[J] &&
            u.Pi[J] >= 1 || continue
            if rand(rng) <= ringp[ring]
                u.Pi[I] += 1
                u.Pi[J] = max(0.0, u.Pi[J] - 1)
                gained = true
            end
        end
        return gained
    end
    for k in Random.randperm(rng, L1 * L2)
        I = CartesianIndices(u.state)[k]
        s = u.state[I]
        if s == CMB_S
            wind = rand(rng) <= (border(I) ? p.w1 : p.w2)
            wind && (u.Pi[I] += 1)
            if gather!(I) || wind
                u.state[I] = CMB_I
                u.since[I] = day
            end
        elseif s == CMB_I
            total = u.Pi[I] + u.Pm[I] + u.Pe[I]
            if total <= 0          # no mealybugs left: cleared
                u.state[I] = CMB_S
                u.since[I] = -1
                u.Pe[I] = u.Pi[I] = u.Pm[I] = 0.0
            else
                if survey && rand(rng) <= f
                    found += 1
                    keep = 120 < day <= 210 && total <= p.m1
                    if !keep
                        u.state[I] = CMB_E
                        u.Pe[I] = u.Pi[I] = u.Pm[I] = 0.0
                        u.Me[I] = u.Mi[I] = u.Md[I] = u.Mm[I] = 0.0
                        visited[I] = true
                        continue
                    end
                end
                rand(rng) <= (border(I) ? p.w1 : p.w2) && (u.Pi[I] += 1)
                gather!(I)
            end
        end
        visited[I] = true
    end
    return found, sampled
end

"""
    cassava_release!(u, p, day) -> (p, larvae)

The release rule of [Promrak2015](@cite): while the last survey found
infested plants, lacewing larvae are released every `interval` days from the
first such survey, `G1` if over half the plants sampled were infested and `G2`
otherwise, split evenly over the infested plants. Returns the updated
management memory and the larvae released.
"""
function cassava_release!(u, p, day::Integer)
    p.found > 0 || return merge(p, (; last_release=0)), 0.0
    due = p.last_release == 0 || (day - p.last_release) % p.interval == 0
    due || return p, 0.0
    infested = findall(==(CMB_I), u.state)
    isempty(infested) && return p, 0.0
    G = p.found > p.sampled / 2 ? p.G1 : p.G2
    for I in infested
        u.Mi[I] += G / length(infested)
    end
    return merge(p, (; last_release=day, released=p.released + G)), G
end

"""
    cassava_management_callback() -> DiscreteCallback

A callback that, after each day's within-plant step, runs
[`cassava_sweep!`](@ref) and [`cassava_release!`](@ref) on the state, keeping
the field's management memory (last release, last survey, larvae released) in
the parameters, with a random number generator seeded from the problem's seed
and the day.
"""
function cassava_management_callback()
    affect! = function (integ)
        day = Int(integ.t)
        rng = Random.Xoshiro(hash((integ.prob.seed, day, :cassava)))
        found, sampled = cassava_sweep!(integ.u, integ.p, day, rng)
        p = sampled > 0 ? merge(integ.p, (; found, sampled)) : integ.p
        p, _ = cassava_release!(integ.u, p, day)
        integ.p = p
        return nothing
    end
    return DiscreteCallback((u, t, integ) -> true, affect!; save_positions=(false, false))
end

"""
    problem_cassava_mealybug(; n=40, interval=60, days=360, seed=0, pmap=(;), kwargs...)
        -> CAProblem

The setup of [Promrak2015](@cite): an `n × n` field of susceptible cassava
plants 1 m apart, free of mealybugs and lacewings, for 360 days, with
lacewing larvae released every `interval` days while surveys find
mealybugs. Nothing lies beyond the field.
"""
function problem_cassava_mealybug(; n::Integer=40, interval::Integer=60, days::Integer=360,
                                  seed::Integer=0, pmap=(;), kwargs...)
    n >= 5 || throw(ArgumentError("the field must be at least 5 × 5, got $n"))
    interval >= 1 || throw(ArgumentError("interval must be at least one day"))
    given = merge((; interval), (; (Symbol(k) => v for (k, v) in pairs(pmap))...))
    return CAProblem(build_cassava_mealybug(), (;), days, given; space=Grid(n, n), seed,
                     callback=cassava_management_callback(), kwargs...)
end

"""
    cassava_yield(u; per_plant=2.25) -> Float64

The estimated yield of the field, in kg: `per_plant` for each plant never
infested for two weeks, a tenth less for those infested two weeks in days 211
to 360, three tenths less for days 121 to 210, nothing for days 1 to 120 or
for removed plants.
"""
function cassava_yield(u; per_plant::Real=2.25)
    share = (1.0, 0.9, 0.7, 0.0)
    return per_plant * sum(u.state[I] == CMB_E ? 0.0 : share[u.loss[I] + 1]
                           for I in eachindex(u.state))
end

_register!(:cassava_mealybug; build=build_cassava_mealybug,
           problem=problem_cassava_mealybug)
