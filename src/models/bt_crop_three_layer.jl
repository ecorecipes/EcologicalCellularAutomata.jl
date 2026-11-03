# Yang et al. (2009): a Bt crop, its target pest, a secondary pest and their natural enemy on three stacked
# lattices. Individuals move and breed in contiguous pairs of cells, enemies eat the pest below them, pests eat
# the plant below them, and individuals age, starve or die of inherited diseases as in the Penna model of He et
# al. (2002). The pairs are a domino tiling: the Margolus partition into 2 x 2 blocks, each block split into
# two horizontal or two vertical pairs at random, so the whole cycle is one block rule.

const _BTC_NONE = (s=Int8(0), age=Int8(0), hunger=Int16(0), genome=UInt32(0))

@inline _btc_pest(c) = (s=c.pest, age=c.pest_age, hunger=c.pest_hunger,
                        genome=c.pest_genome)
@inline function _btc_enemy(c)
    return (s=Int8(c.enemy), age=c.enemy_age, hunger=c.enemy_hunger, genome=c.enemy_genome)
end

# Maturity age and birth interval of a pest (species 1 = B, 2 = C).
@inline _btc_life_pest(p, s) = s == 1 ? (p.R_B, p.r_B) : (p.R_C, p.r_C)
@inline _btc_breeds(x, R, r) = x.age >= R && iszero(mod(x.age - R, r))

# A child: the parent's genome with one more deleterious mutation at a random age (bit) drawn from u.
@inline function _btc_child(x, u)
    return (s=x.s, age=Int8(0), hunger=Int16(0),
            genome=x.genome | (UInt32(1) << min(floor(Int, 32u), 31)))
end

# One pair of cells on one layer: an individual next to an empty cell breeds into it if this is a birth
# cycle, and otherwise moves into it. Pairs of two individuals or of two empty cells do not change.
@inline function _btc_pair(a, b, Rr_a, Rr_b, ua, ub)
    if a.s != 0 && b.s == 0
        return _btc_breeds(a, Rr_a...) ? (a, _btc_child(a, ub)) : (_BTC_NONE, a)
    elseif b.s != 0 && a.s == 0
        return _btc_breeds(b, Rr_b...) ? (_btc_child(b, ua), b) : (b, _BTC_NONE)
    end
    return (a, b)
end

# A layer of a block (top left, bottom left, top right, bottom right) in two horizontal or two vertical pairs.
@inline function _btc_pairs(x, Rr, horizontal, us)
    if horizontal
        a, c = _btc_pair(x[1], x[3], Rr[1], Rr[3], us[1], us[3])
        b, d = _btc_pair(x[2], x[4], Rr[2], Rr[4], us[2], us[4])
    else
        a, b = _btc_pair(x[1], x[2], Rr[1], Rr[2], us[1], us[2])
        c, d = _btc_pair(x[3], x[4], Rr[3], Rr[4], us[3], us[4])
    end
    return (a, b, c, d)
end

@inline _btc_fed(x, h) = (; x..., hunger=Int16(max(x.hunger - h, 0)))
@inline _btc_hungry(x, U) = (; x..., hunger=Int16(x.hunger + U))

# Ageing and the survival check: death at the maximum age, when the inherited diseases expressed up to the
# individual's age reach the threshold T, or when its hunger exceeds hmax.
@inline function _btc_survive(x, T, hmax, max_age)
    x.s == 0 && return x
    age = x.age + 1
    expressed = age >= 32 ? x.genome : x.genome & ((UInt32(1) << age) - UInt32(1))
    (age >= max_age || count_ones(expressed) >= T || x.hunger > hmax) && return _BTC_NONE
    return (; x..., age=Int8(age))
end

# Predation, feeding, survival and plant regrowth in one cell after the moves.
@inline function _btc_cell(plant, pest, enemy, p, u_prey, u_feed)
    if enemy.s != 0
        if pest.s != 0
            P, h = pest.s == 1 ? (p.P_DB, p.h_DB) : (p.P_DC, p.h_DC)
            if u_prey < P
                pest = _BTC_NONE
                enemy = _btc_fed(enemy, h)
            else
                enemy = _btc_hungry(enemy, p.U_D)
            end
        else
            enemy = _btc_hungry(enemy, p.U_D)
        end
    end
    if pest.s != 0
        P, d, h, U = pest.s == 1 ? (p.P_BA, p.d_B, p.h_BA, p.U_B) :
                     (p.P_CA, p.d_C, p.h_CA, p.U_C)
        if plant >= d && u_feed < P
            plant -= d
            pest = _btc_fed(pest, h)
        else
            pest = _btc_hungry(pest, U)
        end
    end
    T, hmax = pest.s == 1 ? (p.T_B, p.hmax_B) : (p.T_C, p.hmax_C)
    pest = _btc_survive(pest, T, hmax, p.max_age)
    enemy = _btc_survive(enemy, p.T_D, p.hmax_D, p.max_age)
    return min(plant + p.regrowth, 1.0), pest, enemy
end

@inline function _btc_update(c, pest, enemy, p, u_prey, u_feed)
    plant, pest, enemy = _btc_cell(c.plant, pest, enemy, p, u_prey, u_feed)
    return _btc_write(c, plant, pest, enemy)
end

@inline function _btc_write(c, plant, pest, enemy)
    return (; c..., plant, pest=pest.s, pest_age=pest.age, pest_hunger=pest.hunger,
            pest_genome=pest.genome, enemy=enemy.s != 0, enemy_age=enemy.age,
            enemy_hunger=enemy.hunger, enemy_genome=enemy.genome)
end

"""
    build_bt_crop_three_layer() -> CAModel

The three-layer model of [Yang2009](@cite), with the individuals of
[He2002](@cite), one cycle a step. Layer M3 holds the plant A (`plant`, in
units of its maximum), layer M2 the target pest B or the secondary pest C
(`pest`: 0 empty, 1 B, 2 C) and layer M1 their natural enemy D (`enemy`).
Every individual has an age, a hunger and a 32-bit Penna genome
(`pest_age`, `pest_hunger`, `pest_genome` and the same for `enemy`). In each
cycle:

1. on M1 and M2 separately, the lattice is cut into pairs of contiguous cells
   (a random domino tiling: the `Margolus` blocks, each split into two
   horizontal or two vertical pairs); an individual next to an empty cell
   moves into it, or, in a birth cycle (age at least its maturity `R`, then
   every `r` cycles), leaves an offspring there with its genome and one more
   deleterious mutation;
2. in every cell, an enemy above a pest eats it with probability `P_DB` (B) or
   `P_DC` (C), and its hunger falls by `h_DB` or `h_DC` (not below 0); a pest
   above at least `d_B` (or `d_C`) of plant eats that much with probability
   `P_BA` (or `P_CA`), and its hunger falls by `h_BA` (or `h_CA`); an
   individual that does not eat grows hungrier by `U_B`, `U_C` or `U_D`;
3. individuals age by one; one dies when its hunger exceeds `hmax_B`,
   `hmax_C` or `hmax_D`, when the diseases its genome expresses up to its age
   reach `T_B`, `T_C` or `T_D`, or at `max_age`; the plant regrows by
   `regrowth`, up to its maximum.

Solve it with the `Margolus()` scheme on a periodic lattice with even sides.
"""
function build_bt_crop_three_layer()
    rule = BlockRule() do blk, ctx
        p = ctx.p
        before = map(_btc_pest, blk)
        pests = _btc_pairs(before, map(x -> _btc_life_pest(p, x.s), before),
                           uniform(ctx, 1) < 0.5, ntuple(i -> uniform(ctx, 10i + 1), 4))
        Rr_D = (p.R_D, p.r_D)
        enemies = _btc_pairs(map(_btc_enemy, blk), (Rr_D, Rr_D, Rr_D, Rr_D),
                             uniform(ctx, 2) < 0.5, ntuple(i -> uniform(ctx, 10i + 2), 4))
        return ntuple(i -> _btc_update(blk[i], pests[i], enemies[i], p,
                                       uniform(ctx, 10i + 3), uniform(ctx, 10i + 4)), 4)
    end
    states = [StateVariable(:plant, 1.0;
                            description="quantity of the plant A, in units of its maximum"),
              StateVariable(:pest, Int8(0);
                            description="0 empty, 1 target pest B, 2 secondary pest C"),
              StateVariable(:pest_age, Int8(0); description="age of the pest"),
              StateVariable(:pest_hunger, Int16(0); description="hunger of the pest"),
              StateVariable(:pest_genome, UInt32(0);
                            description="Penna genome of the pest"),
              StateVariable(:enemy, false; description="a natural enemy D"),
              StateVariable(:enemy_age, Int8(0); description="age of the enemy"),
              StateVariable(:enemy_hunger, Int16(0); description="hunger of the enemy"),
              StateVariable(:enemy_genome, UInt32(0);
                            description="Penna genome of the enemy")]
    parameters = [Parameter(:P_BA, 0.6;
                            description="probability that B feeds on the plant"),
                  Parameter(:P_CA, 0.5;
                            description="probability that C feeds on the plant"),
                  Parameter(:P_DB, 0.6; description="probability that D eats B"),
                  Parameter(:P_DC, 0.6; description="probability that D eats C"),
                  Parameter(:d_B, 0.2; description="plant eaten by B in one meal"),
                  Parameter(:d_C, 0.05; description="plant eaten by C in one meal"),
                  Parameter(:regrowth, 0.1; description="plant regrowth in a cycle"),
                  Parameter(:h_BA, 3; description="fall in B's hunger after a meal"),
                  Parameter(:h_CA, 3; description="fall in C's hunger after a meal"),
                  Parameter(:h_DB, 3; description="fall in D's hunger after eating B"),
                  Parameter(:h_DC, 3; description="fall in D's hunger after eating C"),
                  Parameter(:U_B, 1;
                            description="rise in B's hunger in a cycle without food"),
                  Parameter(:U_C, 1;
                            description="rise in C's hunger in a cycle without food"),
                  Parameter(:U_D, 1;
                            description="rise in D's hunger in a cycle without food"),
                  Parameter(:hmax_B, 2; description="hunger B endures"),
                  Parameter(:hmax_C, 2; description="hunger C endures"),
                  Parameter(:hmax_D, 3; description="hunger D endures"),
                  Parameter(:R_B, 3; description="maturity age of B"),
                  Parameter(:R_C, 3; description="maturity age of C"),
                  Parameter(:R_D, 4; description="maturity age of D"),
                  Parameter(:r_B, 1; description="cycles between births of B"),
                  Parameter(:r_C, 1; description="cycles between births of C"),
                  Parameter(:r_D, 4; description="cycles between births of D"),
                  Parameter(:T_B, 3;
                            description="inherited diseases B survives (Penna threshold)"),
                  Parameter(:T_C, 3;
                            description="inherited diseases C survives (Penna threshold)"),
                  Parameter(:T_D, 7;
                            description="inherited diseases D survives (Penna threshold)"),
                  Parameter(:max_age, 32;
                            description="maximum age (at most 32, the genome's length)")]
    return CAModel(rule; name=:bt_crop_three_layer, states, parameters,
                   description="Yang et al. (2009) Bt crop, target and secondary pests and their enemy on three layers")
end

"""
    problem_bt_crop_three_layer(; n=100, occupancy=0.5, cycles=1000, seed=0, pmap=(;),
                                kwargs...) -> CAProblem

The setup of [Yang2009](@cite) on an `n × n` torus (`n` even; the source uses
500 and 2000 cycles): the plant at its maximum everywhere, an enemy in each
cell of M1 with probability `occupancy`, and a pest in each cell of M2 with
probability `occupancy`, B or C with equal probability. Individuals start
with clean genomes, no hunger and ages drawn uniformly below their maturity
age (and the maximum age). Drawn from `seed`; parameters through `pmap`. Solve with
`solve(prob, Margolus())`.
"""
function problem_bt_crop_three_layer(; n::Integer=100, occupancy::Real=0.5,
                                     cycles::Integer=1000, seed::Integer=0, pmap=(;),
                                     kwargs...)
    (n >= 2 && iseven(n)) || throw(ArgumentError("n must be even and at least 2, got $n"))
    0 <= occupancy <= 1 ||
        throw(ArgumentError("occupancy must lie in [0, 1], got $occupancy"))
    model = build_bt_crop_three_layer()
    given = (; (Symbol(a) => v for (a, v) in pairs(pmap))...)
    p = merge(defaults(model), given)
    p.max_age <= 32 || throw(ArgumentError("max_age must be at most 32, got $(p.max_age)"))
    rng = Random.Xoshiro(hash((seed, 0, :bt_crop_three_layer)))
    enemy = [rand(rng) < occupancy for _ in 1:n, _ in 1:n]
    pest = [rand(rng) < occupancy ? Int8(rand(rng) < 0.5 ? 1 : 2) : Int8(0)
            for _ in 1:n,
                _ in 1:n]
    young(R) = Int8(rand(rng, 0:(min(R, p.max_age) - 1)))
    pest_age = [s == 0 ? Int8(0) : young(s == 1 ? p.R_B : p.R_C) for s in pest]
    enemy_age = [e ? young(p.R_D) : Int8(0) for e in enemy]
    u0 = (; plant=ones(n, n), pest, pest_age, pest_hunger=zeros(Int16, n, n),
          pest_genome=zeros(UInt32, n, n), enemy, enemy_age,
          enemy_hunger=zeros(Int16, n, n),
          enemy_genome=zeros(UInt32, n, n))
    return CAProblem(model, u0, cycles, given; space=Torus(n, n), seed, kwargs...)
end

"""
    bt_crop_summary(u) -> NamedTuple

The source's outputs for the state `u`: the mean quantity of the plant
(`plant`, as a share of its maximum) and the numbers of target pests (`B`),
secondary pests (`C`) and enemies (`D`).
"""
function bt_crop_summary(u)
    return (; plant=sum(u.plant) / length(u.plant), B=count(==(1), u.pest),
            C=count(==(2), u.pest), D=count(u.enemy))
end

_register!(:bt_crop_three_layer; build=build_bt_crop_three_layer,
           problem=problem_bt_crop_three_layer)
