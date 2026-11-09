# The models, each against its source's equations: one-step transition frequencies over many cells or
# seeds against the probabilities written out from the source (never from the package's code paths),
# plus structural checks of states, layouts and helpers. Seeds are fixed.

# |frequency - p| within four binomial standard errors, plus three counts for rare outcomes, where
# the normal approximation fails (one hit where 0.01 are expected happens one time in a hundred)
close_to(freq, p, n) = abs(freq - p) <= 4 * sqrt(p * (1 - p) / n) + 3 / n

@testset "fungus_dispersal" begin
    w, total = fungus_kernel(3.0; radius=4)
    @test size(w) == (9, 9) && w == reverse(w) && w == permutedims(w)
    f(d) = exp(-d / 3) / (2π * 9)
    @test w[5, 5] ≈ -log(1 - f(0)) && w[5, 7] ≈ -log(1 - f(2)) && total ≈ sum(w)
    @test_throws ArgumentError fungus_kernel(-1.0)

    # one infectious host in the middle of a 5 x 5 field: after one step each healthy host at
    # distance d is latent with probability F (1 - (1 - f(d))^c)
    model = build_fungus_dispersal(; radius=4)
    host = zeros(Int8, 5, 5)
    host[3, 3] = 2
    n = 3000
    hits = zeros(Int, 5, 5)
    for s in 1:n
        prob = CAProblem(model, (; host), 1, (F_permissive=0.5, load=40.0);
                         space=Grid(5, 5),
                         seed=s)
        hits .+= solve(prob).u[end].host .== 1
    end
    for I in CartesianIndices(hits)
        I == CartesianIndex(3, 3) && continue
        d = hypot(Tuple(I - CartesianIndex(3, 3))...)
        @test close_to(hits[I] / n, 0.5 * (1 - (1 - f(d))^40), n)
    end

    # a host infected in the step from day t is infectious from day t + latency
    prob = CAProblem(model, (; host), 6, (F_permissive=1.0, load=1e6, latency=4);
                     space=Grid(5, 5), seed=1)
    sol = solve(prob)
    @test sol.u[2].host[3, 4] == 1 && sol.u[2].clock[3, 4] == 1
    @test sol.u[4].host[3, 4] == 1 && sol.u[5].host[3, 4] == 2
    @test all(sol.u[end].host .>= 1)

    # the temperate season has no load far from its peak; the tropical one has
    temperate = solve(problem(:fungus_dispersal; seasonality=:temperate, tspan=10, seed=1))
    @test fungus_infected(temperate.u[end]) == 9
    tropical = solve(problem(:fungus_dispersal; seasonality=:tropical, tspan=30, seed=1))
    @test fungus_infected(tropical.u[end]) > 9

    @test count(fungus_hosts(60; layout=:bands)) == 1800
    @test fungus_hosts(60; layout=:bands)[1, 16] && !fungus_hosts(60; layout=:bands)[1, 15]
    @test !any(fungus_hosts(10)) && all(fungus_hosts(10; layout=:resistant))
    @test 0.4 < count(fungus_hosts(60; layout=:random)) / 3600 < 0.6
    @test_throws ArgumentError fungus_hosts(10; layout=:checkerboard)
    @test_throws ArgumentError build_fungus_dispersal(; seasonality=:monsoon)
    p0 = problem(:fungus_dispersal; layout=:resistant)
    @test p0.u0.host[30:32, 30:32] == fill(2, 3, 3) && count(p0.u0.host .== 2) == 9 &&
          all(p0.u0.resistant)
end

@testset "rhizoctonia_contact" begin
    z = rhizoctonia_pathozone([(10.0, 10.0), (10.0, 13.0), (55.0, 190.0)])
    @test z[1, 1] ≈ 1 + (1 + 0.49 * 3) * exp(-0.277 * 3)
    @test z[3, 10] ≈ (1 + 0.49 * 5) * exp(-0.277 * 5) && count(!iszero, z) == 2
    @test_throws ArgumentError rhizoctonia_pathozone([(120.0, 10.0)])
    @test length(rhizoctonia_inoculum(Random.Xoshiro(1); units=7)) == 7

    # secondary infection only (no inoculum), one infectious neighbour, at day τ_p:
    # P(S -> L) = 1 - exp(-C_s exp(-σ_s x_d²) exp(-d_s τ_p) dt)
    p = defaults(build_rhizoctonia_contact())
    Cs = p.a_s * p.d_s * exp(p.d_s * p.τ_s)
    rs = Cs * exp(-p.σ_s * p.x_d^2) * exp(-p.d_s * p.τ_p)
    expected1 = 1 - exp(-rs * p.dt)
    expected4 = 1 - exp(-4rs * p.dt)
    n = 4000
    u = fill(Int8(2), 3, 3)
    u[2, 2] = 0
    line = Int8[0 2 0]
    one_hits = any_hits = four_hits = 0
    for s in 1:n
        cov = (; pathozone=zeros(1, 3))
        sol = solve(CAProblem(build_rhizoctonia_contact(), (; plant=line), 1;
                              space=Grid(1, 3),
                              seed=s, covariates=cov))
        one_hits += sol.u[end][1, 1] == 1
        cov3 = (; pathozone=zeros(3, 3))
        four_hits += solve(CAProblem(build_rhizoctonia_contact(), (; plant=u), 1;
                                     space=Grid(3, 3), seed=s, covariates=cov3)).u[end][2,
                                                                                        2] ==
                     1
        any_hits += solve(CAProblem(build_rhizoctonia_contact(; secondary=:any),
                                    (; plant=u),
                                    1; space=Grid(3, 3), seed=s,
                                    covariates=cov3)).u[end][2, 2] == 1
    end
    @test close_to(one_hits / n, expected1, n)
    @test close_to(four_hits / n, expected4, n)
    @test close_to(any_hits / n, expected1, n)

    # primary infection only: a unit at the plant itself, at day τ_p
    Cp = p.a_p * p.d_p * exp(p.d_p * p.τ_p)
    expected_p = 1 - exp(-Cp * exp(-p.d_p * p.τ_p) * p.dt)
    hits = count(1:n) do s
        sol = solve(CAProblem(build_rhizoctonia_contact(), (;), 1; space=Grid(1, 1), seed=s,
                              covariates=(; pathozone=ones(1, 1))))
        return sol.u[end][1, 1] == 1
    end
    @test close_to(hits / n, expected_p, n)

    sol = solve(problem(:rhizoctonia_contact; seed=3))
    @test rhizoctonia_days(sol)[1] ≈ 4.898 && rhizoctonia_days(sol)[end] ≈ 4.898 + 15.2
    y = rhizoctonia_infected.(sol.u)
    @test issorted(y) && y[1] == 0 && 0 < y[end] <= 50
    @test_throws ArgumentError build_rhizoctonia_contact(; secondary=:all)
end

@testset "domany_kinzel_viral" begin
    R, π, p, q = 6, 0.3, 0.2, 0.05
    for (r1, r2) in ((0, 0), (1, 0), (2, 3), (6, 6), (6, 1))
        probs = [viral_transition(r, r1, r2; π, p, q, R) for r in 0:R]
        @test sum(probs) ≈ 1
        # every cell of a ring alternating r1, r2 has the parents (r1, r2) or (r2, r1)
        L = 20_000
        u0 = Int8[isodd(i) ? r1 : r2 for i in 1:L]
        prob = CAProblem(build_domany_kinzel_viral(; R), (; r=u0), 1, (; π, p, q);
                         space=Torus(L), seed=5)
        u = solve(prob).u[end]
        for r in 0:R
            @test close_to(count(==(r), u) / L, probs[r + 1], L)
        end
    end
    # the directed bond percolation limit: p = q = 0 and only class R
    @test viral_transition(R, R, 0; π, p=0, q=0, R) ≈ 1 - π^R
    @test viral_transition(R, R, R; π, p=0, q=0, R) ≈ 1 - π^(2R)

    u = Int8[0, 2, 4, 0]
    @test viral_density(u) == 0.5 && viral_mean_ability(u) == 3 && viral_diversity(u) == 1
    @test isnan(viral_mean_ability(zeros(Int8, 3))) &&
          isnan(viral_diversity(zeros(Int8, 3)))
    @test problem(:domany_kinzel_viral; L=10, saveat=5).kwargs == (saveat=5,)
    @test_throws ArgumentError build_domany_kinzel_viral(; R=0)
    @test_throws ArgumentError problem(:domany_kinzel_viral; L=1)
end

@testset "bt_cotton_refuge" begin
    model = build_bt_cotton_refuge()
    p = defaults(model)
    n = 4000
    # on a 3 x 3 torus the centre's Moore neighbours are the other eight cells
    function centre_outcomes(u, pmap)
        out = zeros(Int, 3)
        for s in 1:n
            v = solve(CAProblem(model, (; plant=u), 1, pmap; space=Torus(3, 3), seed=s)).u[end]
            out[v[2, 2] + 1] += 1
        end
        return out ./ n
    end
    # an empty cell with two S. frugiperda and one H. armigera neighbour
    u = Int8[1 1 2; 0 0 0; 0 0 0]
    w1, w2 = 2 * p.fecundity_Sf, 1 * p.fecundity_Ha
    f1, f2 = p.emergence * w1 / (w1 + w2), p.emergence * w2 / (w1 + w2)
    freq = centre_outcomes(u, (;))
    @test close_to(freq[2], f1, n) && close_to(freq[3], (1 - f1) * f2, n)
    # an S. frugiperda larva at low density with one H. armigera neighbour, unsprayed and sprayed
    u = Int8[2 0 0; 0 1 0; 0 0 0]
    meet = 1 / 8
    survive = (1 - p.mortality) * (meet * p.sL_Sf_Ha + (1 - meet) * p.sL_Sf)
    @test close_to(centre_outcomes(u, (;))[2], survive, n)
    @test close_to(centre_outcomes(u, (; spray=true))[2], survive * (1 - p.kill), n)
    @test close_to(centre_outcomes(u, (; spray=true, resistant_Sf=true))[2], survive, n)
    @test close_to(centre_outcomes(u, (; spray=true, resistant_Sf=true, escape=0.5))[2],
                   survive * (1 - p.kill * 0.5), n)
    # an H. armigera larva at high density (four occupied neighbours, three of them S. frugiperda)
    u = Int8[1 1 1; 2 2 0; 0 0 0]
    meet = 3 / 8
    survive = (1 - p.mortality) * (meet * p.sH_Ha_Sf + (1 - meet) * p.sH_Ha)
    @test close_to(centre_outcomes(u, (;))[3], survive, n)

    # the spray callback reads the state after each generation
    prob = problem(:bt_cotton_refuge; n=20, initial=100, scenario=:insecticide, tspan=1,
                   control_level=0.0)
    integ = init(prob)
    step!(integ)
    @test integ.p.spray
    integ = init(problem(:bt_cotton_refuge; n=20, initial=100, tspan=1, control_level=0.0))
    step!(integ)
    @test !integ.p.spray
    @test problem(:bt_cotton_refuge; n=20, initial=10, scenario=:Ha_resistant).p.resistant_Ha
    @test !problem(:bt_cotton_refuge; n=20, initial=10, scenario=:Ha_resistant).p.resistant_Sf
    u0 = problem(:bt_cotton_refuge; n=50, initial=200, seed=2).u0
    @test count(==(1), u0) == 200 && count(==(2), u0) == 200
    @test bt_cotton_dominance(Int8[1, 1, 2, 0]) == (Sf=2 / 3, Ha=1 / 3)
    @test all(isnan, bt_cotton_dominance(zeros(Int8, 4)))
    @test bt_cotton_survival(:table1).sL_Ha_Sf == 0.90
    @test_throws ArgumentError bt_cotton_survival(:table2)
    @test_throws ArgumentError problem(:bt_cotton_refuge; scenario=:organic)
    @test_throws ArgumentError problem(:bt_cotton_refuge; n=10, initial=60)
end

@testset "coffee_rust" begin
    n = 4000
    # a susceptible farm with N infected Moore neighbours: 0.1 (N + 1) / 10 (archived), N / 8 (paper)
    for (N, rule, expected) in ((3, :archived, 0.1 * 4 / 10), (8, :archived, 0.09),
                                (3, :paper, 3 / 8))
        coffee = fill(Int8(1), 3, 3)
        coffee[[1, 2, 3, 4, 6, 7, 8, 9][1:N]] .= 2
        model = build_coffee_rust(; ca_rule=rule)
        hits = count(1:n) do s
            u = solve(CAProblem(model, (; coffee), 1; space=Grid(3, 3), seed=s)).u[end]
            return u.coffee[2, 2] == 2
        end
        @test close_to(hits / n, expected, n)
    end
    # a spore that cannot move infects next to where it was released
    coffee = Int8[2 0 1]
    @test coffee_rust_spores!(coffee, [0.0 0.5 0.0], Random.Xoshiro(1);
                              spore_infection=1.0) ==
          1
    @test coffee == Int8[2 0 2]
    coffee = Int8[2 0 1]
    @test coffee_rust_spores!(coffee, [0.0 0.5 0.0], Random.Xoshiro(1);
                              spore_infection=0.0) ==
          0
    # a spore walks 1000 steps along a corridor of five matrix cells from its second cell; by
    # parity it ends on cell 2, 4 or 6 with probabilities 1:2:1 (the degrees), so it lands next to
    # the far farm with probability 1/4
    coffee = Int8[2 0 0 0 0 0 1]
    canopy = [0.0 1.0 1.0 1.0 1.0 1.0 0.0]
    far = count(1:2000) do s
        c = copy(coffee)
        coffee_rust_spores!(c, canopy, Random.Xoshiro(s); credit=1000.0,
                            spore_infection=1.0)
        return c[1, 7] == 2
    end
    @test close_to(far / 2000, 1 / 4, 2000)
    @test coffee_rust_spores!(copy(coffee), canopy, Random.Xoshiro(1); credit=0.5,
                              spore_infection=1.0) == 0

    a = rand(Random.Xoshiro(2), 40, 40)
    @test count(==(0), classify_landscape(a, [0.25, 0.75])) == 400
    rng = Random.Xoshiro(3)
    lc = random_cluster_landscape(rng, 30, 0.2)
    le = random_element_landscape(rng, 30, 50)
    @test extrema(lc) == (0.0, 1.0) && extrema(le) == (0.0, 1.0)
    @test length(unique(le)) == 50
    # nearest seeds, ties included, as SciPy's distance_transform_edt gives them (SciPy 1.17)
    mask = falses(7, 9)
    foreach(S -> mask[S...] = true, [(1, 2), (4, 3), (4, 6), (5, 3), (6, 3), (6, 5)])
    R, C = EcologicalCellularAutomata._nearest_seed(mask)
    @test R == [1 1 1 1 1 4 4 4 4; 1 1 1 1 4 4 4 4 4; 1 4 4 4 4 4 4 4 4; 4 4 4 4 4 4 4 4 4;
                5 5 5 5 6 4 4 4 4; 6 6 6 6 6 6 6 4 4; 6 6 6 6 6 6 6 6 6]
    @test C == [2 2 2 2 2 6 6 6 6; 2 2 2 2 6 6 6 6 6; 2 3 3 3 6 6 6 6 6; 3 3 3 3 6 6 6 6 6;
                3 3 3 3 5 6 6 6 6; 3 3 3 3 5 5 5 6 6; 3 3 3 3 5 5 5 5 5]
    # against a search of every seed: each cell gets one of its nearest seeds
    m = rand(Random.Xoshiro(6), 40, 40) .< 0.05
    R, C = EcologicalCellularAutomata._nearest_seed(m)
    seeds = findall(m)
    @test all(CartesianIndices(m)) do I
        m[R[I], C[I]] &&
            (R[I] - I[1])^2 + (C[I] - I[2])^2 ==
            minimum(sum(abs2, Tuple(I - S)) for S in seeds)
    end
    coffee, canopy = coffee_rust_landscape(Random.Xoshiro(4); n=60)
    @test 0.2 < count(==(1), coffee) / 3600 < 0.3
    @test all(canopy[coffee .== 1] .== 0) && all(0 .<= canopy .<= 1)
    @test_throws ArgumentError random_cluster_landscape(rng, 10, 1.5)
    @test_throws ArgumentError build_coffee_rust(; ca_rule=:beta)

    prob = problem(:coffee_rust; n=40, buffer=5, tspan=20, seed=1)
    @test haskey(prob.kwargs, :callback)
    @test !haskey(problem(:coffee_rust; n=40, buffer=5, spores=false).kwargs, :callback)
    sol = solve(prob)
    @test coffee_rust_prevalence(sol.u[end]; buffer=5) >=
          coffee_rust_prevalence(sol.u[1]; buffer=5) > 0
    @test coffee_rust_rate(sol; buffer=5) >= 0
    # the spores draw from the problem's seed: the same seed gives the same run
    @test solve(prob).u[end].coffee == sol.u[end].coffee

    # archived spores (LandscapeSim.py), on a 16 x 16 landscape, 0-based (r, c): a spore is kept
    # only if r < (16 & c) < 16; 16 & c is 0 on the landscape, so only spores off the top edge stay
    coffee = fill(Int8(1), 16, 16)
    canopy = zeros(16, 16)
    coffee[4, 6] = Int8(2)                # (3, 5): its only free neighbour (4, 5) is dropped
    coffee[5, 6] = 0
    canopy[5, 6] = 0.5
    @test coffee_rust_spores_archived!(copy(coffee), canopy, Random.Xoshiro(1);
                                       spore_infection=1.0) == 0
    # (0, 5) releases off the edge into row -1, where no move is open: it steps by (1, 1) onto
    # coffee at (0, 6), stops, and infects one of that cell's susceptible neighbours with
    # probability 1/2 (the target is relative to the last spore, here the same one)
    coffee = fill(Int8(1), 16, 16)
    coffee[1, 6] = Int8(2)
    hits = [coffee_rust_spores_archived!(copy(coffee), canopy, Random.Xoshiro(k))
            for k in 1:4000]
    @test all(in((0, 1)), hits) && close_to(count(==(1), hits) / 4000, 0.5, 4000)
    # with archived spores the first farm may lie in the buffer, as in that code
    @test problem(:coffee_rust; n=40, buffer=5, spores=:archived).kwargs[:callback] isa
          DiscreteCallback
    @test_throws ArgumentError problem(:coffee_rust; n=40, buffer=5, spores=:walk)
    @test_throws ArgumentError problem(:coffee_rust; n=40, buffer=5, start=:edge)
end

@testset "huffaker_mites" begin
    p = defaults(build_huffaker_mites())
    state(mite; energy=nothing, post=fill(false, size(mite))) = (; mite=Int8.(mite),
                                                                 fertility=zeros(Int16,
                                                                                 size(mite)),
                                                                 energy=energy === nothing ?
                                                                        Int16.(ifelse.(mite .==
                                                                                       1,
                                                                                       10,
                                                                                       mite)) :
                                                                        Int16.(energy),
                                                                 post)
    # a predator next to a prey eats it, moving into its cell and gaining energy
    u = state([2 1]; energy=[1 10])
    huffaker_step!(u, p, Random.Xoshiro(1))
    @test u.mite == Int8[0 2] && u.energy[1, 2] == 2 && u.fertility[1, 2] == 1
    # a predator with one life point that moves without eating dies where it stood
    u = state([2 0]; energy=[1 0])
    huffaker_step!(u, p, Random.Xoshiro(1))
    @test u.mite == Int8[0 0]
    # a blocked predator loses nothing
    u = state([2 2]; energy=[1 1])
    huffaker_step!(u, p, Random.Xoshiro(1))
    @test u.mite == Int8[2 2] && u.energy == Int16[1 1]
    # a prey that moves at its fertility threshold leaves an offspring behind
    u = state([1 0])
    huffaker_step!(u, p, Random.Xoshiro(1))
    @test u.mite == Int8[1 1] && u.fertility == Int16[0 0] && u.energy == Int16[10 10]
    # a predator below its threshold moves without breeding
    u = state([2 0 1]; energy=[3 0 10])
    huffaker_step!(u, merge(p, (; F_predator=5)), Random.Xoshiro(1))
    @test count(==(2), u.mite) == 1
    # posts: with the code's rule any empty post is a possible move
    to_post(posts) = count(1:2000) do s
        u = state([1 0 0 0 0]; post=Bool[0 0 0 0 1])
        huffaker_step!(u, merge(p, (; F_prey=100)), Random.Xoshiro(s); posts)
        return u.mite[1, 5] == 1
    end
    @test to_post(:none) == 0
    @test close_to(to_post(:code) / 2000, 0.5, 2000)
    @test to_post(:text) == 0          # the prey is not on a post
    from_post = count(1:2000) do s
        u = state([0 1 0 0 0]; post=Bool[0 1 0 0 1])
        huffaker_step!(u, merge(p, (; F_prey=100)), Random.Xoshiro(s); posts=:text)
        return u.mite[1, 5] == 1
    end
    @test close_to(from_post / 2000, 1 / 3, 2000)
    @test_throws ArgumentError huffaker_step!(state([1 0]), p, Random.Xoshiro(1);
                                              posts=:all)

    @test count_peaks([sin(2π * t / 50) + 1 for t in 0:499]) == 10
    @test count_peaks(ones(100)) == 0 && count_peaks(collect(1.0:100)) == 0
    @test count_peaks(Float64[]) == 0
    osc = [100 + 50sin(2π * t / 50) for t in 0:499]
    @test huffaker_outcome(osc, fill(5, 500)) == :no_extinction
    @test huffaker_outcome(vcat(osc[1:400], zeros(100)), zeros(500)) == :both_extinct
    @test huffaker_outcome(vcat(collect(1.0:50), zeros(450)), zeros(500)) ==
          :both_extinct_immediately
    @test huffaker_outcome(vcat(osc[1:400], zeros(100)), fill(3, 500)) == :prey_extinct
    @test huffaker_outcome(osc, vcat(fill(3, 400), zeros(100))) == :predator_extinct
    @test huffaker_outcome(fill(10.0, 500), zeros(500)) == :predator_extinct_immediately

    prob = problem(:huffaker_mites; experiment=3, seed=1, tspan=5)
    @test size(prob.u0.mite) == (30, 36) && count(prob.u0.post) == 18
    @test count(==(1), prob.u0.mite) == 120 && count(==(2), prob.u0.mite) == 27
    @test size(problem(:huffaker_mites; experiment=1).u0.mite) == (20, 24)
    @test problem(:huffaker_mites; experiment=4).p.F_predator == 20
    @test problem(:huffaker_mites; experiment=4).u0.energy[findfirst(==(2),
                                                                     problem(:huffaker_mites;
                                                                             experiment=4).u0.mite)] ==
          2
    sol = solve(prob)
    prey, pred = huffaker_counts(sol)
    @test length(prey) == 6 && prey[1] == 120 && pred[1] == 27
    @test solve(prob).u[end].mite == sol.u[end].mite
    @test_throws ArgumentError problem(:huffaker_mites; experiment=5)
end

@testset "cvc_citrus" begin
    K = cvc_flight_kernel(; radius=20, samples=10^5)
    @test size(K) == (41, 41) && 0.5 < sum(K) <= 1 && all(>=(0), K)
    @test_throws ArgumentError cvc_flight_kernel(; stability=1)
    # a hand-made kernel: each flight lands on one of the eight neighbours with probability 0.2...
    k = fill(0.2, 3, 3)
    k[2, 2] = 0.0
    model = build_cvc_citrus(; kernel=k, method=:direct)
    diseased = fill(false, 5, 5)
    diseased[3, 3] = true
    n = 4000
    pmin = 0.001
    # ... and at t = 0 there is n_v = 5 + trunc(5 sin(-π/3)) = 1 flight per diseased tree
    pnst0 = pmin + (1 - 2pmin) * (1 + sin(-π / 3)) / 2
    infected = stressed_far = 0
    for s in 1:n
        u = solve(CAProblem(model, (; diseased), 1; space=Torus(5, 5), seed=s)).u[end]
        infected += u.diseased[3, 4]
        stressed_far += u.stressed[1, 1]
    end
    @test close_to(infected / n, 0.2, n)
    @test close_to(stressed_far / n, 1 - pnst0, n)
    # a stressed tree is never infected
    stressed = fill(true, 5, 5)
    hits = count(s -> solve(CAProblem(model, (; diseased, stressed), 1; space=Torus(5, 5),
                                      seed=s)).u[end].diseased[3, 4], 1:500)
    @test hits == 0
    @test_throws ArgumentError build_cvc_citrus(; kernel=fill(1.0, 3, 3))

    prob = problem(:cvc_citrus; grove=:SJ71, L=40, kernel=K, tspan=3, seed=2)
    @test size(prob.u0.diseased) == (80, 80) && !any(prob.u0.diseased[41:end, :])
    @test prob.p.p_min == 0.35 && count(prob.mask) == 1600
    sol = solve(prob)
    @test !any(sol.u[end].diseased[41:end, :]) && !any(sol.u[end].diseased[:, 41:end])
    @test cvc_incidence(sol.u[end]) >= cvc_incidence(sol.u[1]) > 0
    @test cvc_field_curve(:SJ01, 20) ≈ 0.218 atol = 0.001
    @test_throws ArgumentError cvc_field_curve(:SJ99, 1)
    @test_throws ArgumentError problem(:cvc_citrus; grove=:SJ02)
    @test_throws ArgumentError problem(:cvc_citrus; L=10, kernel=K)
end

@testset "coffee_ant_fungus" begin
    n = 4000
    # on a 3 x 3 torus the centre's Moore neighbours are the other eight trees
    function centre(u, model, pmap)
        out = zeros(Int, 3)
        for s in 1:n
            v = solve(CAProblem(model, (; tree=u), 1, pmap; space=Torus(3, 3), seed=s)).u[end]
            out[v[2, 2] + 1] += 1
        end
        return out ./ n
    end
    one = build_coffee_ant_fungus(; phases=:one)
    p = defaults(one)
    u = Int8[1 1 1; 0 0 0; 0 0 0]                  # empty centre, three nests around
    @test close_to(centre(u, one, (;))[2], p.s0 + 3p.s1, n)
    u = Int8[2 2 1; 0 1 0; 0 0 0]                  # a nest with two infected neighbours
    @test close_to(centre(u, one, (;))[3], min(1, p.t0 + 2p.t1), n)
    @test close_to(centre(u, one, (; t1=0.2))[3], 0.4, n)
    u = Int8[1 0 0; 0 1 0; 0 0 0]                  # no infected neighbour: background infection
    @test close_to(centre(u, one, (; a=0.3))[3], 0.3, n)
    u = Int8[0 0 0; 0 2 0; 0 0 0]                  # an infected nest dies with probability v
    @test close_to(centre(u, one, (;))[1], p.v, n)
    phorid = build_coffee_ant_fungus(; mortality=:phorid, phases=:one)
    u = Int8[1 1 0; 0 1 0; 0 0 0]                  # a nest with two neighbours under phorid attack
    @test close_to(centre(u, phorid, (;))[1], p.m0 + 2p.m1, n)

    # two phases: the even step only colonises, the odd step only infects and kills
    two = build_coffee_ant_fungus()
    u = Int8[2 2 0; 0 1 0; 0 0 2]
    sure = (s0=1.0, v=1.0, t1=1.0)
    sol = solve(CAProblem(two, (; tree=u), 2, sure; space=Torus(3, 3), seed=1))
    @test all(!=(0), sol.u[2]) && count(==(2), sol.u[2]) == 3   # step 0: everything colonised
    @test all(sol.u[3][u .== 2] .== 0)                          # step 1: the infected nests died
    @test sol.u[3][2, 2] == 2                                   # and the nest was infected
    @test ant_iterations(10) == 5 && ant_iterations(10; phases=:one) == 10

    @test nest_density(Int8[0 1; 2 0]) == 0.5
    clumped = zeros(Int8, 6, 6)
    clumped[1:3, 1:3] .= 1
    @test mean_variance_ratio(clumped) < 1
    @test isnan(mean_variance_ratio(ones(Int8, 6, 6)))
    prob = problem(:coffee_ant_fungus; iterations=10, seed=2)
    @test prob.tspan == (0, 20) && size(prob.u0) == (90, 120) && !any(==(2), prob.u0)
    @test problem(:coffee_ant_fungus; iterations=10, phases=:one).tspan == (0, 10)
    @test_throws ArgumentError build_coffee_ant_fungus(; mortality=:birds)
    @test_throws ArgumentError build_coffee_ant_fungus(; phases=:three)
end

@testset "pine_wilt_casim" begin
    n = 4000
    @test pine_wilt_weights(5) == [1, 0.8, 0.5, 0.3, 0.1]
    @test_throws ArgumentError pine_wilt_weights(6)
    # on a 5 x 5 torus the centre's rings hold exactly 8 and 16 cells:
    # two infected in ring 1 and four in ring 2 give p = 2/8 + 0.8 * 4/16 = 0.45 with Moore 2
    cover = fill(Int8(1), 5, 5)
    cover[3, 3] = 2
    cover[2, 3] = 4
    cover[3, 4] = 3
    cover[1, 1] = cover[1, 5] = cover[5, 1] = cover[5, 5] = 4
    model = build_pine_wilt_casim(; weights=pine_wilt_weights(2))
    hits = count(1:n) do s
        u = solve(CAProblem(model, (; cover), 1, (; α=0.0); space=Torus(5, 5), seed=s)).u[end]
        return u.cover[3, 3] == 4
    end
    @test close_to(hits / n, 0.45, n)
    # with incubation a new infection is latent with probability p_latent
    model1 = build_pine_wilt_casim(; incubation=true)
    c1 = fill(Int8(4), 3, 3)
    c1[2, 2] = 2
    latent = count(s -> solve(CAProblem(model1, (; cover=c1), 1, (; α=0.0);
                                        space=Torus(3, 3),
                                        seed=s)).u[end].cover[2, 2] == 3, 1:n)
    @test close_to(latent / n, 0.8, n)
    # a latent pine stays latent with probability p_stay, for at most max_latent years
    c2 = fill(Int8(1), 3, 3)
    c2[2, 2] = 3
    stays = count(s -> solve(CAProblem(model1, (; cover=c2), 1; space=Torus(3, 3), seed=s)).u[end].cover[2,
                                                                                                         2] ==
                       3, 1:n)
    @test close_to(stays / n, 0.1, n)
    forced = solve(CAProblem(model1, (; cover=c2, years=fill(Int8(2), 3, 3)), 1,
                             (; p_stay=1.0);
                             space=Torus(3, 3), seed=1)).u[end]
    @test forced.cover[2, 2] == 4
    # a symptomatic pine is felled with probability α, and otherwise dies after four years
    c3 = fill(Int8(1), 3, 3)
    c3[2, 2] = 4
    felled = count(s -> solve(CAProblem(model1, (; cover=c3), 1, (; α=0.3);
                                        space=Torus(3, 3),
                                        seed=s)).u[end].cover[2, 2] == 1, 1:n)
    @test close_to(felled / n, 0.3, n)
    sol = solve(CAProblem(model1, (; cover=c3), 4, (; α=0.0); space=Torus(3, 3), seed=1))
    @test [u.cover[2, 2] for u in sol.u] == [4, 4, 4, 4, 5]

    land = pine_landscape(Random.Xoshiro(1); n=50)
    @test 0.55 < count(==(2), land) / 2500 < 0.65 && all(in((0, 1, 2)), land)
    prob = problem(:pine_wilt_casim; n=60, symptomatic=100, latent=20, seed=1)
    @test count(==(4), prob.u0.cover) == 100 && count(==(3), prob.u0.cover) == 20
    @test pine_wilt_infected(prob.u0) == 120
    @test all(b -> b isa Reflecting, boundaries(prob.space))
    @test_throws ArgumentError problem(:pine_wilt_casim; cover=fill(3, 4, 4))
    @test_throws ArgumentError problem(:pine_wilt_casim; n=20)
end

@testset "grapevine_leafroll" begin
    n = 4000
    model = build_grapevine_leafroll()
    p = defaults(model)
    # a healthy vine with one infective neighbour along its trellis row and one across rows
    state = zeros(Int8, 3, 3)
    state[1, 2] = 3          # along the row (first axis)
    state[2, 3] = 4          # across rows
    hits = count(s -> solve(CAProblem(model, (; state), 1, (; α=0.2); space=Grid(3, 3),
                                      seed=s)).u[end].state[2, 2] ==
                      1, 1:n)
    @test close_to(hits / n, 1 - exp(-(0.2 + p.β)), n)
    # exposed-detectable to infective-moderate at the rate of the vine's age class
    state = fill(Int8(2), 1, 1)
    for (age, L) in ((10, p.L_young), (100, p.L_mature), (300, p.L_old))
        moved = count(s -> solve(CAProblem(model, (; state, age=fill(Int16(age), 1, 1)), 1;
                                           space=Grid(1, 1), seed=s)).u[end].state[1, 1] ==
                           3,
                      1:n)
        @test close_to(moved / n, 1 - exp(-1 / L), n)
    end
    # the undetectable period: a drawn triangular duration, or the source's monthly cdf
    @test EcologicalCellularAutomata._triangular(0.0, 4, 12, 18) == 4 &&
          EcologicalCellularAutomata._triangular(1.0, 4, 12, 18) == 18
    @test EcologicalCellularAutomata._triangular_cdf(12, 4, 12, 18) ≈ 8 / 14
    sol = solve(CAProblem(model, (; state=fill(Int8(1), 1, 1), clock=fill(Int16(3), 1, 1)),
                          4;
                          space=Grid(1, 1)))
    @test [u.state[1, 1] for u in sol.u] == [1, 1, 1, 2, 2]

    # management: roguing the trigger, testing its neighbours, protecting young replants
    u = (state=Int8[0 2 0; 1 3 0; 0 2 0], age=fill(Int16(100), 3, 3),
         clock=fill(Int16(2), 3, 3), rogued=falses(3, 3) |> collect,
         tested=falses(3, 3) |> collect)
    @test grapevine_manage!(deepcopy(u), :Im) == (1, 0)
    v = deepcopy(u)
    # (1,2) and (3,2) are tested, and both are detectable, so three vines are rogued
    @test grapevine_manage!(v, :ImNS) == (3, 2)
    @test v.state[2, 2] == 0 && v.state[1, 2] == 0 && v.tested[1, 2] && v.tested[3, 2] &&
          !v.tested[2, 1]
    v = deepcopy(u)
    @test grapevine_manage!(v, :ImNSEW) == (3, 4)
    @test v.state[2, 1] == 1                             # undetectable: tests negative
    @test grapevine_manage!(deepcopy(u), :ImNS; delay=3) == (0, 0)
    young = merge(u, (; age=fill(Int16(10), 3, 3)))
    @test grapevine_manage!(deepcopy(young), :Im) == (0, 0)
    @test grapevine_manage!(deepcopy(young), :Im; protect_replants=false) == (1, 0)
    @test grapevine_manage!(deepcopy(u), :ImO) == (0, 0)
    @test grapevine_manage!(deepcopy(u), :ImM_NS) == (3, 2)
    @test_throws ArgumentError grapevine_manage!(deepcopy(u), :roguing)

    # the economics: without disease, every vine earns from 36 months
    nodisease = solve(problem(:grapevine_leafroll; infected=0.0, months=100))
    expected = sum(0.9959^t * 5720 * 0.4266 for t in 36:99)
    @test grapevine_npv(nodisease) ≈ expected
    @test grapevine_half_life(nodisease) === nothing
    @test grapevine_revenue((state=Int8[0 4], age=Int16[40 40])) ≈ 0.4266 * (1 + 0.25 * 0.9)
    @test_throws ArgumentError grapevine_npv(solve(problem(:grapevine_leafroll; months=4);
                                                   saveat=2))
    prob = problem(:grapevine_leafroll; initial=:Im, seed=1)
    @test count(==(3), prob.u0.state) == 114 && size(prob.u0.state) == (130, 44)
    @test_throws ArgumentError problem(:grapevine_leafroll; initial=:Ih)
end

@testset "cassava_mealybug" begin
    model = build_cassava_mealybug()
    p = defaults(model)
    # the within-plant equations, one day, against the source's difference equations
    u0 = (state=fill(Int8(1), 1, 1), Pe=fill(10.0, 1, 1), Pi=fill(30.0, 1, 1),
          Pm=fill(5.0, 1, 1), Me=fill(4.0, 1, 1), Mi=fill(1.0, 1, 1), Md=fill(2.0, 1, 1),
          Mm=fill(3.0, 1, 1))
    v = solve(CAProblem(model, u0, 1; space=Grid(1, 1))).u[end]
    @test v.Pi[1] ≈ 30 + p.r1 * p.α1 * 10 - p.α2 * 30 - 20
    @test v.Pm[1] ≈ 5 + p.r2 * p.α2 * 30 - p.α3 * 5 - 5
    @test v.Pe[1] ≈ 10 + p.r3 * p.α4 * p.ν1 * 5 - p.α1 * 10 - 10
    δ1 = 0.0521 * min(1, 45 / 60)
    @test v.Md[1] ≈ 2 + p.s2 * p.γ2 * δ1 * 1 - p.γ3 * 2
    @test v.Me[1] ≈ 4 + p.s4 * p.ν2 * 3 - p.γ1 * 4
    # yield loss: infested two weeks in days 1-120 loses everything
    w = solve(CAProblem(model, merge(u0, (; since=fill(Int16(1), 1, 1))), 20;
                        space=Grid(1, 1))).u[end]
    @test w.loss[1] == 3

    n = 3000
    blank() = (state=zeros(Int8, 7, 7), Pe=zeros(7, 7), Pi=zeros(7, 7), Pm=zeros(7, 7),
               Me=zeros(7, 7), Mi=zeros(7, 7), Md=zeros(7, 7), Mm=zeros(7, 7),
               since=fill(Int16(-1), 7, 7), loss=zeros(Int8, 7, 7))
    q = merge(p, (; n1=0.0, n2=0.0, n3=0.0))
    border = inner = 0
    for s in 1:n
        u = blank()
        cassava_sweep!(u, q, 5, Random.Xoshiro(s))
        border += u.state[1, 1] == 1
        inner += u.state[4, 4] == 1
    end
    @test close_to(border / n, p.w1, n) && close_to(inner / n, p.w2, n)
    # an infested neighbour with an instar passes it on with probability n1 if not yet visited
    gained = 0
    for s in 1:n
        u = blank()
        u.state[4, 5] = 1
        u.Pi[4, 5] = 5.0
        cassava_sweep!(u, merge(p, (; w1=0.0, w2=0.0)), 5, Random.Xoshiro(s))
        gained += u.state[4, 4] == 1
    end
    # the neighbour is visited first in half the sweeps, and then cannot give
    @test close_to(gained / n, p.n1 / 2, n)
    u = blank()
    u.state[4, 5] = 1
    u.Pi[4, 5] = 0.5
    cassava_sweep!(u, merge(p, (; w1=0.0, w2=0.0, n1=1.0)), 5, Random.Xoshiro(1))
    @test u.state[4, 4] == 0
    # a plant whose total is not positive is cleared
    u = blank()
    u.state[4, 4] = 1
    u.Pi[4, 4] = -0.2
    cassava_sweep!(u, merge(p, (; w1=0.0, w2=0.0)), 5, Random.Xoshiro(1))
    @test u.state[4, 4] == 0 && u.Pi[4, 4] == 0
    # survey days: every infested plant sampled (survey of all plants) is removed in days 1-120,
    # and kept in days 121-210 when it has at most m1 mealybugs
    every = merge(p, (; w1=0.0, w2=0.0, n1=0.0, n2=0.0, n3=0.0, survey_plants=49))
    u = blank()
    u.state[4, 4] = 1
    u.Pi[4, 4] = 3.0
    @test cassava_sweep!(u, every, 31, Random.Xoshiro(1)) == (1, 49)
    @test u.state[4, 4] == 2
    u = blank()
    u.state[4, 4] = 1
    u.Pi[4, 4] = 3.0
    @test cassava_sweep!(u, every, 31 + 14 * 7, Random.Xoshiro(1))[1] == 1  # day 129
    @test u.state[4, 4] == 1
    @test cassava_sweep!(u, every, 32, Random.Xoshiro(1)) == (0, 0)
    # releases: G1 when over half the sampled plants are infested, split over infested plants
    u = blank()
    u.state[1, 1] = u.state[2, 2] = 1
    p1, G = cassava_release!(u, merge(p, (; found=40, sampled=60)), 31)
    @test G == p.G1 && u.Mi[1, 1] == p.G1 / 2 && p1.last_release == 31 &&
          p1.released == p.G1
    p2, G = cassava_release!(u, merge(p1, (; found=10)), 45)
    @test G == 0                                     # not due: interval 60
    p3, G = cassava_release!(u, merge(p1, (; found=10)), 91)
    @test G == p.G2
    p4, G = cassava_release!(u, merge(p1, (; found=0)), 91)
    @test G == 0 && p4.last_release == 0
    @test cassava_yield((state=Int8[0 1 2 1], loss=Int8[0 1 0 3])) ≈ 2.25 * 1.9

    prob = problem(:cassava_mealybug; n=20, interval=14, days=40, seed=1)
    integ = init(prob; save_everystep=false)
    sol = solve!(integ)
    @test integ.p.interval == 14 && cassava_yield(sol.u[end]) <= 2.25 * 400
    @test_throws ArgumentError problem(:cassava_mealybug; n=3)
end

@testset "telenomus_release" begin
    model = build_telenomus_release()
    p = defaults(model)
    z(k) = ntuple(_ -> Int16(0), Val(k))
    nomass = StinkBugEggMass(0, 0, 0, 0)
    blank(n) = (nymphs=fill(z(30), n, n), adults=fill(z(36), n, n),
                masses=fill(ntuple(_ -> nomass, Val(60)), n, n), wasps=fill(z(26), n, n))
    put(t, i, k) = Base.setindex(t, Int16(k), i)
    # within-cell step: 400 cells, each with 100 adults aged 20, 100 nymphs aged 10, 100 wasps
    # aged 5, 50 wasps aged 26, an intact egg mass about to hatch and a parasitised one emerging
    u0 = blank(20)
    mass5 = StinkBugEggMass(11, 5, 0, 0)
    mass14 = StinkBugEggMass(0, 0, 11, 14)
    for I in eachindex(u0.adults)
        u0.adults[I] = put(z(36), 20, 100)
        u0.nymphs[I] = put(z(30), 10, 100)
        u0.wasps[I] = put(put(z(26), 5, 100), 26, 50)
        u0.masses[I] = Base.setindex(Base.setindex(u0.masses[I], mass5, 1), mass14, 2)
    end
    v = solve(CAProblem(model, u0, 1; space=Grid(20, 20; boundary=Reflecting()), seed=3)).u[end]
    n = 400 * 100
    # the first day: adults do not age, and survive with probability 1 - 0.03
    @test all(c -> c[21] == 0, v.adults)
    @test close_to(sum(c -> c[20], v.adults) / n, 1 - 0.03, n)
    # nymphs age a day and survive with probability 1 - 0.014
    @test close_to(sum(c -> c[11], v.nymphs) / n, 1 - 0.014, n)
    # wasps die with probability 0.0566, and all at age 26
    @test close_to(sum(c -> c[5], v.wasps) / n, 1 - 0.0566, n) &&
          all(c -> c[26] == 0, v.wasps)
    # an intact mass hatches floor(11 x 0.5 x 0.63) = 3 nymphs; parasitoids are removed at 14
    @test all(c -> c[1] == 3, v.nymphs)
    @test all(c -> c[1] == nomass && c[2] == nomass, v.masses)
    # with offspring, they give floor(11 x 0.7 x 0.8243) = 6 wasps
    w = solve(CAProblem(model, u0, 1, (; offspring=true);
                        space=Grid(20, 20; boundary=Reflecting()), seed=3)).u[end]
    @test all(c -> c[1] == 6, w.wasps)
    # the second day: adults age, and those of 36 die
    u1 = blank(2)
    u1.adults[1] = put(put(z(36), 36, 10), 10, 10)
    v = solve(CAProblem(model, u1, 2; space=Grid(2, 2; boundary=Reflecting()), seed=1)).u[end]
    @test v.adults[1][36] == 0 && v.adults[1][10] == 0 && 0 < v.adults[1][11] <= 10
    # nymphs of 30 become adults while the cell has room for 300
    u1 = blank(2)
    u1.adults[1] = put(z(36), 5, 295)
    u1.nymphs[1] = put(z(30), 30, 20)
    v = solve(CAProblem(model, u1, 1, (; mort_adult=0.0);
                        space=Grid(2, 2; boundary=Reflecting()))).u[end]
    @test v.adults[1][1] == 5 && v.nymphs[1][30] == 0

    # the source's dispersal tables are truncated negative binomials of mean 0.5
    @test all(isapprox.(telenomus_dispersal_steps(0.5, 2.8433), p.nymph_steps; atol=1e-4))
    @test all(isapprox.(telenomus_dispersal_steps(0.5, 3.1728), p.adult_steps; atol=1e-4))

    # one stink bug move: stays with probability P0^2, moves one cell right with (P1 - P0) / 2 ...
    q = merge(p, (; colonisation_days=0, p_oviposit=1.0))
    runs = 4000
    stay = right = laid = 0
    for s in 1:runs
        u = blank(21)
        u.adults[11, 11] = put(z(36), 5, 1)
        telenomus_step!(u, q, 1, Random.Xoshiro(s))
        stay += u.adults[11, 11][5] == 1
        right += u.adults[12, 11][5] == 1 && u.adults[11, 11][5] == 0
    end
    P = p.adult_steps
    @test close_to(stay / runs, P[1]^2, runs)
    @test close_to(right / runs, (P[2] - P[1]) / 2 * P[1], runs)
    # ... and a mature female lays an egg mass of 11 with probability p_oviposit
    q = merge(p, (; colonisation_days=0))
    for s in 1:runs
        u = blank(21)
        u.adults[11, 11] = put(z(36), 20, 1)
        telenomus_step!(u, q, 1, Random.Xoshiro(s))
        laid += count(c -> c[1] == StinkBugEggMass(11, 1, 0, 0), u.masses)
    end
    @test close_to(laid / runs, p.p_oviposit, runs)
    # nymphs up to age 5 stay
    u = blank(21)
    u.nymphs[11, 11] = put(z(30), 5, 40)
    telenomus_step!(u, q, 1, Random.Xoshiro(1))
    @test u.nymphs[11, 11][5] == 40

    # a wasp of age 2 flies to the only egg mass within 8 m and lays 12 eggs: 11 in it, and one
    # in a mass next to it; then it ages a day
    u = blank(21)
    u.wasps[11, 11] = put(z(26), 2, 1)
    u.masses[16, 11] = Base.setindex(u.masses[16, 11], StinkBugEggMass(11, 2, 0, 0), 1)
    u.masses[17, 11] = Base.setindex(u.masses[17, 11], StinkBugEggMass(5, 3, 2, 4), 1)
    q2 = merge(q, (; radius=5))
    telenomus_step!(u, q2, 1, Random.Xoshiro(1))
    @test u.wasps[16, 11][3] == 1 && u.wasps[11, 11] == z(26)
    @test u.masses[16, 11][1] == StinkBugEggMass(0, 0, 11, 1)
    # the mean age of 2 immatures of age 4 and 1 new one is 3
    @test u.masses[17, 11][1] == StinkBugEggMass(4, 3, 3, 3)
    # the source's fecundity: 12 eggs on the second day, 1 on the 25th
    @test floor(Int, p.fecundity * exp(-p.fecundity_decay)) == 12
    @test floor(Int, p.fecundity * exp(-p.fecundity_decay * 24)) == 1
    # with no eggs within reach, wasps move within the disc of radius 8
    far = 0
    for s in 1:200
        u = blank(41)
        u.wasps[21, 21] = put(z(26), 3, 1)
        telenomus_step!(u, q, 1, Random.Xoshiro(s))
        I = findfirst(c -> c[4] == 1, u.wasps)
        far = max(far, isqrt((I[1] - 21)^2 + (I[2] - 21)^2))
    end
    @test far == 8

    # releases: 4 points 25 m apart in the central hectare, or 4 strips of 100 cells
    pts = telenomus_release_cells(merge(p, (; mode=:points)), 200, 200)
    @test sort(unique(first.(pts))) == [63, 88, 113, 138] == sort(unique(last.(pts)))
    strips = telenomus_release_cells(p, 200, 200)
    @test length(strips) == 400 && sort(unique(last.(strips))) == 51:150
    for (mode, per_ha, spacing) in
        ((:points, 5000, 25), (:points, 10000, 50), (:strips, 3000, 5))
        u = blank(200)
        telenomus_step!(u, merge(q, (; releases=1, mode, per_ha, spacing)), 4,
                        Random.Xoshiro(1))
        @test sum(c -> c[1], u.wasps) == per_ha && maximum(c -> c[1], u.wasps) <= 550
    end
    # colonists arrive in the central 81 x 81 m, aged 1 to 13
    u = blank(200)
    telenomus_step!(u, p, 3, Random.Xoshiro(1))
    @test sum(c -> sum(Int, c), u.adults) == 800
    occupied = findall(c -> sum(c) > 0, u.adults)
    @test all(I -> all(60 .<= Tuple(I) .<= 140), occupied)
    @test all(c -> sum(c[14:end]) == 0, u.adults)

    # a small field runs, with releases, and is reproducible
    small = problem(:telenomus_release; n=40, releases=2, mode=:points, spacing=10,
                    pmap=(; release_area=20, colonised=10, colonists=400), days=20, seed=4)
    s1 = solve(small; saveat=20)
    s2 = solve(small; saveat=20)
    @test s1.u[end] == s2.u[end] && telenomus_counts(s1.u[end]).adults > 0
    @test telenomus_et(s1.u[end]; area=20) > 0
    @test_throws ArgumentError problem(:telenomus_release; mode=:rows)
    @test_throws ArgumentError problem(:telenomus_release; n=50)
    @test_throws ArgumentError telenomus_dispersal_steps(0.5, 0)
end

@testset "bark_beetle_forest" begin
    model = build_bark_beetle_forest()
    p = defaults(model)
    blank(n1, n2=n1) = (tree=fill(false, n1, n2), age=zeros(Int16, n1, n2),
                        mass=zeros(n1, n2),
                        beetles=zeros(n1, n2), fate=zeros(Int8, n1, n2),
                        contaminated=fill(false, n1, n2))
    function plant!(u, I, age)
        u.tree[I...] = true
        u.age[I...] = age
        u.mass[I...] = 28 * min(age, 70) + 14 * clamp(age - 70, 0, 70)
        return u
    end
    rulestep(u, pmap=(;); seed=1) = solve(CAProblem(model, u, 1, pmap;
                                                    space=Grid(size(u.tree)...),
                                                    seed)).u[end]
    # the yearly step of a tree: three blocks of 1800 trees aged 2, 50 and 151, which die with
    # probability 0.05 (reaching 3), 0.01 and 0.2 (beyond 150)
    u = blank(90, 60)
    for I in CartesianIndices(u.tree)
        plant!(u, Tuple(I), I[1] <= 30 ? 2 : I[1] <= 60 ? 50 : 151)
    end
    v = rulestep(u)
    alive(rows) = count(v.tree[rows, :]) / (30 * 60)
    @test close_to(alive(1:30), 1 - 0.05, 1800)
    @test close_to(alive(31:60), 1 - 0.01, 1800)
    @test close_to(alive(61:90), 1 - 0.2, 1800)
    @test all(I -> !v.tree[I] || v.age[I] == u.age[I] + 1, eachindex(v.tree))
    # in a drought, mortality is multiplied by 1.2
    v = rulestep(u, (; droughts=(0,)); seed=2)
    @test close_to(count(v.tree[31:60, :]) / 1800, 1 - 0.012, 1800)
    # growth: 28 kg a year until 70, 14 until 140, nothing after; halved in a drought
    u = blank(1, 4)
    foreach(((k, a),) -> plant!(u, (1, k), a), enumerate((69, 70, 139, 140)))
    q = (; mort_adult=0.0)
    v = rulestep(u, q)
    @test v.mass[1, :] - u.mass[1, :] == [28, 14, 14, 0]
    v = rulestep(u, merge(q, (; droughts=(0,))))
    @test v.mass[1, :] - u.mass[1, :] == [14, 7, 7, 0]
    # contaminated ground kills trees older than 3; felled and killed trees are removed
    u = blank(1, 4)
    foreach(((k, a),) -> plant!(u, (1, k), a), enumerate((2, 3, 10, 10)))
    u.contaminated[1, 1:2] .= true
    u.fate[1, 3], u.fate[1, 4] = 1, 2
    u.beetles[1, 4] = 150
    v = rulestep(u, (; mort_sapling=0.0, mort_adult=0.0))
    @test v.tree[1, :] == [true, false, false, false] && v.age[1, 1] == 3
    @test v.beetles[1, 4] == 0 && all(==(0), v.fate)

    # spreading: a 6 x 6 lattice of trees 7 m apart; a tree with n others within the square of
    # half-width 10 spreads with probability 0 (n = 0), 1/2 (1-3), 1 (4-9), 1/2 (10-50), 0 (>50).
    # A sapling lands uniformly in the 13 x 13 square around its parent, failing only on the
    # parent itself (others are 7 m away; spacing 0 and no site limit)
    ps(n) = n == 0 ? 0.0 : n <= 3 ? 0.5 : n <= 9 ? 1.0 : n <= 50 ? 0.5 : 0.0
    pos = 8:7:43
    trees = [(i, j) for i in pos for j in pos]
    nbrs(t) = count(s -> s != t && maximum(abs.(s .- t)) <= 10, trees)
    expected = sum(t -> ps(nbrs(t)), trees) * 168 / 169
    q = merge(p, (; spacing=0.0, site_limit=-1))
    runs = 300
    total = 0
    for s in 1:runs
        u = blank(50)
        foreach(t -> plant!(u, t, 20), trees)
        bark_beetle_step!(u, q, 2, Random.Xoshiro(s))
        total += count(u.tree) - length(trees)
    end
    @test sort(unique(nbrs.(trees))) == [3, 5, 8]
    @test abs(total / runs - expected) <= 4 * sqrt(expected / runs) + 0.1
    # nothing spreads in odd years, or from an isolated tree, or where the site limit is 0
    u = blank(50)
    foreach(t -> plant!(u, t, 20), trees)
    @test count(bark_beetle_step!(deepcopy(u), q, 3, Random.Xoshiro(1)).tree) == 36
    @test count(bark_beetle_step!(deepcopy(u), merge(q, (; site_limit=0)), 2,
                                  Random.Xoshiro(1)).tree) == 36
    lone = plant!(blank(30), (15, 15), 20)
    @test count(bark_beetle_step!(lone, q, 2, Random.Xoshiro(1)).tree) == 1
    # a sapling never takes root within 2 m of another tree: only the initial trees, aged 70 or
    # more after 60 years, may stand closer together
    sol = solve(problem(:bark_beetle_forest; n=80, years=60, seed=3); saveat=60)
    w = sol.u[end]
    T = findall(w.tree)
    @test length(T) > 100
    @test all(((a, b),) -> a == b || (w.age[a] >= 70 && w.age[b] >= 70) ||
                           sum(abs2, Tuple(a) .- Tuple(b)) > 4, Iterators.product(T, T))

    # harvest: in year 101 trees older than 90 in rows 1-30 of 300 are felled, with their
    # beetles; in year 102 rows 31-60; none after year 250
    q = merge(p, (; harvest=true))
    u = blank(300, 2)
    for i in (1, 30, 31, 60)
        plant!(u, (i, 1), 91)
        plant!(u, (i, 2), 90)
    end
    u.beetles[1, 1] = 50
    bark_beetle_step!(u, q, 101, Random.Xoshiro(1))
    @test u.fate[[1, 30, 31, 60], 1] == [1, 1, 0, 0] && all(==(0), u.fate[:, 2])
    @test u.beetles[1, 1] == 0
    bark_beetle_step!(u, q, 102, Random.Xoshiro(1))
    @test u.fate[[31, 60], 1] == [1, 1]
    u = plant!(blank(300, 1), (1, 1), 100)
    @test bark_beetle_step!(u, q, 251, Random.Xoshiro(1)).fate[1, 1] == 0

    # infestation: 100 draws of 20 beetles on the trees in rows and columns 100-120
    q = merge(p, (; infestation=5, beetle_spreads=0, beetle_kill=Inf))
    u = blank(130)
    foreach(t -> plant!(u, t, 20),
            ((100, 100), (110, 115), (120, 120), (121, 110), (99, 105)))
    bark_beetle_step!(u, q, 5, Random.Xoshiro(1))
    @test sum(u.beetles) == 2000 && u.beetles[121, 110] == 0 && u.beetles[99, 105] == 0
    @test count(>(0), u.beetles) == 3 && all(x -> x % 20 == 0, u.beetles)
    u = plant!(blank(130), (10, 10), 20)
    @test sum(bark_beetle_step!(u, q, 5, Random.Xoshiro(1)).beetles) == 0

    # one spread: an infested tree multiplies its beetles by 10 and shares the same number
    # equally among the 4 other trees within the square of half-width 10, rounded down to whole
    # beetles; a tree 11 m away gets none. Trees reaching 200 die with their beetles
    u = blank(41)
    foreach(t -> plant!(u, t, 20),
            ((21, 21), (31, 31), (11, 21), (21, 30), (25, 25), (21, 32)))
    u.beetles[21, 21] = 3
    q = merge(p, (; beetle_spreads=1))
    v = bark_beetle_step!(deepcopy(u), q, 1, Random.Xoshiro(1))
    @test v.beetles[21, 21] == 30
    @test [v.beetles[31, 31], v.beetles[11, 21], v.beetles[21, 30], v.beetles[25, 25]] ==
          fill(7.0, 4)
    @test v.beetles[21, 32] == 0
    v = bark_beetle_step!(deepcopy(u), merge(q, (; whole_beetles=false)), 1,
                          Random.Xoshiro(1))
    @test v.beetles[25, 25] == 7.5
    u.beetles[21, 21] = 20
    v = bark_beetle_step!(deepcopy(u), q, 1, Random.Xoshiro(1))
    @test v.fate[21, 21] == 2 && v.beetles[21, 21] == 0 && v.beetles[25, 25] == 50
    # three spreads a year: a lone tree with one beetle has 1000 by the end of the year, and dies
    u = plant!(blank(30), (15, 15), 20)
    u.beetles[15, 15] = 1
    v = bark_beetle_step!(u, merge(p, (; beetle_kill=Inf)), 1, Random.Xoshiro(1))
    @test v.beetles[15, 15] == 1000
    u.beetles[15, 15] = 1
    @test bark_beetle_step!(u, p, 1, Random.Xoshiro(1)).fate[15, 15] == 2

    # the setup: 20 distinct trees aged 10 (280 kg) within 20 m of the centre; contamination of
    # rows 90-120 and columns 120-180
    prob = problem(:bark_beetle_forest; seed=5, contamination=true)
    u0 = prob.u0
    T = findall(u0.tree)
    @test length(T) == 20 && all(I -> u0.age[I] == 10 && u0.mass[I] == 280, T)
    @test all(I -> (I[1] - 150.5)^2 + (I[2] - 150.5)^2 <= 400, T)
    @test findall(u0.contaminated) == vec(CartesianIndices((90:120, 120:180)))
    # counts, growth rate, reproducibility and errors
    u = blank(3)
    plant!(u, (1, 1), 10)
    plant!(u, (2, 2), 100)
    plant!(u, (3, 3), 100)
    u.fate[3, 3] = 1
    u.beetles[1, 1] = 5
    @test bark_beetle_counts(u) ==
          (trees=2, mass=0.28 + 2.38, felled=2.38, killed=0, infested=1, beetles=5.0)
    @test bark_beetle_growth_rate([1.02^t for t in 0:300]) ≈ 2.0
    @test_throws ArgumentError bark_beetle_growth_rate(ones(100))
    small = problem(:bark_beetle_forest; n=60, years=40, infestation=30, seed=2,
                    pmap=(; infest_rows=(20, 40), infest_cols=(20, 40)))
    @test solve(small; saveat=40).u[end] == solve(small; saveat=40).u[end]
    @test_throws ArgumentError problem(:bark_beetle_forest; trees=5000)
    @test_throws ArgumentError problem(:bark_beetle_forest; droughts=(1.5,))
    @test_throws ArgumentError problem(:bark_beetle_forest; n=0)
end

@testset "potato_moth_farmers" begin
    model = build_potato_moth_farmers()
    p = defaults(model)
    # the source's Sharpe–DeMichele survival of eggs, larvae and pupae (its Table 3, R = 1.987)
    sdm(T, a, b, c, d, e, f) = a * T / 298.16 * exp(b * (1 / 298.16 - 1 / T)) /
                               (1 + exp(c / 1.987 * (1 / d - 1 / T)) +
                                exp(e / 1.987 * (1 / f - 1 / T)))
    S(t) = sdm(t + 273.15, 0.822, -758.5, -212100, 281.9, 405200, 303.8) *
           sdm(t + 273.15, 0.758, -180.2, -475700, 282.7, 1298000, 301.5) *
           sdm(t + 273.15, 0.9, -73.72, -1263000, 286.5, 1095000, 306.3)
    @test all(t -> potato_moth_survival(t, p) ≈ S(t), (5.0, 12.0, 13.0, 18.0, 23.0))
    @test potato_moth_survival(18.0, p) > 0.5 && potato_moth_survival(12.0, p) < 1e-4

    land = potato_moth_landscape()
    @test size(land.temperature) == (80, 80) && land.temperature[50, 41] == 12
    @test count(>(0), land.population) == 37 && land.population[50, 41] == 3000
    blank(n1, n2, T) = (immature=zeros(n1, n2), source=zeros(Int8, n1, n2),
                        temperature_now=fill(Float64(T), n1, n2),
                        temperature=fill(Float64(T), n1, n2),
                        population=zeros(Int16, n1, n2),
                        farmers=zeros(UInt64, n1, n2), carriers=zeros(UInt64, n1, n2))
    q = (; size=(9, 9), temperature_sd=0.0, farmers=0)
    step1(u, q) = solve(CAProblem(model, u, 1, q; space=Grid(size(u.immature)...),
                                  seed=1)).u[end]
    # one generation at 18 °C: survivors (1 - crude) S I with crude = 1 - exp(-0.205 x 2), half
    # gravid, 10% killed; 30% of them leave, 7.5% to each neighbour (shares off the lattice stay);
    # then logistic growth K G / (G + (K - G) e^-5)
    I0 = 400.0
    G0 = 0.45 * I0 * exp(-0.41) * S(18.0)
    logistic(G) = 1000G / (G + (1000 - G) * exp(-5))
    u = blank(9, 9, 18)
    u.immature[5, 5] = I0
    v = step1(u, q)
    @test v.immature[5, 5] ≈ logistic(0.7G0)
    @test v.immature[4, 5] ≈ v.immature[5, 6] ≈ logistic(0.075G0)
    @test count(>(0), v.immature) == 5 && v.source[4, 5] == 0
    u = blank(9, 9, 18)
    u.immature[1, 1] = I0
    @test step1(u, q).immature[1, 1] ≈ logistic((0.7 + 2 * 0.075) * G0)
    @test step1(u, merge(q, (; neighbourhood=false))).immature[1, 1] ≈ logistic(G0)
    # with no moths left after mortality, a cell takes the source of a random neighbour with
    # gravid females: between a red and a blue one, each with probability 1/2
    u = blank(3, 400, 18)
    u.immature[1, :] .= I0
    u.immature[3, :] .= I0
    u.source[1, :] .= 2
    u.source[3, :] .= 3
    v = step1(u, merge(q, (; size=(3, 400))))
    @test close_to(count(==(2), v.source[2, :]) / 400, 0.5, 400)
    # the next generation's temperature is normal about the mean with sd 0.1
    v = step1(blank(80, 80, 15), merge(q, (; size=(80, 80), temperature_sd=0.1)))
    m = sum(v.temperature_now) / 6400
    @test abs(m - 15) < 4 * 0.1 / 80
    @test abs(sqrt(sum(abs2, v.temperature_now .- m) / 6399) - 0.1) < 0.005
    # without farmers, moths reach in t generations every cell within t steps of the first
    prob = problem(:potato_moth_farmers; farmers=0, generations=5)
    w = solve(prob; saveat=5).u[end]
    @test findall(>(0), w.immature) ==
          [I for I in CartesianIndices((80, 80)) if sum(abs.(Tuple(I) .- (50, 41))) <= 5]
    @test potato_moth_counts(w).infested == 61 == potato_moth_counts(w).biological

    # farmers: two, on the same moth-free village cell; farmer 1 carries moths
    put(u, I, f; carrier=false) = (u.farmers[I...] |= UInt64(1) << (f - 1);
                                   carrier && (u.carriers[I...] |= UInt64(1) << (f - 1));
                                   u)
    fq = merge(p, (; farmers=2, moves=1, knowledge=(30.0, 80.0), size=(80, 80)))
    function village_state()
        u = blank(80, 80, 15)
        u.population .= Int16.(land.population)
        return u
    end
    u = put(put(village_state(), (65, 19), 1; carrier=true), (65, 19), 2)
    potato_moth_farmers_step!(u, fq, Random.Xoshiro(1))
    # the carrier left 90 immatures marked red; farmer 2 (knowledge 80 < 100) then picked them
    # up, if it acted after farmer 1; each then left for another village
    @test u.immature[65, 19] == 90 && u.source[65, 19] == 2
    @test sum(count_ones, u.farmers) == 2 && u.farmers[65, 19] == 0
    @test all(I -> land.population[I] > 0, findall(!=(0), u.farmers))
    # under the threshold rule, only farmers with knowledge below it pick up moths
    u = put(put(village_state(), (65, 19), 1), (65, 19), 2)
    u.immature[65, 19] = 10.0
    potato_moth_farmers_step!(u, merge(fq, (; threshold=50.0)), Random.Xoshiro(1))
    @test sum(count_ones, u.carriers) == 1 && u.carriers[findfirst(!=(0), u.carriers)] == 1
    # under the probability rule, a farmer of knowledge 30 picks them up with probability 0.7
    runs = 2000
    caught = 0
    one = merge(fq, (; farmers=1, knowledge=(30.0,)))
    for s in 1:runs
        u = put(village_state(), (65, 19), 1)
        u.immature[65, 19] = 10.0
        potato_moth_farmers_step!(u, merge(one, (; infection=:probability)),
                                  Random.Xoshiro(s))
        caught += sum(count_ones, u.carriers)
    end
    @test close_to(caught / runs, 0.7, runs)
    # a farmer at Simiatug, where there are moths, goes (with probability p_nearest = 1) to the
    # village j minimising population_j x distance_j; from a cell without moths, anywhere
    villages = [I for I in CartesianIndices((80, 80)) if land.population[I] > 0]
    others = filter(!=(CartesianIndex(50, 41)), villages)
    target = argmin(J -> land.population[J] * hypot((Tuple(J) .- (50, 41))...), others)
    nq = merge(one, (; p_nearest=1.0))
    u = put(village_state(), (50, 41), 1; carrier=true)
    u.immature[50, 41] = 10.0
    potato_moth_farmers_step!(u, nq, Random.Xoshiro(1))
    @test findfirst(!=(0), u.farmers) == target
    hits = 0
    for s in 1:runs
        u = put(village_state(), (1, 1), 1)
        potato_moth_farmers_step!(u, nq, Random.Xoshiro(s))
        hits += u.farmers[50, 41] != 0
    end
    @test close_to(hits / runs, 1 / 37, runs)
    @test_throws ArgumentError potato_moth_farmers_step!(put(village_state(), (1, 1), 2),
                                                         fq,
                                                         Random.Xoshiro(1))

    # the setup: 90 immatures at Simiatug; farmers placed, carrying nothing; knowledge drawn
    prob = problem(:potato_moth_farmers; farmers=20, seed=3)
    @test prob.u0.immature[50, 41] == 90 && count(>(0), prob.u0.immature) == 1
    @test sum(count_ones, prob.u0.farmers) == 20 && all(==(0), prob.u0.carriers)
    @test all(κ -> κ in 0:99, prob.p.knowledge)
    κ = problem(:potato_moth_farmers; farmers=4, knowledge=(50, 200), seed=1).p.knowledge
    @test all(x -> 0 <= x <= 100, κ) && any(x -> x in (0.0, 100.0), κ)
    @test problem(:potato_moth_farmers; farmers=2, knowledge=[10, 20]).p.knowledge ==
          (10.0, 20.0)
    s1 = solve(problem(:potato_moth_farmers; generations=8, seed=2); saveat=8).u[end]
    s2 = solve(problem(:potato_moth_farmers; generations=8, seed=2); saveat=8).u[end]
    @test s1 == s2 && potato_moth_counts(s1).carriers > 0
    @test potato_moth_counts(s1; threshold=1).infested <= potato_moth_counts(s1).infested
    @test potato_moth_invasion_time(problem(:potato_moth_farmers; farmers=0)) == 13
    @test potato_moth_invasion_time(problem(:potato_moth_farmers; farmers=0, generations=5)) ===
          nothing
    @test_throws ArgumentError problem(:potato_moth_farmers; farmers=3)
    @test_throws ArgumentError problem(:potato_moth_farmers; infection=:always)
    @test_throws ArgumentError problem(:potato_moth_farmers; farmers=2, knowledge=[1, 2, 3])
    @test_throws ArgumentError potato_moth_invasion_time(problem(:potato_moth_farmers);
                                                         fraction=0)
end

@testset "mpb_multiscale" begin
    model = build_mpb_multiscale()
    p = defaults(model)
    blank(n) = (density=zeros(n, n), blocked=fill(false, n, n))
    # hexagons: odd columns sit lower. The six neighbours of an odd-column cell (5, 5) and of an
    # even-column one (5, 6), as in the source's code
    hexes(i, j) = isodd(j) ?
                  [(i - 1, j - 1), (i - 1, j), (i - 1, j + 1), (i, j - 1), (i, j + 1),
                   (i + 1, j)] :
                  [(i - 1, j), (i, j - 1), (i, j + 1), (i + 1, j - 1), (i + 1, j + 1),
                   (i + 1, j)]
    # complexity (6N - Σ n_a)/N: 6 for one cell, 5 for two neighbours, 18/7 for a cell and its ring
    u = blank(12)
    u.density[5, 5] = 0.5
    @test mpb_complexity(u) == 6
    u.density[4, 5] = 0.5
    @test mpb_complexity(u) == 5
    u = blank(12)
    for c in [(5, 5); hexes(5, 5)]
        u.density[c...] = 0.5
    end
    @test mpb_complexity(u) ≈ 18 / 7
    @test mpb_radius(u) ≈ sqrt(6 * sqrt(3) * 7 / π)
    @test isnan(mpb_complexity(blank(5)))

    # the microscale: from x beetles, one step of length 1/N gains one with probability
    # b x (1 - x/C)/N and loses one with d x/N; with N = 1, M = 1, b = 0.2, d = 0.1, x = 2
    # (π = 0.02, C = 100): 3 with 0.392, 1 with 0.2
    q = (; birth=0.2, death=0.1, micro_steps=1, replicates=1)
    u = blank(60)
    u.density .= 0.02
    v = solve(CAProblem(model, u, 1, q; space=Grid(60, 60), seed=1)).u[end]
    @test close_to(count(≈(0.03), v.density) / 3600, 0.2 * 2 * 0.98, 3600)
    @test close_to(count(≈(0.01), v.density) / 3600, 0.1 * 2, 3600)
    # fewer than one beetle die out; M runs are averaged
    u.density .= 0.005
    @test all(==(0), solve(CAProblem(model, u, 1, q; space=Grid(60, 60))).u[end].density)
    u.density .= 0.02
    v = solve(CAProblem(model, u, 1, merge(q, (; replicates=2)); space=Grid(60, 60))).u[end]
    @test all(x -> any(y -> x ≈ y, 0.01:0.005:0.03), v.density)
    # periodic birth: b (1 + a sin(2π t/T)), here at t = 1 with T = 4, so b (1 + 2/3)
    v = solve(CAProblem(model, u, 1, merge(q, (; schedule=:periodic, period=4));
                        space=Grid(60, 60), seed=2)).u[end]
    @test close_to(count(≈(0.03), v.density) / 3600, 0.2 * 5 / 3 * 2 * 0.98, 3600)

    # the macroscale. With beta near 0 every cell with a density difference is drawn. A lone cell
    # of density π and its six neighbours: the mean is π/7; the cell contributes
    # γ1 (π - π/7); each neighbour receives γ2 γ1 (6π/7)/7, and the cell keeps π - γ1 (6π/7) plus
    # γ1 (6π/7)/7 (the source's code) or γ2 times that (its article)
    π0, γ1, γ2 = 0.7, 0.6, 0.5
    sure = merge(p, (; beta=1e-12, gamma1=γ1, gamma2=γ2, phi=1.0, rho=0.0))
    u = blank(12)
    u.density[5, 5] = π0
    mpb_disperse!(u, sure, Random.Xoshiro(1))
    pool = γ1 * 6π0 / 7
    @test u.density[5, 5] ≈ π0 - pool + pool / 7
    @test all(c -> u.density[c...] ≈ γ2 * pool / 7, hexes(5, 5))
    @test count(>(0), u.density) == 7
    u = blank(12)
    u.density[5, 6] = π0
    mpb_disperse!(u, merge(sure, (; gamma2_contributors=true)), Random.Xoshiro(1))
    @test u.density[5, 6] ≈ π0 - pool + γ2 * pool / 7
    @test all(c -> u.density[c...] ≈ γ2 * pool / 7, hexes(5, 6))
    # a receiving cell keeps its beetles with probability φ + ρ atan(π)
    runs = 600
    kept = 0
    for s in 1:runs
        u = blank(12)
        u.density[5, 5] = π0
        mpb_disperse!(u, merge(sure, (; phi=0.2, rho=0.5)), Random.Xoshiro(s))
        kept += count(c -> u.density[c...] > 0, hexes(5, 5))
    end
    @test close_to(kept / (6runs), 0.2 + 0.5 * atan(γ2 * pool / 7), 6runs)
    # inclusion: a lone cell and each neighbour are drawn with probability w d/(1 + w) = π/2
    # (w = 1, d = π); the cell changes only if drawn with a neighbour
    changed = 0
    for s in 1:runs
        u = blank(12)
        u.density[5, 5] = 0.8
        mpb_disperse!(u, merge(p, (; phi=1.0, rho=0.0)), Random.Xoshiro(s))
        changed += u.density[5, 5] != 0.8
    end
    @test close_to(changed / runs, 0.4 * (1 - 0.6^6), runs)
    # absorbing edges: what reaches cells beyond the lattice is lost
    u = blank(12)
    u.density[1, 1] = π0
    mpb_disperse!(u, sure, Random.Xoshiro(1))
    @test findall(>(0), u.density) == CartesianIndex.([(1, 1), (2, 1), (1, 2)])
    @test sum(u.density) ≈ π0 - pool + pool / 7 + 2 * γ2 * pool / 7
    # blocked cells take no part
    u = blank(12)
    u.density[5, 5] = π0
    u.blocked[4, 5] = true
    mpb_disperse!(u, sure, Random.Xoshiro(1))
    @test u.density[4, 5] == 0 && count(>(0), u.density) == 6
    # with linkage 2, distant clusters disperse independently; with one pool (linkage 0), not
    function far(linkage, both)
        u = blank(20)
        u.density[5, 5] = π0
        both && (u.density[15, 15] = 0.3)
        mpb_disperse!(u, merge(sure, (; linkage)), Random.Xoshiro(1))
        return u.density[5, 5]
    end
    @test far(2, true) == far(2, false) && far(0, true) != far(0, false)
    # wind: low-density cells upwind of the neighbourhood's centre are rarely drawn
    upwind = 0
    for s in 1:runs
        u = blank(12)
        u.density[5, 5] = π0
        mpb_disperse!(u, merge(sure, (; wind=:east, wind_factor=0.0)), Random.Xoshiro(s))
        upwind += u.density[5, 4] > 0 || u.density[4, 4] > 0
    end
    @test upwind == 0

    # the setup and the source's patterns
    prob = problem(:mpb_multiscale)
    @test size(prob.u0.density) == (50, 50) && count(>(0), prob.u0.density) == 45
    @test prob.u0.density[21, 25] == 0.2 && prob.u0.density[26, 22] == 0.7
    @test count(>(0), mpb_pattern(:sweep)) == 37 && count(>(0), mpb_pattern(:endings)) == 53
    @test count(>(0), mpb_pattern(:behaviour; n=30)) == 45
    mask = fill(false, 50, 50)
    mask[25, 25] = true
    @test problem(:mpb_multiscale; blocked=mask).u0.density[25, 25] == 0
    small = problem(:mpb_multiscale; n=30, steps=4, seed=3)
    s1 = solve(small; saveat=4).u[end]
    @test s1 == solve(small; saveat=4).u[end] && count(>(0), s1.density) > 0
    @test all(x -> 0 <= x <= 1, s1.density)
    @test_throws ArgumentError problem(:mpb_multiscale; pattern=:spiral)
    @test_throws ArgumentError problem(:mpb_multiscale; pmap=(; linkage=1))
    @test_throws ArgumentError problem(:mpb_multiscale; pmap=(; schedule=:yearly))
    @test_throws ArgumentError problem(:mpb_multiscale; pattern=fill(2.0, 50, 50))
end

@testset "western_x_orchard" begin
    model = build_western_x_orchard()
    p = defaults(model)
    orchard(S) = (state=Int8.(S), age=fill(Int16(60), size(S)), clock=zeros(Int16, size(S)))
    month(u, q; seed=1) = solve(CAProblem(model, u, 1, q; space=Grid(size(u.state)...),
                                          seed)).u[end]
    # inoculation: an unvisited tree with k visible or high neighbours, with probability
    # alpha + beta k a year, so 1 - (1 - alpha - beta k)^(1/12) a month; here k = 3 along a
    # row below a row of visible trees, and k = 0 in the row below that
    S = fill(1, 3, 2000)
    S[1, :] .= 4
    v = month(orchard(S), p)
    q3 = 1 - (1 - (0.015 + 3 * 0.2))^(1 / 12)
    @test close_to(count(==(2), v.state[2, 2:1999]) / 1998, q3, 1998)
    @test close_to(count(==(2), v.state[3, :]) / 2000, 1 - (1 - 0.015)^(1 / 12), 2000)
    v = month(orchard(S), merge(p, (; rates=:monthly)); seed=2)
    @test close_to(count(==(2), v.state[2, 2:1999]) / 1998, 0.615, 1998)
    # progression: after a whole year in a state, the next with probability 1/2; not before;
    # high is final; an empty cell stays empty
    u = orchard(fill(2, 40, 40))
    u.clock .= 11
    v = month(u, p)
    @test close_to(count(==(3), v.state) / 1600, 0.5, 1600) &&
          all(==(0), v.clock[v.state .== 3])
    u.clock .= 5
    @test all(==(2), month(u, p).state)
    u = orchard(fill(5, 4, 4))
    u.clock .= 11
    @test all(==(5), month(u, p).state)
    @test all(==(0), month(orchard(zeros(Int, 3, 3)), p).state)
    @test month(orchard(fill(1, 3, 3)), merge(p, (; alpha=0.0))).age[2, 2] == 61

    # management: with 6 of 600 trees symptomatic (1%), the threshold rule cuts them all in a
    # decision month, and replants them as unvisited trees of age 0; with 5 it does nothing
    function managed(nsymptomatic, strategy; month=24, kw...)
        u = (state=fill(Int8(1), 20, 30), age=fill(Int16(60), 20, 30),
             clock=zeros(Int16, 20, 30), cut=falses(20, 30) |> collect,
             replanted=falses(20, 30) |> collect, tested=falses(20, 30) |> collect)
        u.state[1:nsymptomatic] .= 4
        q = merge(p, (; strategy, interval=24, first_decision=24, kw...))
        return western_x_manage!(u, q, month, Random.Xoshiro(1))
    end
    v = managed(6, :threshold_cut)
    @test count(v.cut) == 6 && all(==(0), v.state[1:6]) && !any(v.replanted)
    v = managed(6, :threshold_replant)
    @test all(==(1), v.state[1:6]) && all(==(0), v.age[1:6]) && count(v.replanted) == 6
    @test count(managed(5, :threshold_cut).cut) == 0
    @test count(managed(6, :threshold_cut; month=36).cut) == 0
    # testing: around a symptomatic tree, the trees within the radius are tested (all of them
    # with test_share = 1); a visited or infected one is found with probability `sensitivity`
    function tested(sensitivity, seed)
        u = (state=fill(Int8(1), 20, 30), age=fill(Int16(60), 20, 30),
             clock=zeros(Int16, 20, 30), cut=fill(false, 20, 30),
             replanted=fill(false, 20, 30), tested=fill(false, 20, 30))
        u.state[10, 10] = 4
        u.state[9, 10] = 2
        u.state[11, 11] = 3
        u.state[10, 13] = 2
        q = merge(p,
                  (; strategy=:test_cut, interval=24, first_decision=48, test_share=1.0,
                   sensitivity, radius=1))
        return western_x_manage!(u, q, 48, Random.Xoshiro(seed))
    end
    v = tested(1.0, 1)
    @test count(v.tested) == 8 && v.state[9, 10] == 0 && v.state[11, 11] == 0
    @test v.state[10, 10] == 0 && v.state[10, 13] == 2 && count(v.cut) == 3
    found = sum(s -> count(tested(0.85, s).cut) - 1, 1:500)
    @test close_to(found / 1000, 0.85, 1000)

    # economics, from the source's yield y(a) and cost c(a) a month by age in months: the
    # disease-free orchard of 600 trees over 25 years
    y(a) = a < 24 ? 0 : a < 36 ? 0.9 : a < 48 ? 1.8 : a < 60 ? 2.2 : 2.7
    c(a) = a < 12 ? 2.0 : a < 24 ? 0.6 : a < 36 ? 1.9 : a < 48 ? 2.9 : a < 60 ? 3.4 : 3.9
    sol = solve(problem(:western_x_orchard; disease=false); saveat=1)
    e = western_x_economics(sol)
    δ = 1.05^(-1 / 12)
    @test sum(e.yield) ≈ 600 * sum(y, 1:300)
    @test e.yield[3] ≈ 600 * sum(y, 25:36)
    @test e.npv ≈ 600 * sum(k -> δ^k * (2.35y(k) - c(k)), 1:300)
    # losses scale the yield; at 25% symptomatic nothing is harvested; cuts, plants and tests cost
    u = orchard(fill(2, 20, 30))
    fake(u) = (t=[0, 1],
               u=[u,
                  merge(u,
                        (; cut=fill(false, 20, 30), replanted=fill(false, 20, 30),
                         tested=fill(false, 20, 30)))])
    @test western_x_economics(fake(u)).yield[1] ≈ 600 * 0.9 * 2.7
    S = fill(1, 20, 30)
    S[1:150] .= 5
    @test western_x_economics(fake(orchard(S))).yield[1] == 0
    S[150] = 1
    @test western_x_economics(fake(orchard(S))).yield[1] ≈ (451 + 149 * 0.25) * 2.7
    w = fake(orchard(fill(1, 20, 30)))
    w.u[2].cut[1:2] .= true
    w.u[2].replanted[1] = true
    w.u[2].tested[1:3] .= true
    @test western_x_economics(w).net[1] ≈ 600 * (2.35 * 2.7 - 3.9) - 2 * 33 - 14.3 - 3 * 300

    # the setup: 2-8% of a newly planted orchard visited; none without disease
    for s in 1:5
        u0 = problem(:western_x_orchard; seed=s).u0
        @test 12 <= count(==(2), u0.state) <= 48 && all(==(0), u0.age)
    end
    @test count(==(2), problem(:western_x_orchard; initial=0.1).u0.state) == 60
    @test count(==(2), problem(:western_x_orchard; disease=false).u0.state) == 0
    prob = problem(:western_x_orchard; strategy=:test_replant, years=6, seed=4)
    @test prob.p.first_decision == 48 && prob.p.interval == 24
    s1 = solve(prob; saveat=1)
    @test s1.u[end] == solve(prob; saveat=1).u[end] && length(s1.u) == 73
    @test_throws ArgumentError problem(:western_x_orchard; strategy=:spray)
    @test_throws ArgumentError problem(:western_x_orchard; initial=(0.5, 0.2))
    @test_throws ArgumentError problem(:western_x_orchard; pmap=(; rates=:weekly))
    @test_throws ArgumentError western_x_economics(solve(prob; saveat=12))
end

@testset "sudden_oak_death" begin
    # the kernel: a spore from a cell's centre travels in a uniform direction a half-Cauchy
    # distance α tan(πu/2); the landing probabilities match a direct simulation
    k = sod_kernel(20.57; radius=2)
    rng = Random.Xoshiro(7)
    n = 400_000
    hits = zeros(Int, 5, 5)
    for _ in 1:n
        d = 20.57 * tan(π * rand(rng) / 2)
        θ = 2π * rand(rng)
        i, j = round(Int, d * cos(θ) / 250), round(Int, d * sin(θ) / 250)
        abs(i) <= 2 && abs(j) <= 2 && (hits[i + 3, j + 3] += 1)
    end
    @test all(abs.(hits ./ n .- k) .<= 4 .* sqrt.(k .* (1 .- k) ./ n) .+ 1e-4)
    @test k ≈ reverse(k; dims=1) ≈ permutedims(k) && 0.89 < k[3, 3] < 0.93
    long = sod_kernel(9504; radius=3)
    @test sum(long) < 0.1 &&
          2 / π * atan(125 / 9504) < long[4, 4] < 2 / π * atan(125 * sqrt(2) / 9504)
    # the temperature index of Davidson et al. (2005), at least 0 (the cubic is negative at 35 °C)
    @test sod_temperature_index(15) ≈ -0.066 + 0.056 * 15
    @test sod_temperature_index(20) ≈ -0.066 + 1.12 - 0.0036 * 25 - 0.0003 * 125
    @test sod_temperature_index(35) == 0

    # one week: a cell of 100 infected, sporulating units at the centre of a 9 × 9 torus, and
    # cells of 100 susceptible units around it, with weather index 1: the spore units landing on
    # a cell at offset o are Poisson with mean β I (γ P_short(o) + (1 - γ) P_long(o)), and each
    # infects with probability S/nmax = 1
    short = sod_kernel(20.57; radius=2)
    model = build_sudden_oak_death(; short, long, method=:direct)
    p = merge(defaults(model), (; weather=ones(52)))
    function week(seed; I=Int8(100), forest=Int8(1), t0=0, q=p)
        u = (S=fill(Int8(100), 9, 9), I=zeros(Int8, 9, 9), forest=fill(forest, 9, 9),
             sporulating=fill(false, 9, 9))
        u.S[5, 5], u.I[5, 5], u.sporulating[5, 5] = 0, I, true
        prob = CAProblem(model, u, (t0, t0 + 1), q; space=Torus(9, 9), seed)
        return solve(prob).u[end]
    end
    Pshort(o) = maximum(abs, o) <= 2 ? short[3 + o[1], 3 + o[2]] : 0.0
    λ(o) = 4.4 * 100 * (0.9947 * Pshort(o) + 0.0053 * long[4 + o[1], 4 + o[2]])
    runs = 2000
    new = [week(s).I[6, 5] for s in 1:runs]
    @test abs(sum(new) / runs - λ((1, 0))) <= 4 * sqrt(λ((1, 0)) / runs)
    far = [week(s).I[8, 5] for s in 1:runs]
    @test abs(sum(far) / runs - λ((3, 0))) <= 4 * sqrt(λ((3, 0)) / runs) + 0.01
    # without long-range dispersal, nothing reaches 3 cells away (the short kernel's radius is 2)
    nolong = merge(p, (; long_factor=0.0, management_week=1))
    @test all(s -> week(s; q=nolong).I[8, 5] == 0, 1:200)
    # new infections never exceed the susceptible units
    @test all(<=(100), week(1; q=merge(p, (; beta=1e6))).I)
    @test all(==(100), week(1; q=merge(p, (; beta=1e6))).I[4:6, 4:6])
    # seasons: mixed evergreen forest takes no part in January (week 1), redwood-tanoak does;
    # neither does from October (week 40 is in October)
    @test sum(Int, week(1; forest=Int8(2)).I) == 100
    @test sum(Int, week(1; forest=Int8(1)).I) > 100
    @test sum(Int, week(1; forest=Int8(1), t0=39).I) == 100
    # a cell sporulates in the coming week if infected and in season then
    @test week(1; t0=37).sporulating[5, 5] && !week(1; t0=38).sporulating[5, 5]

    # the setup: half the cells of the landscape hold hosts, three sites with one infected unit,
    # the torus beyond the landscape masked; weather of 41 years, the future drawn by scenario
    prob = problem(:sudden_oak_death; n=40, seed=2)
    u0 = prob.u0
    hosts = count(>(0), u0.forest[1:40, 1:40])
    @test size(u0.S) == (80, 80) && all(==(0), u0.forest[41:80, :])
    @test 0.4 < hosts / 1600 < 0.6 && sum(Int, u0.I) == 3
    @test all(<=(100), u0.S .+ u0.I) && length(prob.p.weather) == 52 * 41
    W = [problem(:sudden_oak_death; n=40, seed=2, weather=w).p.weather
         for w in
             (:favourable, :unfavourable)]
    @test sum(W[1][(52 * 19 + 1):end]) > sum(W[2][(52 * 19 + 1):end])
    @test W[1][1:(52 * 19)] == W[2][1:(52 * 19)]
    m = problem(:sudden_oak_death; n=40, management=(2009, 0.1)).p
    @test m.management_week == 52 * 19 + 1 && m.long_factor == 0.1
    s1 = solve(problem(:sudden_oak_death; n=20, years=1990:1991, observed=1990:1991,
                       seed=3);
               saveat=104)
    @test s1.u[end] == solve(problem(:sudden_oak_death; n=20, years=1990:1991,
                                     observed=1990:1991, seed=3); saveat=104).u[end]
    @test sod_infected_area(s1.u[end]) == 0.0625 * count(>(0), s1.u[end].I[1:20, 1:20])
    @test_throws ArgumentError problem(:sudden_oak_death; weather=:wet)
    @test_throws ArgumentError problem(:sudden_oak_death; observed=1995:2000)
    @test_throws ArgumentError sod_kernel(-1.0; radius=3)
end

@testset "emerald_ash_borer" begin
    model = build_emerald_ash_borer()
    p = defaults(model)
    blank(N; phloem=200.0) = (phloem=fill(phloem, N, N), adults=zeros(N, N),
                              larvae1=zeros(N, N),
                              larvae2_odd=zeros(N, N), larvae2_even=zeros(N, N),
                              border=fill(false, N, N), phloem0=fill(phloem, N, N))
    year(u; t0=0) = solve(CAProblem(model, u, (t0, t0 + 1), p;
                                    space=Grid(size(u.phloem)...))).u[end]
    # growth (Eq. 3): r0 N larvae while the phloem allows A/C; they eat C each; a share
    # -0.194 ln(x) + 0.986 of them, x larvae per m² of phloem left, take two years
    u = blank(3)
    u.adults[2, 2] = 50.0
    v = year(u)
    new = 10.35 * 50
    left = 200 - 0.0113 * new
    p2 = -0.194 * log((new + 0.001) / (left + 0.0001)) + 0.986
    @test v.phloem[2, 2] ≈ left && v.larvae1[2, 2] ≈ (1 - p2) * new
    @test v.larvae2_even[2, 2] ≈ p2 * new && v.larvae2_odd[2, 2] == 0 && v.adults[2, 2] == 0
    # in the next (odd) year, the two-year larvae go to the other slot and the earlier ones stay
    w = year(merge(v, (; adults=copy(u.adults))); t0=1)
    @test w.larvae2_odd[2, 2] > 0 && w.larvae2_even[2, 2] ≈ p2 * new
    # too many adults for the phloem: A/C larvae, no phloem left, and all develop in a year
    u.adults[2, 2] = 1e6
    v = year(u)
    @test v.phloem[2, 2] == 0 && v.larvae1[2, 2] ≈ 200 / 0.0113 && v.larvae2_even[2, 2] == 0

    # dispersal (Eq. 1): the adults of a cell spread in proportion to exp(-b D), D in metres
    # between cell centres (12.5 m for the cell itself), normalised over the 101 × 101 window
    N = 121
    u = blank(N)
    c = 61
    u.larvae1[c, c] = 1000.0
    eab_disperse!(u, merge(p, (; behaviour=:standard)), 0)
    K = [exp(-0.037 * 50 * (i == j == 0 ? 0.25 : hypot(i, j))) for i in -50:50, j in -50:50]
    @test u.adults[c + 1, c] ≈ 1000 * K[52, 51] / sum(K)
    @test u.adults[c, c] ≈ 1000 * K[51, 51] / sum(K) && sum(u.adults) ≈ 1000
    @test u.adults[c + 51, c] == 0
    # avoidance: a cell without phloem is `exit` cells further away
    u = blank(N)
    u.larvae1[c, c] = 1000.0
    u.phloem[c + 1, c] = 0.0
    eab_disperse!(u, p, 0)
    Kx = copy(K)
    Kx[52, 51] *= exp(-0.037 * 50 * 36)
    @test u.adults[c + 1, c] ≈ 1000 * Kx[52, 51] / sum(Kx)
    # a round-off of phloem is phloem unless it is at most `bare`
    u.phloem[c + 1, c] = 1e-12
    eab_disperse!(u, p, 0)
    @test u.adults[c + 1, c] ≈ 1000 * K[52, 51] / sum(K)
    eab_disperse!(u, merge(p, (; bare=1e-9)), 0)
    @test u.adults[c + 1, c] ≈ 1000 * Kx[52, 51] / sum(Kx)
    # transposed sources (the source's code): the beetles of cell (i, j) spread from (j, i)
    u = blank(N)
    u.larvae1[c - 10, c + 5] = 1000.0
    eab_disperse!(u, merge(p, (; behaviour=:standard, transposed=true)), 0)
    @test u.adults[c + 5, c - 10] ≈ 1000 * K[51, 51] / sum(K) && sum(u.adults) ≈ 1000
    @test u.adults[c + 6, c - 10] ≈ 1000 * K[52, 51] / sum(K)
    v = (; (k => (x isa Matrix ? x[:, 1:(N - 1)] : x) for (k, x) in pairs(blank(N)))...)
    @test_throws ArgumentError eab_disperse!(v, merge(p, (; transposed=true)), 0)
    # beetles reaching the border are lost; the parity of the year decides which two-year larvae
    # emerge
    u = blank(N)
    u.border[1:2, :] .= true
    u.larvae2_odd[3, c] = 100.0
    eab_disperse!(u, p, 1)
    @test 0 < sum(u.adults) < 100 && all(==(0), u.adults[1:2, :])
    u.adults .= 0
    eab_disperse!(u, p, 0)
    @test sum(u.adults) == 0

    # environments of 8,405,000 m² on 205 × 205 cells
    for kind in (:homogeneous, :clumped, :random)
        @test sum(eab_environment(kind)) ≈ 8_405_000
    end
    clumped = eab_environment(:clumped)
    @test clumped[103, 103] > 0 && all(==(0), clumped[103, 106:110]) &&
          clumped[103, 113] > 0
    r = eab_environment(:random; rng=Random.Xoshiro(1))
    @test r[103, 103] > 0 && length(unique(round.(r; digits=6))) == 17
    # the setup: founders at the centre, dispersed for the first year, in a border of 2 cells
    prob = problem(:emerald_ash_borer; n=101, years=1)
    u0 = prob.u0
    @test size(u0.phloem) == (105, 105) && count(u0.border) == 105^2 - 101^2
    @test u0.larvae1[53, 53] == 100 && sum(u0.adults) ≈ 100 && all(==(200), u0.phloem)
    # the first year: 100 founders give 1035 larvae, within 150 m (the source's Fig. 9)
    s = eab_spread(solve(problem(:emerald_ash_borer; years=1); saveat=1).u[end])
    @test s.beetles ≈ 1035 && s.radius == 150
    @test eab_spread(u0).radius == 0
    @test_throws ArgumentError problem(:emerald_ash_borer; behaviour=:stress)
    @test_throws ArgumentError eab_environment(:clumped; n=100)
    @test_throws ArgumentError problem(:emerald_ash_borer; environment=:sparse)
end

@testset "eab_management" begin
    model = build_eab_management()
    p = defaults(model)
    C = 0.0113
    blank(N; phloem=200.0, large=0.0) = (phloem=fill(phloem, N, N), adults=zeros(N, N),
                                         larvae1=zeros(N, N), larvae2_odd=zeros(N, N),
                                         larvae2_even=zeros(N, N), border=fill(false, N, N),
                                         phloem0=fill(phloem, N, N),
                                         treated=fill(true, N, N),
                                         large=fill(large, N, N))
    year(u, q; t=0) = solve(CAProblem(model, u, (t, t + 1), q;
                                      space=Grid(size(u.phloem)...))).u[end]
    total(v) = v.larvae1[2, 2] + v.larvae2_odd[2, 2] + v.larvae2_even[2, 2]
    # girdling (the source's code): a tree of g m² takes y g / (y g + A - g) of the larvae, with
    # y = max(18.25 (x + 1)^-0.533, 1) at x = N / A, and is removed with them
    g = merge(p, (; treatment=:girdling))
    u = blank(3; large=100.0)
    u.adults[2, 2] = 50.0
    v = year(u, g)
    N = 10.35 * 50
    x = N / 200.0001
    y = 18.25 * (x + 1)^-0.533
    kept = N - N * 10.32 * y / (10.32 * y + 200 - 10.32)
    @test total(v) ≈ kept && v.phloem[2, 2] ≈ 200 - C * kept - 10.32
    @test v.large[2, 2] ≈ 100 - 10.32 - (100 - 10.32 + 0.0001) / 200.01 * C * kept
    @test total(year(u, merge(g, (; attraction_decline=false)))) ≈
          N - N * 10.32 * 18.25 / (10.32 * 18.25 + 200 - 10.32)
    # the ratio never falls below 1, nor does a cell girdle more than its large trees
    u = blank(3; large=4.0)
    u.adults[2, 2] = 50.0
    v = year(u, merge(g, (; girdle_a=0.5)))
    @test total(v) ≈ N - N * 4 / 200 && v.large[2, 2] == 0
    # before `start` nothing is treated
    u = blank(3; large=100.0)
    u.adults[2, 2] = 50.0
    @test total(year(u, merge(g, (; start=1)))) ≈ N
    # insecticide: larvae laid in the large trees die (x times their area), and feeding is
    # reduced by the unprotected share, as in the source's code; the large trees are not eaten
    u = blank(3; large=50.0)
    u.adults[2, 2] = 1e6
    v = year(u, merge(p, (; treatment=:insecticide)))
    N0 = 200 / C
    N = N0 - N0 / 200.0001 * 50
    @test total(v) ≈ N && v.phloem[2, 2] ≈ 200 - (150.0001 / 200.0001) * C * N
    @test v.large[2, 2] == 50
    # ash reduction, in the year `start` only: the large trees left after feeding are removed,
    # and the source's code then sets every cell's large-tree phloem to its phloem
    u = blank(3; large=100.0)
    u.adults[2, 2] = 1.0
    r = merge(p, (; treatment=:ash_reduction, start=2))
    v = year(u, r; t=2)
    eaten = C * 10.35
    left = 100 - (100.0001 / 200.01) * eaten
    @test v.phloem[2, 2] ≈ 200 - eaten - left && v.large[2, 2] == v.phloem[2, 2]
    @test year(u, r; t=3).phloem[2, 2] ≈ 200 - eaten

    # the stress bias (the source's code): 7 in every treated cell under girdling, else
    # 7 lost / (max(A0, 3) / 3), kept within 1 and 7
    u = blank(3)
    u.phloem[1, 1] = 190.0
    u.phloem[2, 2] = 200 * (1 - 0.1)
    u.phloem[3, 3] = 100.0
    A = eab_attraction(u, p)
    @test A[1, 1] ≈ 1.05 && A[2, 2] ≈ 2.1 && A[3, 3] == 7 && A[1, 2] == 1
    @test all(==(7), eab_attraction(u, g))
    # dispersal weights each receiving cell by its attraction, normalised source by source over
    # the 81 x 81 window; cells with no more than 0.9 adults disperse none
    M = 101
    c = 51
    u = blank(M)
    u.treated .= false
    u.larvae1[c, c] = 1000.0
    u.larvae1[5, 5] = 0.9
    u.phloem[c + 1, c] = 100.0
    eab_manage!(u, merge(p, (; stress_bias=true)), 1)
    K = [exp(-0.037 * 50 * (i == j == 0 ? 0.25 : hypot(i, j))) for i in -40:40, j in -40:40]
    W = copy(K)
    W[42, 41] *= 7     # half its phloem lost, more than a third
    @test u.adults[c + 1, c] ≈ 1000 * W[42, 41] / sum(W) && sum(u.adults) ≈ 1000
    @test u.adults[c + 41, c] == 0
    # under ash reduction, adults landing in a treated cell lose the share on its large trees
    u = blank(M; large=100.0)
    u.larvae1[c, c] = 1000.0
    eab_manage!(u, merge(p, (; treatment=:ash_reduction)), 0)
    @test sum(u.adults) ≈ 1000 * (1 - 100 / 200.0001)

    # with no treatment and no stress bias this is the base model, exactly: the two extra layers
    # do not feed back, and the management callback disperses as the base model's does
    plain = problem(:eab_management; n=30, years=8)
    base = build_emerald_ash_borer()
    names = Tuple(s.name for s in base.states)
    pbase = (; (k => v for (k, v) in pairs(plain.p) if haskey(defaults(base), k))...)
    baseprob = CAProblem(base, NamedTuple{names}(plain.u0), 8, pbase; space=plain.space,
                         seed=plain.seed, callback=eab_callback())
    @test map(u -> NamedTuple{names}(u), solve(plain).u) == solve(baseprob).u
    treated = problem(:eab_management; n=30, years=8, treatment=:girdling)
    @test map(u -> NamedTuple{names}(u), solve(treated).u) != solve(baseprob).u
    # the setup: the origin is the cell (n ÷ 2 + 1, n ÷ 2 + 1); without treatment the 100
    # founders give 1035 larvae within 150 m after a year; large trees under 6.32 m² do not count
    s = eab_spread(solve(problem(:eab_management; years=1)).u[end])
    @test s.beetles ≈ 1035 && s.radius == 150
    prob = problem(:eab_management; treatment=:ash_reduction, years=1)
    @test count(prob.u0.treated) == count(i^2 + j^2 <= 36 for i in -6:6, j in -6:6)
    @test prob.u0.treated[51, 51] && prob.u0.treated[57, 51] && !prob.u0.treated[58, 51]
    @test all(==(100), problem(:eab_management; years=1).u0.large)
    @test all(==(0), problem(:eab_management; share=0.03, years=1).u0.large)
    @test eab_beyond(problem(:eab_management; years=1).u0, -1) == 100
    @test_throws ArgumentError problem(:eab_management; treatment=:burning)
    @test_throws ArgumentError problem(:eab_management; share=1.5)
    @test_throws ArgumentError problem(:eab_management; phloem=ones(3, 3))
    @test_throws ArgumentError problem(:eab_management; pmap=(; behaviour=:resource_bias))
end

@testset "ipydisp_dispersal" begin
    # haversine on a sphere of 6371 km: one degree of latitude is 6371 pi / 180 km
    @test ipydisp_distance(0, 0, 1, 0) ≈ 6371 * pi / 180
    @test ipydisp_distance(60, 10, 60, 11) ≈ 2 * 6371 * asin(cosd(60) * sind(0.5))
    # the grid: cells of cellsize x 0.00833 degrees from the north-west corner, row by row,
    # the last row and column cut at the box
    g = ipydisp_grid(; west=0, east=0.02, south=0, north=0.01, cellsize=1)
    @test length(g.lat) == 6 && g.row == [1, 1, 1, 2, 2, 2] && g.col == [1, 2, 3, 1, 2, 3]
    @test g.lon[1] ≈ 0.00833 / 2 && g.lat[1] ≈ 0.01 - 0.00833 / 2
    @test g.lon[3] ≈ (2 * 0.00833 + 0.02) / 2 && g.lat[4] ≈ (0.01 - 0.00833) / 2
    @test length(ipydisp_grid(; west=0, east=0.02, south=0, north=0.01,
                              inside=(la, lo) -> lo < 0.01).lat) == 2
    # neighbours: other cells closer than the travel distance
    nb = ipydisp_neighbours([0.0, 0.0, 0.0], [0.0, 0.1, 0.25], 12)
    @test nb == [[2], [1], Int[]]
    # constraints: lower < value <= upper, every one met
    @test ipydisp_suitable([([0.1, 0.3, 0.5, 0.9], 0.3, 0.9), ([1, 1, 1, 0], 0, 1)], 4) ==
          [false, false, true, false]
    @test_throws ArgumentError ipydisp_suitable([([0.1], 0, 1)], 2)

    # one step on a path 1 - 2 - 3 - 4 - 5 with cell 3 unsuitable: cells next to an exposed or
    # colonised cell become colonised if suitable, else exposed, and exposed cells spread on
    space = GraphSpace(5, [(1, 2), (2, 3), (3, 4), (4, 5)])
    u0 = (; status=Int8[0, 2, 0, 0, 0], suitable=[true, true, false, true, true])
    sol = solve(CAProblem(build_ipydisp_dispersal(), u0, 3, (;); space); saveat=1)
    @test [u.status for u in sol.u] ==
          [Int8[0, 2, 0, 0, 0], Int8[2, 2, 1, 0, 0], Int8[2, 2, 1, 2, 0],
           Int8[2, 2, 1, 2, 2]]
    @test ipydisp_counts(sol.u[3]) == (; unexposed=1, exposed=1, colonised=3)

    # the example: about 9.3 km cells joined to their eight neighbours by 15 km of travel
    ex = ipydisp_example()
    @test maximum(length, ex.neighbours) == 8
    @test ipydisp_suitable(ex.constraints, length(ex.grid.lat)) ==
          (ex.constraints[1][1] .> 0.3)
    prob = problem(:ipydisp_dispersal; steps=60)
    @test count(==(2), prob.u0.status) == 1
    final = solve(prob).u[end]
    @test ipydisp_counts(final).unexposed == 0
    @test all(final.status[final.suitable] .== 2) &&
          all(final.status[.!final.suitable] .== 1)
    @test_throws ArgumentError problem(:ipydisp_dispersal; grid=ex.grid)
    @test_throws ArgumentError problem(:ipydisp_dispersal; start=[10^6])
    @test_throws ArgumentError ipydisp_grid(; west=1, east=0, south=0, north=1)
end

@testset "dune_formation" begin
    model = build_dune_formation()
    p = defaults(model)
    quiet = merge(p, (; hydro_every=10^6, iterations=10^6))
    bed(m, n, z; h=0.02) = (topo=fill(Int32(round(z / h)), m, n), pioneer=zeros(m, n),
                            stabilizer=zeros(m, n), balance=zeros(Int32, m, n),
                            zref=zeros(Int32, m, n), sea=zeros(m, n))
    # transport capacity (Eq. 3): Q = h L (pe0/pd0) n = 25 m³/m/yr on a dry bare bed, for slabs
    # of 0.1 m hopping 1 cell and of 0.02 m hopping 5; counted crossing the step at column 200
    for (h, L) in ((0.1, 1), (0.02, 5))
        q = merge(quiet, (; h, L, rho_star=Inf))
        total = 0
        for t in 1:20
            u = bed(40, 400, 2.9; h)
            u.topo[:, 200:end] .= round(Int32, 3.1 / h)
            total += dune_step!(u, q, t, Random.Xoshiro(t))
        end
        @test abs(total * h * 50 / (40 * 20) - h * L * 0.5 / 0.1 * 50) < 2.5
    end
    # a slab hopping 5 cells lands part-way only where vegetation is at least rho_star, with
    # probability ρ (1 - pd0) + pd0 (Eq. 5); here ρ = 0.5 two cells downwind of its source
    function trap(ρ; runs=2000)
        landed = 0
        for s in 1:runs
            u = bed(1, 20, -1.0)
            u.topo[1, 5] = 50
            u.pioneer[1, 7] = ρ
            dune_step!(u, merge(quiet, (; pe0=1.0)), 1, Random.Xoshiro(s))
            landed += u.topo[1, 7] > -50
        end
        return landed / runs
    end
    @test close_to(trap(0.5), 0.5 * 0.9 + 0.1, 2000)
    @test trap(0.2) == 0
    # shadow: cells in the 15° lee of a step cannot be eroded; avalanches leave no cell whose
    # steepest drop exceeds the angle of repose (30° bare: 28 slabs of 0.02 m orthogonally, 40
    # along a diagonal), as in DuBeVeg's code
    u = bed(1, 60, 1.0)
    u.topo[1, 10] = round(Int32, 3.0 / 0.02)
    shade = falses(1, 60)
    EcologicalCellularAutomata._dune_shadow!(shade, u.topo, p)
    lee = (2.0 - 0.02) / tand(15)
    @test all(shade[1, 11:(10 + floor(Int, lee))]) && !shade[1, 10 + ceil(Int, lee) + 1]
    u = bed(9, 9, 1.0)
    u.topo[5, 5] = round(Int32, 5.0 / 0.02)
    dune_step!(u, merge(quiet, (; pe0=0.0)), 1, Random.Xoshiro(1))
    limit, diagonal = floor(Int, tand(30) / 0.02), floor(Int, sqrt(2) * tand(30) / 0.02)
    @test (limit, diagonal) == (28, 40)
    @test sum(u.topo) == 80 * 50 + 250
    steepest_within(u, I) = all(J -> !checkbounds(Bool, u.topo, J) ||
                                     u.topo[I] - u.topo[J] <
                                     maximum(u.topo[I] - u.topo[K]
                                             for K in I .+ CartesianIndex.(-1:1, (-1:1)')
                                             if checkbounds(Bool, u.topo, K)) ||
                                     u.topo[I] - u.topo[J] <=
                                     (J[1] != I[1] && J[2] != I[2] ? diagonal : limit),
                                I .+ CartesianIndex.(-1:1, (-1:1)'))
    @test all(I -> steepest_within(u, I), CartesianIndices(u.topo))
    # only the steepest drop counts: a drop of 30 orthogonally (over 28) stays when the steepest,
    # 35 along a diagonal, is within 40; with it at 41 one slab goes down the diagonal
    u = bed(2, 2, 0.0)
    u.topo .= Int32[100 70; 70 65]
    EcologicalCellularAutomata._dune_avalanche!(u, p, Random.Xoshiro(1))
    @test u.topo == Int32[100 70; 70 65]
    u.topo .= Int32[100 70; 70 59]
    EcologicalCellularAutomata._dune_avalanche!(u, p, Random.Xoshiro(1))
    @test u.topo == Int32[99 70; 70 60] && u.balance == Int32[-1 0; 0 1]
    # the sea (Eq. 4): below the water level, each slab above the reference surface goes with
    # probability (1 - Rv ρ)(wave + Pinun); with no dissipation and a wave strength of 1 (0.64 m
    # of water, 32 slabs, times wave_energy 1/32), certainly when bare, half of them under full
    # plant cover with Rv = 0.5; the plants themselves stay until the year's end
    sea = merge(p,
                (; Wdiss=0.0, water_levels=[0.64], wave_energy=1 / 32, pe0=0.0,
                 iterations=10^6))
    u = bed(200, 3, 0.2)
    u.pioneer[:, 2] .= 1.0
    dune_step!(u, sea, 2, Random.Xoshiro(1))
    @test all(==(0), u.topo[:, 1])
    @test close_to(sum(10 .- u.topo[:, 2]) / 2000, 0.5, 2000)
    # slab by slab, not the whole column at once: a binomial(10, 1/2) is 0 or 10 once in 512
    @test count(c -> c == 0 || c == 10, u.topo[:, 2]) <= 5
    @test all(==(1), u.pioneer[:, 2])
    @test all(==(1), u.sea[:, 1]) && all(==(0.5), u.sea[:, 2])
    # slabs missing below the reference surface appear with the same probability
    u = bed(200, 2, 0.2)
    u.zref .= 20
    u.pioneer .= 1.0
    dune_step!(u, sea, 2, Random.Xoshiro(2))
    @test close_to(sum(u.topo .- 10) / 4000, 0.5, 4000) && all(<=(20), u.topo)
    # the year's probabilities add, to no more than 1
    u = bed(1, 2, 0.2)
    u.pioneer .= 1.0
    for (k, t) in enumerate((2, 2, 2))
        dune_step!(u, sea, t, Random.Xoshiro(k))
        @test u.sea[1, 1] == min(1.0, 0.5k)
    end
    # DuBeVeg's wave strength: (WE - Wdiss D) times wave_energy times the water level in slabs,
    # D the sum of 1/depth (m) over the cells crossed, no less than wave_min; under full cover
    # here, 0.44 m of water over a flat bed
    sea = merge(p, (; water_levels=[0.64], pe0=0.0, iterations=10^6))
    u = bed(1, 60, 0.2)
    u.pioneer .= 1.0
    dune_step!(u, sea, 2, Random.Xoshiro(1))
    @test u.sea[1, :] ≈ [min(1, 0.5 * max((1 - 0.012 * (x - 1) / 0.44) * 32, 0.1))
                         for x in 1:60]
    @test u.sea[1, 1] == 1 && u.sea[1, 60] ≈ 0.05
    # at the year's end plant cover is multiplied by one minus the sea's probabilities that year
    u = bed(1, 2, 2.0)
    u.stabilizer .= 0.5
    u.sea[1, 1] = 0.4
    dune_step!(u, merge(p, (; pe0=0.0, hydro_every=10^6)), 50, Random.Xoshiro(1))
    @test u.stabilizer[1, 1] ≈ (0.5 + 0.05) * 0.6 && u.stabilizer[1, 2] ≈ 0.55 &&
          all(==(0), u.sea)
    # plants: a year with 0.6 m of burial adds 0.2 to pioneer cover and takes some from the
    # stabilizer (beyond its range); 0.1 m of erosion kills neither
    u = bed(1, 3, 2.0)
    u.pioneer .= 0.3
    u.stabilizer .= 0.3
    u.balance[1, 1] = round(Int32, 0.6 / 0.02)
    u.balance[1, 2] = round(Int32, -0.1 / 0.02)
    dune_step!(u, merge(p, (; pe0=0.0, hydro_every=10^6)), 50, Random.Xoshiro(1))
    @test u.pioneer[1, 1] ≈ 0.5 && u.stabilizer[1, 1] ≈ 0.3 - (0.6 - 0.2) / (2.8 - 0.2)
    @test u.pioneer[1, 2] ≈ 0.3 - (1 - (-0.1 + 1.4) / 1.6) &&
          u.stabilizer[1, 2] ≈ 0.3 + 0.05 * 0.55 / 0.65
    @test all(==(0), u.balance)
    # a gain is at most the cover left free at the year's start; a loss is not capped
    u = bed(1, 1, 2.0)
    u.pioneer .= 0.5
    u.stabilizer .= 0.45
    u.balance[1, 1] = round(Int32, 0.6 / 0.02)
    dune_step!(u, merge(p, (; pe0=0.0, hydro_every=10^6)), 50, Random.Xoshiro(1))
    @test u.pioneer[1, 1] ≈ 0.55 && u.stabilizer[1, 1] ≈ 0.45 - (0.6 - 0.2) / (2.8 - 0.2)
    # each species spreads from its own: certain lateral spread and no new establishment, a
    # bare cell beside pioneers gets pioneers (0.6 m of burial, 0.2), and one beside stabilizers
    # gets stabilizers under 0.1 m of erosion, where the balance is not zero
    grow = merge(p, (; pe0=0.0, hydro_every=10^6, p_lateral=1.0, p_establish=0.0))
    u = bed(1, 2, 2.0)
    u.pioneer[1, 1] = 0.5
    u.balance[1, 2] = round(Int32, 0.6 / 0.02)
    dune_step!(u, grow, 50, Random.Xoshiro(1))
    @test u.pioneer[1, 2] ≈ 0.2 && u.stabilizer[1, 2] == 0
    u = bed(1, 2, 2.0)
    u.stabilizer[1, 1] = 0.5
    u.balance[1, 2] = round(Int32, -0.1 / 0.02)
    dune_step!(u, grow, 50, Random.Xoshiro(1))
    @test u.stabilizer[1, 2] ≈ 0.05 * 0.55 / 0.65 && u.pioneer[1, 2] == 0
    # ...and not from the other species: beside only stabilizers, burial that suits pioneers
    # brings none
    u = bed(1, 2, 2.0)
    u.stabilizer[1, 1] = 0.5
    u.balance[1, 2] = round(Int32, 0.6 / 0.02)
    dune_step!(u, grow, 50, Random.Xoshiro(1))
    @test u.pioneer[1, 2] == 0 && u.stabilizer[1, 2] == 0
    # lateral spread with its yearly probability
    u = bed(1000, 2, 2.0)
    u.pioneer[:, 1] .= 0.5
    u.balance[:, 2] .= round(Int32, 0.6 / 0.02)
    dune_step!(u, merge(grow, (; p_lateral=0.2)), 50, Random.Xoshiro(1))
    @test close_to(count(>(0), u.pioneer[:, 2]) / 1000, 0.2, 1000)
    # establishing anew: pioneers anywhere burial suits them from 0 to 20 m; stabilizers only
    # where the balance is zero
    settle = merge(p, (; pe0=0.0, hydro_every=10^6, p_lateral=0.0, p_establish=0.1))
    u = bed(2000, 1, 2.0)
    u.balance .= round(Int32, 0.6 / 0.02)
    dune_step!(u, settle, 50, Random.Xoshiro(1))
    @test close_to(count(>(0), u.pioneer) / 2000, 0.1, 2000) && all(==(0), u.stabilizer)
    u = bed(2000, 1, 2.0)
    dune_step!(u, settle, 50, Random.Xoshiro(2))
    @test close_to(count(>(0), u.stabilizer) / 2000, 0.1, 2000)
    u = bed(50, 1, 21.0)
    dune_step!(u, settle, 50, Random.Xoshiro(3))
    @test all(==(0), u.stabilizer) && all(==(0), u.pioneer)
    # wave run-up (Stockdon et al.) and the source's profile
    @test dune_runup(0.8, 0.01) == 0
    H, Lw = 2.931 * 2 - 2.637, 46.74 * 2 - 30.59
    @test dune_runup(2.0, 0.01) ≈ 0.043 * sqrt(H * Lw)
    prof = dune_profile(; beach=100, h=0.02)
    z = prof.elevation .* 0.02
    @test length(z) == 60 + 100 + 120 && z[60] < 0 && z[160] ≈ 2.5
    @test 14.5 < maximum(z) <= 15 && z[end] == 5 && prof.reference[200] == 125
    prob = problem(:dune_formation; scenario=:narrow, alongshore=4, years=1, seed=2)
    u0 = prob.u0
    cover = u0.pioneer .+ u0.stabilizer
    @test all(cover[u0.topo .> 300] .== 1) && all(cover[u0.topo .< 200] .== 0)
    @test length(prob.p.water_levels) == 25 && prob.p.h == 0.02 && prob.p.L == 5
    @test problem(:dune_formation; representation=:default, years=1, alongshore=2).p.L == 1

    # the water levels measured at Scheveningen, against the statistics of the source's Fig. 5d
    # read from the figure: median about 1.5 m, 25-75% 1.36-1.72 m, 90% 1.97 m, highest about 3.2 m
    observed = scheveningen_water_levels()
    ranked = sort(observed)
    level_at(share) = ranked[round(Int, share * length(ranked))]
    average(x) = sum(x) / length(x)
    @test length(observed) == 25 * 50
    @test observed[1:2] == [1.39, 1.13] && observed[end] == 2.32
    @test abs(level_at(0.5) - 1.5) < 0.05 && abs(level_at(0.25) - 1.36) < 0.05
    @test abs(level_at(0.75) - 1.72) < 0.05 &&
          abs(level_at(0.9) - 1.97) < 0.1
    @test abs(maximum(observed) - 3.2) < 0.1
    # higher in winter: the first and last cycles of the year against the middle ones
    season(k) = (k - 1) % 25 + 1
    winter = [l for (k, l) in enumerate(observed) if season(k) <= 3 || season(k) >= 23]
    summer = [l for (k, l) in enumerate(observed) if 9 <= season(k) <= 15]
    @test average(winter) > average(summer) + 0.3
    # the record is the default, repeats past its 50 years, and the series is the seed's to draw only
    # when asked; drawn from the problem's generator, so a seed still gives the same one
    @test scheveningen_water_levels(; years=2) == observed[1:50]
    @test scheveningen_water_levels(; years=51)[1251:1275] == observed[1:25]
    @test_throws ArgumentError scheveningen_water_levels(; years=-1)
    default_levels(; kw...) = problem(:dune_formation; alongshore=2, years=2, kw...).p.water_levels
    @test default_levels() == observed[1:50] == default_levels(; seed=9)
    @test default_levels(; levels=:observed) == observed[1:50]
    @test default_levels(; levels=:synthetic, seed=1) ==
          dune_water_levels(Random.Xoshiro(hash((1, 0, :dune))); years=2)
    @test default_levels(; levels=:synthetic, seed=1) !=
          default_levels(; levels=:synthetic, seed=2)
    @test_throws ArgumentError problem(:dune_formation; levels=:tidal)
    @test problem(:dune_formation; alongshore=2, years=1, pmap=(; water_levels=[1.0])).p.water_levels ==
          [1.0]
    m = dune_measures(u0)
    @test 161 < m.toe < 170 && m.width == 281 - m.toe
    @test_throws ArgumentError problem(:dune_formation; scenario=:steep)
    @test_throws ArgumentError problem(:dune_formation; representation=:fast)

    # the prograding beach (Teixeira et al. 2023, Fig. 2b): the narrow beach on a domain 200 m
    # longer out to sea, whose reference surface builds out 4 m a year with a berm at 2 m
    prob = problem(:dune_formation; scenario=:prograding, alongshore=3, years=50, seed=1)
    refs = prob.p.references
    n = size(prob.u0.topo, 2)
    @test n == 260 + 100 + 120 && length(refs) == 51
    @test refs[1] == prob.u0.zref[1, :]
    shore(z) = findfirst(>=(0), z)
    @test [shore(refs[y + 1]) for y in 0:10:50] == [261 - 4y for y in 0:10:50]
    # a 200 m berm at 2 m (100 slabs), next to the lower beach's top cell, also at 2 m
    @test count(==(100), refs[51]) == 201 && all(refs[51][(end - 119):end] .== 125)
    @test refs[51][(end - 139):end] == refs[1][(end - 139):end]      # the upper beach stays
    # the year's end replaces the reference surface in every row
    small = problem(:dune_formation; scenario=:prograding, alongshore=2, years=2, seed=1)
    u = deepcopy(small.u0)
    for t in 1:50
        dune_step!(u, small.p, t, Random.Xoshiro(t))
        t < 50 && @test u.zref[1, :] == small.p.references[1]
    end
    @test u.zref[1, :] == u.zref[2, :] == small.p.references[2]
    @test_throws ArgumentError problem(:dune_formation; scenario=:prograding, rate=-1)
end

@testset "landform profiles" begin
    # pieces join when the elevations meet; joining is associative with empty pieces as
    # identities, so profiles compose like paths between elevations
    a = shoreface_segment(30)
    b = slope_segment(40, 0, 2)
    c = flat_segment(10, 2)
    d = slope_segment(5, 2, 2.5) * ridge_segment(20, 2.5, 10, 4)
    @test ((a * b) * c * d).heights == (a * (b * (c * d))).heights
    @test [p.kind for p in (a * b * c * d).pieces] ==
          [:shoreface, :beach, :berm, :beach, :dune]
    @test (flat_segment(0, 0.0) * b) === b && (b * flat_segment(0, 2.0)) === b
    @test (a * b).z0 == a.z0 && (a * b).z1 == 2 && length(a * b) == 70
    @test_throws ArgumentError b * a
    @test_throws ArgumentError slope_segment(10, 0, 1) * flat_segment(5, 2)
    # level pieces at one elevation add, as the berm of a prograding beach does
    @test (flat_segment(3, 2) * flat_segment(4, 2)).heights == flat_segment(7, 2).heights
    @test dune_landform(; beach=100, berm=12).heights ==
          (shoreface_segment(60) * slope_segment(80, 0, 2) * flat_segment(5, 2) *
           flat_segment(7, 2) * slope_segment(20, 2, 2.5) * ridge_segment(100, 2.5, 15, 5) *
           flat_segment(20, 5; kind=:hinterland)).heights
    # the pieces' formulas: a Dean shoreface, straight slopes, a ridge peaking at its middle
    @test a.heights[end] ≈ -0.12 * 0.5^(2 / 3) && a.heights[1] ≈ -0.12 * 29.5^(2 / 3)
    @test b.heights == [2k / 40 for k in 1:40]
    r = ridge_segment(20, 2.5, 10, 4)
    @test r.heights[10] ≈ 10 && r.heights[20] ≈ 4 && r.heights[5] ≈ 2.5 + 7.5 * sin(π / 4)^2
    # the reference surface levels the dune and the land behind it at the dune's foot, piece by
    # piece, so it respects joins
    land = d * flat_segment(6, 4; kind=:hinterland)
    @test reference_profile(a * b * c * land).heights ==
          (a * b * c * reference_profile(land)).heights
    @test reference_profile(land).heights == vcat(d.pieces[1].heights, fill(2.5, 26))
    @test reference_profile(a * b) === a * b ||
          reference_profile(a * b).heights == (a * b).heights
    @test profile_heights(a * b, 10) == (a * b).heights[61:70]
    @test profile_heights(b, 40; h=0.5) == round.(Int32, b.heights ./ 0.5)
    @test_throws ArgumentError profile_heights(b, 41)
    @test_throws ArgumentError ridge_segment(1, 0, 1, 0)
    @test sprint(show, a * b) ==
          "ProfileSegment(shoreface * beach: $(a.z0) → 2.0 m over 70 m)"
    @test profile_kinds(c * land, 8) == [:dune, :dune, fill(:hinterland, 6)...]
    @test_throws ArgumentError profile_kinds(b, 41)
end

@testset "dune landforms" begin
    # other landforms from the same pieces, on the same domain as the source's
    fore = dune_landform(; beach=500)
    double = dune_landform(; form=:double_ridge, beach=500)
    wall = dune_landform(; form=:seawall, beach=100)
    breach = dune_landform(; form=:breach, beach=100)
    @test length(double) == length(fore) == 680 && length(wall) == length(breach) == 280
    @test [p.kind for p in double.pieces] ==
          [:shoreface, :beach, :dune, :slack, :dune, :hinterland]
    @test [p.kind for p in wall.pieces] == [:shoreface, :beach, :seawall, :promenade]
    # the double ridge keeps the foredune where it was, behind an embryo dune and a slack
    @test double.heights[561:end] == fore.heights[561:end]
    @test maximum(double.heights[361:400]) ≈ 7 && all(double.heights[401:560] .== 2.5)
    @test maximum(breach.heights[161:260]) ≈ 5 &&
          breach.heights[1:160] == wall.heights[1:160]
    @test wall.heights[161:end] == [5.0; fill(7.5, 119)]
    # the sea's reference surface levels the embryo dune too, at the dune foot
    @test all(reference_profile(double).heights[361:end] .== 2.5)
    @test reference_profile(wall).heights == wall.heights
    @test_throws ArgumentError dune_landform(; form=:cliff)
    @test_throws ArgumentError dune_landform(; form=:double_ridge, beach=100)
    @test_throws ArgumentError dune_landform(; form=:seawall, wall=2)
    @test_throws ArgumentError dune_landform(; form=:breach, berm=10)

    # the seawall and promenade are non-erodible, bare and their own reference surface
    prob = problem(:dune_formation; scenario=:seawall, alongshore=3, years=2, seed=1)
    u0 = prob.u0
    @test all(u0.hard[:, 161:end] .== u0.topo[:, 161:end])
    @test all(u0.hard[:, 1:160] .== typemin(Int32))
    @test all(u0.pioneer[:, 161:end] .== 0) && all(u0.stabilizer[:, 161:end] .== 0)
    @test u0.zref[:, 161:end] == u0.topo[:, 161:end]
    # wind, avalanches and storms never take the wall or promenade below their surface
    p = merge(prob.p, (; water_levels=[2.5, 3.0]))
    u = deepcopy(u0)
    for t in 1:100
        dune_step!(u, p, t, Random.Xoshiro(t))
        @test all(u.topo .>= u.hard)
    end
    # ...whereas the same surface, erodible, loses sand to the wind
    v = merge(deepcopy(u0), (; hard=fill(typemin(Int32), size(u0.topo))))
    for t in 1:100
        dune_step!(v, p, t, Random.Xoshiro(t))
    end
    @test any(v.topo[:, 163:end] .< u0.topo[:, 163:end])
    # a state without a hard layer runs as before
    w = Base.structdiff(deepcopy(u0), NamedTuple{(:hard,)})
    @test dune_step!(w, p, 1, Random.Xoshiro(1)) isa Integer

    # the blowout breaches the foredune in the middle fifth of the rows
    prob = problem(:dune_formation; scenario=:blowout, alongshore=10, years=1, seed=1)
    crest = vec(maximum(prob.u0.topo; dims=2)) .* 0.02
    @test all(crest[5:6] .≈ 5) && all(crest[[1:4; 7:10]] .> 14.5)
    @test all(prob.u0.hard .== typemin(Int32))
    # the double ridge: plants on the embryo dune (above 6 m) as on the foredune
    prob = problem(:dune_formation; scenario=:double_ridge, alongshore=2, years=1, seed=1)
    cover = prob.u0.pioneer .+ prob.u0.stabilizer
    @test all(cover[:, 375:385] .== 1) && all(prob.u0.pioneer[:, 375:385] .== 1)

    # a landform of one's own: one profile or one per row; the run-up gradient from its beach
    own = problem(:dune_formation; landform=dune_landform(; beach=100), alongshore=3,
                  years=1,
                  seed=4)
    narrow = problem(:dune_formation; scenario=:narrow, alongshore=3, years=1, seed=4)
    @test own.u0.topo == narrow.u0.topo && own.u0.zref == narrow.u0.zref
    @test own.u0.pioneer == narrow.u0.pioneer && own.p.slope ≈ 2.5 / 100
    rows = [isodd(r) ? fore : double for r in 1:4]
    mixed = problem(:dune_formation; landform=rows, alongshore=4, years=1)
    @test mixed.u0.topo[2, :] == profile_heights(double, 680; h=0.02)
    @test_throws ArgumentError problem(:dune_formation; landform=rows, alongshore=3)
    @test_throws ArgumentError problem(:dune_formation; landform=[fore, wall], alongshore=2)
    flat = shoreface_segment(60) * flat_segment(20, 0; kind=:hinterland)
    @test_throws ArgumentError problem(:dune_formation; landform=flat)
    @test problem(:dune_formation; landform=flat, pmap=(; slope=0.02), alongshore=1,
                  years=1).p.slope == 0.02
end

@testset "dune wind" begin
    E = EcologicalCellularAutomata
    # the shear stress over a sinusoidal ridge (Kroy et al. 2002, Eq. 12; Schwämmle & Herrmann
    # 2005, Eq. 3): for h = H cos(kx x + ky y), τx = A H (kx²/k cos - B kx|kx|/k sin) and
    # τy = A H kx ky / k cos, with k = |(kx, ky)|
    M, N, H, A, B = 32, 64, 0.7, 4.0, 0.25
    for (a, b) in ((3, 2), (-2, 5), (4, 0))
        kx, ky = 2π * a / N, 2π * b / M
        k = hypot(kx, ky)
        φ = [kx * (x - 1) + ky * (y - 1) for y in 1:M, x in 1:N]
        τx, τy = E._dune_shear(H .* cos.(φ), A, B)
        @test τx ≈ A * H .* (kx^2 / k .* cos.(φ) .- B * kx * abs(kx) / k .* sin.(φ)) atol = 1e-12
        @test τy ≈ A * H * kx * ky / k .* cos.(φ) atol = 1e-12
    end
    # the wind is unchanged over a bed at its reference surface, and off by default
    p = merge(defaults(build_dune_formation()), (; speedup=4.0))
    @test defaults(build_dune_formation()).speedup == 0
    u = (; topo=fill(Int32(100), 6, 40), zref=fill(Int32(100), 6, 40))
    w = dune_wind(u, p)
    @test all(w.along .≈ 1) && all(abs.(w.across) .< 1e-12)
    # over a foredune breached in the middle rows the wind speeds up towards the crest, strongest
    # just upwind of it, and turns into the gap upwind of the crest and out of it in the lee
    prob = problem(:dune_formation; scenario=:blowout, alongshore=20, years=1, seed=1)
    w = dune_wind(prob.u0, merge(prob.p, (; speedup=4.0)))
    crest = argmax(prob.u0.topo[1, :])
    @test argmax(w.along[1, :]) < crest && maximum(w.along[1, :]) > 2
    @test w.along[1, crest - 45] < 1                    # slower at the foot of the dune
    @test w.across[9, crest - 20] > 0.5 && w.across[12, crest - 20] < -0.5
    @test w.across[9, crest + 20] < -0.5 && w.across[12, crest + 20] > 0.5
    @test w.across ≈ -reverse(w.across; dims=1)        # the gap is central: mirror symmetry

    # transport in a given wind: a bare dry cell gives up a slab with probability
    # min(1, pe0 along); a flying slab moves a row towards `across` with probability
    # |across| / along at each cell (here it lands at the first, L = 1 and pd0 = 1)
    q = merge(defaults(build_dune_formation()), (; L=1, pd0=1.0, rho_star=Inf))
    # all cells but the first of the middle row lie below the groundwater surface, so only that
    # one gives up slabs
    bed(m, n) = (topo=fill(Int32(100), m, n), pioneer=zeros(m, n), stabilizer=zeros(m, n),
                 balance=zeros(Int32, m, n),
                 zref=[(r, x) == (cld(m, 2), 1) ? Int32(0) : Int32(200)
                       for r in 1:m, x in 1:n])
    picked, moved, runs = 0, 0, 4000
    for s in 1:runs
        u = bed(3, 2)
        wind = (; along=fill(1.6, 3, 2), across=fill(0.4, 3, 2))
        E._dune_transport!(u, q, Random.Xoshiro(s), falses(3, 2), fill(3, 3), wind)
        picked += u.topo[2, 1] < 100
        moved += u.topo[2, 1] < 100 && u.topo[3, 2] > 100
    end
    @test close_to(picked / runs, 0.8, runs) && close_to(moved / picked, 0.25, picked)
    u = bed(1, 2)                                       # at the domain's edge it stays
    E._dune_transport!(u, merge(q, (; pe0=1.0)), Random.Xoshiro(1), falses(1, 2), [3],
                       (; along=ones(1, 2), across=fill(5.0, 1, 2)))
    @test u.topo == Int32[99 101]
end

@testset "tuta_absoluta_asia" begin
    # the rule (Guimapi et al. 2016, 2020): a land cell next to an exposed or invaded cell becomes
    # invaded in a month with T >= 22 and RH >= 55, else exposed; exposed cells become invaded in a
    # suitable month; both pass the invasion on; sea cells never change
    model = build_tuta_absoluta_asia()
    n = 7
    land = trues(n, n)
    land[:, 6] .= false                   # a strait one cell wide
    state = fill(Int8(0), n, n)
    state[4, 1] = 2
    reached = fill(Int16(-1), n, n)
    invaded = fill(Int16(-1), n, n)
    reached[4, 1] = invaded[4, 1] = 0
    T = fill(25.0, n, n, 2)
    H = fill(60.0, n, n, 2)
    H[:, 3, 1] .= 40.0                    # column 3 too dry in month 1 only
    prob = CAProblem(model, (; state, land, reached, invaded), 4, (;); space=Grid(n, n),
                     covariates=(; temperature=T, humidity=H))
    us = solve(prob).u
    @test us[2].state[4, 2] == 2 && us[2].invaded[4, 2] == 1 && us[2].state[4, 3] == 0
    @test us[3].state[4, 3] == 2 && us[3].reached[4, 3] == 2      # month 2 is humid enough
    @test us[2].state[1, 1] == 0 && us[2].state[3, 2] == 2         # Moore neighbours, one ring a month
    @test all(==(0), us[end].state[:, 6]) && all(==(-1), us[end].reached[:, 6])
    @test all(==(0), us[end].state[:, 7])                          # the strait is not crossed
    # a dry month exposes, and the exposed cells carry the invasion on
    H2 = copy(H)
    H2[:, 2:3, :] .= 40.0
    u = solve(remake(prob; covariates=(; temperature=T, humidity=H2))).u
    @test u[2].state[4, 2] == 1 && u[3].state[4, 3] == 1 && u[4].state[4, 4] == 2
    @test u[4].invaded[4, 4] == 3 && u[end].invaded[4, 2] == -1
    H3 = copy(H)
    H3[:, 2, 1] .= 40.0                                            # dry in month 1, humid in month 2
    u = solve(remake(prob; covariates=(; temperature=T, humidity=H3))).u
    @test u[2].state[4, 2] == 1 && u[3].state[4, 2] == 2 && u[3].reached[4, 2] == 1 &&
          u[3].invaded[4, 2] == 2
    # a radius of 2 crosses the strait
    prob2 = CAProblem(build_tuta_absoluta_asia(; radius=2),
                      (; state, land, reached, invaded), 4,
                      (;); space=Grid(n, n), covariates=(; temperature=T, humidity=H))
    @test solve(prob2).u[end].state[4, 7] == 2
    @test_throws ArgumentError build_tuta_absoluta_asia(; radius=0)

    # the bundled climate: CRU CL 2.0 on 1/3-degree cells of Asian land
    c = tuta_climate()
    @test size(c.temperature) == (267, 462, 12) && count(c.land) == 49628
    @test all(isnan, c.temperature[.!c.land, :]) && all(!isnan, c.humidity[c.land, :])
    @test all(0 .<= c.humidity[c.land, :] .<= 100)
    I = tuta_cell(38.4, 27.1)                                      # Izmir: hot dry summers
    @test c.temperature[I, 7] > 22 && c.humidity[I, 7] < 55 && c.temperature[I, 1] < 12
    @test c.lat[I[1]] ≈ 78 - (I[1] - 0.5) / 3 && c.lon[I[2]] ≈ 26 + (I[2] - 0.5) / 3
    @test_throws ArgumentError tuta_cell(0.0, 0.0)
    @test_throws ArgumentError tuta_cell(10.0, 160.0)              # the open Pacific
    # the setup: the origin invaded, and the months starting at `start_month`
    prob = problem(:tuta_absoluta_asia; years=1, start_month=7)
    @test prob.u0.state[I] == 2 && count(==(2), prob.u0.state) == 1 &&
          prob.u0.reached[I] == 0
    @test isequal(prob.covariates.temperature[:, :, 1], c.temperature[:, :, 7])
    @test tuta_arrival(prob.u0, 38.4, 27.1) == (; reached=0.0, invaded=0.0)
    @test tuta_arrival(prob.u0, 28.6, 77.2) == (; reached=nothing, invaded=nothing)
    @test_throws ArgumentError problem(:tuta_absoluta_asia; start_month=13)
end

@testset "rootworm_rotation" begin
    # the crop probabilities keep the shares of maize (a) and of continuous maize among maize (b):
    # iterate the chain on years in maize (0, 1, 2, ...) to its steady state
    for policy in 0:3, (a, b, Q) in ((0.4, 0.5, 1.0), (0.3, 0.6, 0.5), (0.5, 0.45, 2.0))
        pr = rootworm_probabilities(a, b, Q, policy)
        pr === nothing && continue
        @test pr.pF ≈ Q * pr.pC
        cap = policy == 0 ? 40 : policy + 2
        x = zeros(cap + 1)
        x[1] = 1
        for _ in 1:2000
            y = zeros(cap + 1)
            y[2] += pr.pE * x[1]
            y[1] += (1 - pr.pE) * x[1]
            y[3] += pr.pF * x[2]
            y[1] += (1 - pr.pF) * x[2]
            for s in 3:(cap + 1)
                stay = s == cap + 1 ? (policy == 0 ? pr.pC : 0.0) : pr.pC
                s < cap + 1 && (y[s + 1] += stay * x[s])
                s == cap + 1 && (y[s] += stay * x[s])
                y[1] += (1 - stay) * x[s]
            end
            x = y
        end
        @test 1 - x[1] ≈ a atol = 1e-6
        @test sum(x[3:end]) / (1 - x[1]) ≈ b atol = 1e-6
    end
    @test rootworm_probabilities(0.6, 0.1, 1.0, 0) === nothing        # pE > 1
    @test_throws ArgumentError rootworm_probabilities(0.4, 0.5, 1.0, 4)

    # one year of crop change: continuous maize stays with probability pC, first-year with pF,
    # another crop is sown with maize with pE; mandatory rotation ends maize after policy + 2 years
    model = build_rootworm_rotation()
    pr = rootworm_probabilities(0.4, 0.5, 1.0, 1)
    p = merge(defaults(model), pr, (; policy=1))
    n = 60
    for (c0, q, c1) in ((2, pr.pC, 3), (1, pr.pF, 2), (0, pr.pE, 1))
        u0 = (; crop=fill(Int8(c0), n, n), pop=zeros(n, n))
        u = solve(CAProblem(model, u0, 1, p; space=Grid(n, n; boundary=Periodic()),
                            seed=c0)).u[end]
        @test close_to(count(==(c1), u.crop) / n^2, q, n^2) && all(in((0, c1)), u.crop)
    end
    u0 = (; crop=fill(Int8(3), n, n), pop=zeros(n, n))
    u = solve(CAProblem(model, u0, 1, p; space=Grid(n, n; boundary=Periodic()))).u[end]
    @test all(==(0), u.crop)                                            # rotated after 3 years

    # the dispersal range in fields, and growth and dispersal without noise
    @test rootworm_rings((; range=2000.0, fieldsize=15.0)) == 5
    @test rootworm_rings((; range=1000.0, fieldsize=100.0)) == 1
    @test rootworm_rings((; range=100.0, fieldsize=100.0)) == 1
    q = merge(defaults(model),
              (; growth_sd=0.0, emigration_sd=0.0, year_effect=0.0,
               range=2000.0, fieldsize=100.0))     # 2 rings
    k = rootworm_rings(q)
    @test k == 2
    m = 11
    crop = fill(Int8(1), m, m)
    crop[6, 6] = 2
    crop[1, :] .= 0
    pop = zeros(m, m)
    pop[6, 6] = 25.0
    pop[3, 3] = 7.0                          # first-year maize: lost when the year's growth comes
    u = (; crop=copy(crop), pop=copy(pop))
    rootworm_step!(u, q, Random.Xoshiro(1))
    grown = 25 * ((1 - 4) / (40 - 10) * (25 - 10) + 4)                  # rate 2.5 between the limits
    out = 0.6 * grown                                                    # above 10 adults per plant
    W = 8 * 2 + 16 * 1                                                   # all maize within 2 fields
    @test u.pop[6, 6] ≈ grown - out && u.pop[5, 6] ≈ out * 2 / W && u.pop[4, 6] ≈ out / W
    @test u.pop[3, 3] == 0 && u.pop[6, 9] == 0 && sum(u.pop) ≈ grown
    # emigration at low density, and none from a field with no maize within range
    u = (; crop=fill(Int8(0), m, m), pop=zeros(m, m))
    u.crop[6, 6] = 2
    u.pop[6, 6] = 0.3
    rootworm_step!(u, q, Random.Xoshiro(1); grow=false)
    @test u.pop[6, 6] == 0.3
    u.crop[6, 7] = 1
    rootworm_step!(u, q, Random.Xoshiro(1); grow=false)
    @test u.pop[6, 7] ≈ 0.2 * 0.3 && u.pop[6, 6] ≈ 0.8 * 0.3

    # the setup: rootworms only in continuous maize; risk counts maize fields above the threshold
    prob = problem(:rootworm_rotation; n=40, seed=3)
    @test all(prob.u0.pop[prob.u0.crop .== 0] .== 0)
    @test 0.3 < count(>(0), prob.u0.crop) / 40^2 < 0.5
    r = rootworm_risk((; crop=Int8[2 1; 0 2], pop=[2.0 0.5; 5.0 0.2]))
    @test r.all ≈ 1 / 3 && r.first_year == 0 && r.continuous ≈ 1 / 2
    sol = solve(problem(:rootworm_rotation; n=20, years=3, seed=1))
    @test rootworm_output(sol) ≈
          (rootworm_risk(sol.u[2]).all + rootworm_risk(sol.u[3]).all) / 2
    @test_throws ArgumentError problem(:rootworm_rotation; a=0.6, rotated=0.9)
    @test_throws ArgumentError problem(:rootworm_rotation; rotated=1.0)
end

@testset "resistance_regional" begin
    model = build_resistance_regional()
    p = defaults(model)
    field(N, q; sprayed=true, K=2e7, n=1) = (; N=fill(N, n, n), q=fill(q, n, n),
                                             sprayed=fill(sprayed, n, n),
                                             rotated=fill(false, n, n), K=fill(K, n, n),
                                             neighbours=fill(Int8(0), n, n))
    step1(u, q) = solve(CAProblem(model, u, 1, q; space=Grid(size(u.N)...))).u[end]
    # one field (Peck & Ellner 1997, Fig. 1): logistic growth N + floor(r (1 - N/K) N), Hardy-Weinberg
    # genotypes, the cost of resistance (1, 1 - sh, 1 - s), then a spray above delta Ks killing rho of
    # susceptible homozygotes and rho (1 - h) of heterozygotes
    q = merge(p, (; r=1.5, s=0.005, h=0.05, rho=0.75, delta=0.45, Ks=2e7))
    N0, q0 = 1e7, 0.01
    N = N0 + floor(1.5 * (1 - N0 / 2e7) * N0)
    GG, Gg, gg = (1 - q0)^2 * N, 2 * (1 - q0) * q0 * N, q0^2 * N
    Gg *= 1 - 0.005 * 0.05
    gg *= 1 - 0.005
    sprayGG, sprayGg = GG * 0.25, Gg * (1 - 0.75 * 0.95)
    v = step1(field(N0, q0), q)
    @test v.N[1] ≈ sprayGG + sprayGg + gg
    @test v.q[1] ≈ (gg + sprayGg / 2) / (sprayGG + sprayGg + gg)
    v = step1(field(N0, q0; sprayed=false), q)                         # a refuge is never sprayed
    @test v.N[1] ≈ GG + Gg + gg
    v = step1(field(N0, q0), merge(q, (; delta=0.9)))                  # below the threshold
    @test v.N[1] ≈ GG + Gg + gg
    @test step1(field(3e7, q0), merge(q, (; r=3.5))).N[1] == 0        # overshoot: no pests left

    # migration on a 3 x 3 region with no growth or selection: a share M leaves, split among the
    # neighbours in the lattice (von Neumann), carrying its field's allele frequency
    u = field(1e6, 0.0; sprayed=false, n=3)
    u.N[2, 2] = 4e6
    u.q[2, 2] = 0.5
    u.neighbours .= Int8[2 3 2; 3 4 3; 2 3 2]
    still = merge(p, (; r=0.0, s=0.0, rho=0.0, M=0.1))
    v = step1(u, still)
    @test v.N[2, 2] ≈ 4e6 - 4e5 + 4 * floor(1e5 / 3)
    @test v.N[1, 2] ≈ 1e6 - 3 * floor(1e5 / 3) + 1e5 + 2 * floor(1e5 / 2)
    @test v.q[1, 2] ≈ 0.5 * 1e5 / v.N[1, 2]
    @test sum(v.N) == sum(u.N)                                         # migrants are conserved
    # read literally, the source's equations send M N to each neighbour
    v = step1(u, merge(still, (; migration=:per_neighbour)))
    @test v.N[2, 2] ≈ 4e6 - 4 * 4e5 + 4 * 1e5
    # rotated fields are emptied every other generation
    w = (; u..., rotated=fill(true, 3, 3))
    us = solve(CAProblem(model, w, 3, still; space=Grid(3, 3))).u
    @test all(>(0), us[2].N) && all(==(0), us[3].N) && all(>(0), us[4].N) == false
    @test_throws ArgumentError build_resistance_regional(; neighbourhood=:hex)

    # the setup: a share lambda treated, carrying capacities, neighbours in the lattice
    prob = problem(:resistance_regional; n=20, lambda=0.55, kappa=0.5, N0=1e7, k=2, psi=0.4)
    u0 = prob.u0
    @test count(u0.sprayed) == 220 && count(u0.rotated) == 160
    @test all(u0.K[u0.sprayed] .== 2e7) && all(u0.K[.!u0.sprayed] .== 1e7)
    @test u0.neighbours[1, 1] == 2 && u0.neighbours[1, 5] == 3 && u0.neighbours[5, 5] == 4
    @test all(9e6 .<= u0.N .<= 1.1e7) && all(==(0.005), u0.q)
    @test problem(:resistance_regional; n=5, neighbourhood=:moore).u0.neighbours[3, 3] == 8
    @test_throws ArgumentError problem(:resistance_regional; lambda=1.5)
    @test_throws ArgumentError problem(:resistance_regional; pmap=(; migration=:all))

    # the measures: mean and variance over fields, Moran's I with edge neighbours
    s = resistance_summary((; q=[1.0 0.0; 0.0 1.0]))
    @test s.mean == 0.5 && s.variance == 0.25 && s.morans_i ≈ -1
    @test isnan(resistance_summary((; q=fill(0.3, 4, 4))).morans_i)
    sol = solve(problem(:resistance_regional; n=10, generations=150, pmap=(; M=0.5),
                        seed=2))
    t = resistance_time(sol)
    @test t !== nothing &&
          resistance_summary(sol.u[t + 1]).mean > 0.5 >=
          resistance_summary(sol.u[t]).mean
end

@testset "bt_crop_three_layer" begin
    model = build_bt_crop_three_layer()
    p = defaults(model)
    n = 100
    cells = n * n
    lattice(; pest=Int8(0), enemy=false, plant=1.0, pest_age=Int8(0), pest_hunger=Int16(0),
            pest_genome=UInt32(0), enemy_hunger=Int16(0)) = (; plant=fill(plant, n, n),
                                                             pest=fill(pest, n, n),
                                                             pest_age=fill(pest_age, n, n),
                                                             pest_hunger=fill(pest_hunger,
                                                                              n, n),
                                                             pest_genome=fill(pest_genome,
                                                                              n, n),
                                                             enemy=fill(enemy, n, n),
                                                             enemy_age=zeros(Int8, n, n),
                                                             enemy_hunger=fill(enemy_hunger,
                                                                               n, n),
                                                             enemy_genome=zeros(UInt32, n,
                                                                                n))
    step1(u, q; seed=1) = solve(CAProblem(model, u, 1, q; space=Torus(n, n), seed),
                                Margolus()).u[end]
    # a full layer cannot move: an enemy above every B eats it with probability P_DB (Yang et al.
    # 2009, Sec. 2.2), and its hunger falls to 0 or rises by U_D
    v = step1(lattice(; pest=Int8(1), enemy=true), merge(p, (; P_DB=0.3)))
    @test close_to(count(==(0), v.pest) / cells, 0.3, cells)
    @test all(v.enemy_hunger[v.pest .== 0] .== 0) && all(v.enemy_hunger[v.pest .== 1] .== 1)
    v = step1(lattice(; pest=Int8(2), enemy=true), merge(p, (; P_DC=0.8)))
    @test close_to(count(==(0), v.pest) / cells, 0.8, cells)
    # a pest above enough plant eats d with probability P (B) and its hunger falls; the plant
    # regrows by `regrowth` up to its maximum; a pest that does not eat grows hungrier by U
    q = merge(p, (; P_BA=0.35, d_B=0.3, regrowth=0.1, h_BA=2, hmax_B=5))
    v = step1(lattice(; pest=Int8(1), plant=0.9, pest_hunger=Int16(2)), q)
    fed = v.plant .≈ 0.7
    @test close_to(count(fed) / cells, 0.35, cells)
    @test all(v.plant[.!fed] .== 1.0) && all(v.pest_hunger[fed] .== 0) &&
          all(v.pest_hunger[.!fed] .== 3)
    v = step1(lattice(; pest=Int8(1), plant=0.25), q)                   # too little to eat
    @test all(v.plant .≈ 0.35) && all(v.pest_hunger .== 1)
    # death: hunger above hmax; inherited diseases expressed up to the new age reaching T; the
    # maximum age
    starving = merge(p, (; P_BA=0.0, hmax_B=3))
    @test all(==(0), step1(lattice(; pest=Int8(1), pest_hunger=Int16(3)), starving).pest)
    @test all(==(1), step1(lattice(; pest=Int8(1), pest_hunger=Int16(2)), starving).pest)
    sick = merge(p, (; P_BA=1.0, T_B=3, R_B=40))
    @test all(==(0),
              step1(lattice(; pest=Int8(1), pest_age=Int8(2), pest_genome=UInt32(0b111)),
                    sick).pest)
    v = step1(lattice(; pest=Int8(1), pest_age=Int8(1), pest_genome=UInt32(0b111)), sick)
    @test all(==(1), v.pest) && all(==(2), v.pest_age)
    @test all(==(0),
              step1(lattice(; pest=Int8(1), pest_age=Int8(31)), merge(sick, (; R_B=99))).pest)

    # pairs: a lone immature pest moves to the cell its domino pairs it with (right or below, each
    # with probability 1/2 at time 0); a mature one in a birth cycle stays and leaves an offspring
    # there, with its genome and one more deleterious mutation
    lone = (; plant=fill(1.0, 4, 4), pest=zeros(Int8, 4, 4), pest_age=zeros(Int8, 4, 4),
            pest_hunger=zeros(Int16, 4, 4), pest_genome=zeros(UInt32, 4, 4),
            enemy=falses(4, 4), enemy_age=zeros(Int8, 4, 4),
            enemy_hunger=zeros(Int16, 4, 4),
            enemy_genome=zeros(UInt32, 4, 4))
    lone.pest[1, 1] = 1
    fed_all = merge(p, (; P_BA=1.0, R_B=3, r_B=2))
    right = 0
    runs = 2000
    for seed in 1:runs
        v = solve(CAProblem(model, lone, 1, fed_all; space=Torus(4, 4), seed), Margolus()).u[end]
        @test count(!=(0), v.pest) == 1 && (v.pest[1, 2] == 1 || v.pest[2, 1] == 1)
        right += v.pest[1, 2] == 1
    end
    @test close_to(right / runs, 0.5, runs)
    w = (; lone..., pest_age=fill(Int8(5), 4, 4), pest_genome=fill(UInt32(0b1), 4, 4))
    v = solve(CAProblem(model, w, 1, fed_all; space=Torus(4, 4), seed=3), Margolus()).u[end]
    child = findfirst(==(1), v.pest .* (v.pest_age .== 1))
    @test count(==(1), v.pest) == 2 && v.pest[1, 1] == 1 && v.pest_age[1, 1] == 6
    @test child in (CartesianIndex(1, 2), CartesianIndex(2, 1))
    @test count_ones(v.pest_genome[child]) in (1, 2) && v.pest_genome[child] & 0b1 == 1
    w = (; w..., pest_age=fill(Int8(4), 4, 4))                       # age 4: not a birth cycle
    v = solve(CAProblem(model, w, 1, fed_all; space=Torus(4, 4), seed=3), Margolus()).u[end]
    @test count(==(1), v.pest) == 1 && v.pest[1, 1] == 0

    # without births, predation or starvation, moves conserve every population
    still = merge(p,
                  (; R_B=99, R_C=99, R_D=99, P_DB=0.0, P_DC=0.0, hmax_B=10^4, hmax_C=10^4,
                   hmax_D=10^4, max_age=32))
    u0 = problem(:bt_crop_three_layer; n=20, seed=4, pmap=still).u0
    u0 = (; u0..., pest_age=zeros(Int8, 20, 20), enemy_age=zeros(Int8, 20, 20))
    sol = solve(CAProblem(model, u0, 25, still; space=Torus(20, 20), seed=4), Margolus())
    s0 = bt_crop_summary(sol.u[1])
    @test all(u -> (bt_crop_summary(u).B, bt_crop_summary(u).C, bt_crop_summary(u).D) ==
                   (s0.B, s0.C, s0.D), sol.u)
    @test sol.u[end].pest != sol.u[1].pest

    # the setup: plants at their maximum, enemies and pests at the given occupancy, B and C equally
    u0 = problem(:bt_crop_three_layer; seed=2).u0
    s = bt_crop_summary(u0)
    @test s.plant == 1 && close_to(s.D / cells, 0.5, cells) &&
          close_to((s.B + s.C) / cells, 0.5, cells) &&
          close_to(s.B / (s.B + s.C), 0.5, s.B + s.C)
    @test all(u0.pest_age[u0.pest .== 1] .< p.R_B) && all(u0.enemy_age[u0.enemy] .< p.R_D)
    @test_throws ArgumentError problem(:bt_crop_three_layer; n=51)
    @test_throws ArgumentError problem(:bt_crop_three_layer; occupancy=1.5)
    @test_throws ArgumentError problem(:bt_crop_three_layer; pmap=(; max_age=40))
    # without feeding on the plant, pests starve before they mature, and their enemies follow
    sol = solve(problem(:bt_crop_three_layer; n=40, cycles=60, seed=1,
                        pmap=(; P_BA=0.0, P_CA=0.0)), Margolus())
    @test bt_crop_summary(sol.u[end]) == (; plant=1.0, B=0, C=0, D=0)
end

@testset "barrier_island" begin
    EZ = EcologicalCellularAutomata
    model = build_barrier_island()
    p = defaults(model)
    MHW = p.MHW
    # the shoreface of Table A1 from the wave climate and sediment (Eq. 15-18): 20.07 m, 0.00822,
    # 5926 m³/m/yr; the alongshore transport constant "∼ 0.06"
    sf = barrier_shoreface()
    @test round(sf.D; digits=2) == 20.07 && round(sf.S_eq; digits=5) == 0.00822 &&
          round(sf.k) == 5926
    @test (p.D_sf, p.S_sf_eq, p.k_sf) == (sf.D, sf.S_eq, sf.k)
    k = EZ._meeb_transport_k()
    @test 0.06 < k < 0.065
    # shoreline diffusivity (Eq. 20-22) at a shoreline angle of 0, the middle of the 181 angles:
    # with waves spread evenly (a = h = 1/2, E = 1/π per radian) the sum over wave angles of 1°
    # of E Ψ, Ψ(φ) = cos^(1/5) φ (cos² φ - 6/5 sin² φ), times k/D Hs^(12/5) T^(1/5) a year
    Ψ(φ) = cos(φ)^0.2 * (cos(φ)^2 - 1.2 * sin(φ)^2)
    even = barrier_diffusivity(; a=0.5, h=0.5)
    by_hand = k / sf.D * 0.98^2.4 * 6.6^0.2 * 365 * 24 * 3600 *
              sum(Ψ(deg2rad(d)) / π * π / 180 for d in -89:90)
    @test length(even) == 181 && isapprox(even[91], by_hand; rtol=1e-9)
    @test barrier_diffusivity()[91] > 0

    # aeolian transport on a bare plateau 2 m high, wind landward: each cell gives up a slab with
    # probability pe0 = 0.10, and with smoothed cover ρ = 0.05 with pe0 (1 - ρ/ρq0) = 0.05
    # (deposition off, so every slab leaves the domain)
    plateau(n; z=2.0, ρ=0.0) = (topo=fill(z, 1, n), veg_eff=fill(ρ, 1, n))
    land = merge(p, (; wind_rose=(1.0, 0.0, 0.0, 0.0)))
    function eroded(ρ; runs=200)
        q = merge(land, (; pd0=0.0, pd1=0.0))
        hits = 0
        for s in 1:runs
            hits += count(==(-1),
                          EZ._meeb_aeolian(plateau(60; ρ), q, Random.Xoshiro(s), MHW, 1))
        end
        return hits / (59runs)
    end
    @test close_to(eroded(0.0), 0.10, 59 * 200)
    @test close_to(eroded(0.05), 0.05, 59 * 200)
    @test eroded(0.1) == 0
    # with certain pick-up and landing, every slab lands one hop downwind, the hop drawn uniformly
    # from 3 to 7 cells each week: the first `hop` cells lose a slab
    hops = zeros(Int, 7)
    for s in 1:1000
        Δ = EZ._meeb_aeolian(plateau(40), merge(land, (; pe0=1.0, pd0=1.0, pd1=1.0)),
                             Random.Xoshiro(s), MHW, 1)
        hops[findfirst(!=(-1), vec(Δ)) - 1] += 1
    end
    @test sum(hops[3:7]) == 1000 && all(close_to(hops[h] / 1000, 0.2, 1000) for h in 3:7)
    # part-way through a hop a slab lands only where ρ ≥ ρv = 0.35, with pd0 + (pd1 - pd0) ρ:
    # 0.38 at ρ = 0.5 two cells downwind of its source (hops are at least 3)
    landed = 0
    for s in 1:2000
        u = plateau(30; ρ=0.2)
        u.veg_eff[1, 5] = 0.0
        u.veg_eff[1, 7] = 0.5
        Δ = EZ._meeb_aeolian(u, merge(land, (; pe0=1.0)), Random.Xoshiro(s), MHW, 1)
        landed += Δ[1, 7] == 1
    end
    @test close_to(landed / 2000, 0.22 + (0.54 - 0.22) * 0.5, 2000)
    # shadow: in the 12° lee of a cell 2 m above a plateau, downwind of the wind
    row = fill(1.0, 1, 40)
    row[1, 10] = 3.0
    lee = 2.0 / tand(12)
    S = EZ._meeb_shadow(row, 1, 12.0)
    @test all(S[1, 11:(10 + floor(Int, lee))]) && !S[1, 11 + floor(Int, lee)] &&
          !any(S[1, 1:10])
    @test all(EZ._meeb_shadow(row, 3, 12.0)[1, 1:9])

    # avalanches (as in DuBeVeg): 20° bare, 18 slabs of 0.02 m orthogonally, 25 along a diagonal;
    # only the steepest drop counts: 24 along the diagonal holds, 26 sends one slab down it (the
    # third column, also above MHW, closes the window in which avalanches act)
    @test (floor(Int, tand(20) / 0.02), floor(Int, sqrt(2) * tand(20) / 0.02)) == (18, 25)
    slabs(A) = (topo=hcat(A, A[:, end]) .* 0.02, veg=zeros(2, 3))
    u = slabs([100 80; 80 76])
    EZ._meeb_avalanche!(u, p, Random.Xoshiro(1), MHW)
    @test u.topo[:, 1:2] ≈ [100 80; 80 76] .* 0.02
    u = slabs([100 80; 80 74])
    EZ._meeb_avalanche!(u, p, Random.Xoshiro(1), MHW)
    @test u.topo[:, 1:2] ≈ [99 80; 80 75] .* 0.02

    # overwash (Eq. 4-8): water 0.5 m above a crest enters at √(2 g R) R per metre of crest
    # (g = 9.8, in m³/h), loses Rin at each cross-shore step, flows down a slope of 0.01 to the one
    # landward neighbour and carries Kow (Q (S + Cs))^m of sand
    E = [0.0 1.0 [2.0 - 0.01 * c for c in 0:9]...]
    _, net = EZ._meeb_overwash(E, [3], MHW, [2.5], zeros(1, 12), zeros(1, 12), p)
    Q0 = sqrt(2 * 9.8 * 0.5) * 0.5 * 3600 - p.Rin
    q(Q) = p.Kow * (Q * (0.01 + p.Cs))^p.mm
    @test net[1, 3] ≈ -q(Q0) && net[1, 4] ≈ q(Q0) - q(Q0 - p.Rin)
    # sharing among the three landward neighbours (Eq. 6-7): downhill in proportion to S^n, n = 1/2;
    # all uphill in proportion to |S|^-n, reduced by 1 - |S|/SQlim, none beyond SQlim = 1.5
    @test all(EZ._meeb_share(100.0, 0.04, 0.01, 0.0, p, true) .≈
              (100 * 0.2 / 0.3, 100 * 0.1 / 0.3, 0.0))
    w = (0.3, 0.6, 1.6) .^ -0.5
    @test all(EZ._meeb_share(100.0, -0.3, -0.6, -1.6, p, true) .≈
              (100 * w[1] / sum(w) * 0.8, 100 * w[2] / sum(w) * 0.6, 0.0))
    # plants slow the flow leaving towards them (Eq. 23): Q (1 - Λ ρ)
    _, net2 = EZ._meeb_overwash(E, [3], MHW, [2.5], fill(0.5, 1, 12), zeros(1, 12), p)
    @test net2[1, 3] ≈ -q(Q0 * (1 - 0.02 * 0.5))

    # the beach (Eq. 2-3): on a beach at the equilibrium slope the flux (βeq - β)(Ztwl - Z)² is
    # zero up to the last wet cell, where the slope to the water level is gentler; the elevation
    # falls by its divergence times Te, so only the cells either side of it change
    beach = merge(p, (; beach_substeps=1))
    E = reshape([MHW + 0.022 * (x - 9) for x in 1:70], 1, 70)
    Rh = MHW + 0.022 * 30.5
    Δ, dV, wet = EZ._meeb_beach!(copy(E), [60], [10.0], MHW, [Rh], beach)
    f39 = (0.022 - 0.5 * (Rh - E[38])) * (Rh - E[39])^2
    @test Δ[1, 38] ≈ -f39 / 2 * p.Te && Δ[1, 40] ≈ f39 / 2 * p.Te &&
          count(!=(0), Δ) == 2 && abs(dV[1]) < 1e-15
    @test findall(vec(wet)) == 1:38

    # runup (Stockdon et al. 2006): setup 0.35 β √(H L0), swash √((0.75 β √(H L0))² +
    # (0.06 √(H L0))²), R2 = 1.1 (setup + swash / 2), L0 = g T² / 2π, on the still water level
    L0 = 9.8 * 8^2 / 2π
    R2 = 1.1 * (0.35 * 0.03 * sqrt(2 * L0) +
                sqrt((0.75 * 0.03 * sqrt(2 * L0))^2 + (0.06 * sqrt(2 * L0))^2) / 2)
    @test EZ._meeb_twl([2.0, 10.0, 1.0, 8.0], [0.03])[1] ≈ 1.0 + R2

    # growth (Fig. 1e): -1 at and beyond the outer vertices, 0 at the second and fourth, the peak
    # at the third, linear between
    v = p.vertices1
    @test [EZ._meeb_growth(x, v, 0.2) for x in (v[1] - 1, v[1], v[2], v[3], v[4], v[5])] ≈
          [-1, -1, 0, 0.2, 0, -1]
    @test EZ._meeb_growth((v[2] + v[3]) / 2, v, 0.2) ≈ 0.1

    # shoreline (Eq. 13-14): at the equilibrium shoreface the toe moves 2 RSLR / S a step (1/25
    # year) and the shoreline not at all; overwash of 10 m³ per metre moves it 2 qow / D landward
    q6 = merge(p, (; RSLR=0.006))
    xs, xt = fill(100.0, 25), fill(100.0 - p.D_sf / p.S_sf_eq, 25)
    EZ._meeb_cross_shore!(xs, xt, zeros(25), zeros(25), q6)
    @test all(xt .≈ 100.0 - p.D_sf / p.S_sf_eq + 2 * 0.006 / 25 / p.S_sf_eq) &&
          all(xs .≈ 100)
    EZ._meeb_cross_shore!(xs, xt, fill(10.0, 25), zeros(25), merge(p, (; k_sf=0.0)))
    @test all(xs .≈ 100 + 2 * 10 / p.D_sf)
    # alongshore diffusion (Eq. 19): a straight shoreline stays; a bump decays, its neighbours
    # rise, the end sections stay
    xs = fill(50.0, 125)
    EZ._meeb_alongshore!(xs, p.diffusivity, p, 0.04)
    @test all(xs .≈ 50)
    xs[51:75] .= 52.0
    EZ._meeb_alongshore!(xs, p.diffusivity, p, 0.04)
    @test xs[51] < 52 && xs[26] > 50 && xs[76] > 50 && xs[1] == 50 && xs[101] == 50
    # regrading a shoreline eroded by a cell: the new shoreline cell drops to MHW at most, and a
    # linear shoreface of slope S runs seaward of it
    topo = reshape([MHW + 0.02 * (c - 61) + 0.06 for c in 1:80], 1, 80)
    EZ._meeb_regrade!(topo, [62.4], [61.3], MHW, [0.01], p)
    @test topo[1, 62] ≈ MHW && all(topo[1, c] ≈ MHW + (c - 62) * 0.01 for c in 1:61)

    # the foredune crest of the idealised barrier is its highest cell, with no gaps
    prob = problem_barrier_island(; alongshore=30, years=1, seed=1)
    u = deepcopy(prob.u0)
    crest, ok = barrier_crest(u.topo, MHW)
    @test all(crest .== argmax(u.topo[1, :])) && all(ok)
    @test size(u.topo) == (30, length(barrier_profile())) && all(u.shoreline .== 61)
    m = barrier_measures(u, MHW)
    @test m.crest ≈ maximum(u.topo) - MHW && m.gaps == 0
    gapped = problem_barrier_island(; scenario=:gap, alongshore=30, years=1).u0.topo
    @test maximum(gapped[13:18, :]) < MHW + 1.7 < maximum(gapped[1, :])

    # high-water events: one in each fortnight with the fortnight's probability (1/3 in the first)
    small = barrier_profile(; offshore=5, beach=20, dune=3.0, dune_width=10, interior=15,
                            bay=10)
    tiny() = (topo=repeat(small', 5), spec1=zeros(5, length(small)),
              spec2=zeros(5, length(small)), veg=zeros(5, length(small)),
              veg_eff=zeros(5, length(small)), balance=zeros(5, length(small)),
              change=zeros(5, length(small)), shoreline=fill(6.0, 5, length(small)),
              toe=fill(6.0 - 2441, 5, length(small)),
              inundated=zeros(Int, 5, length(small)))
    stormy = merge(p, (; storms=[2.0 1.0 1.0 8.0], pe0=0.0))
    events = count(s -> barrier_step!(tiny(), stormy, 1, Random.Xoshiro(s),
                                      Random.Xoshiro(10^6 + s)).event, 1:1500)
    @test close_to(events / 1500, 0.333, 1500)
    @test !barrier_step!(tiny(), stormy, 2, Random.Xoshiro(1), Random.Xoshiro(2)).event
    # a flood over a low barrier kills the plants wherever it reaches (smoothed over 2 m)
    u = tiny()
    u.topo .= repeat(barrier_profile(; offshore=5, beach=20, dune=0.5, dune_width=10,
                                     interior=15, bay=10)', 5)
    u.spec1 .= 0.5
    flood = merge(stormy, (; storms=[2.0 2.0 2.5 8.0], pstorm=ones(25)))
    r = barrier_step!(u, flood, 1, Random.Xoshiro(1), Random.Xoshiro(2))
    @test r.event && r.twl > 2.5 && all(u.spec1[u.inundated .> 0] .== 0) &&
          count(>(0), u.inundated) > length(small)

    # the problem: reproducible, the balance kept over the year or cleared each week
    sol1 = solve(problem_barrier_island(; alongshore=12, years=1, seed=4);
                 save_everystep=false)
    sol2 = solve(problem_barrier_island(; alongshore=12, years=1, seed=4);
                 save_everystep=false)
    @test sol1.u[end].topo == sol2.u[end].topo && !any(isnan, sol1.u[end].topo)
    u = deepcopy(prob.u0)
    barrier_step!(u, merge(prob.p, (; pstorm=zeros(25))), 2, Random.Xoshiro(1),
                  Random.Xoshiro(2))
    @test any(!=(0), u.balance)
    barrier_step!(u, merge(prob.p, (; pstorm=zeros(25), balance_period=:week)), 3,
                  Random.Xoshiro(1), Random.Xoshiro(2))
    @test all(==(0), u.balance)
    @test problem_barrier_island(; storms=:hindcast, alongshore=2).p.storm_series ==
          barrier_storm_series()
    @test size(barrier_storms()) == (9999, 4) && size(barrier_storm_series()) == (222, 4)
    @test_throws ArgumentError problem_barrier_island(; scenario=:cliff)
    @test_throws ArgumentError problem_barrier_island(; storms=:none)
    # the Brier skill score (Eq. 24): 1 for a perfect hindcast, 0 for predicting no change
    obs, base = [1.0 2.0; 3.0 4.0], zeros(2, 2)
    mask = trues(2, 2)
    @test brier_skill_score(obs, obs, base, mask) == 1 &&
          brier_skill_score(base, obs, base, mask) == 0
    # the Florence test site of Fig. 5a: 500 m of North Core Banks, the sea at the first column,
    # overwash where the barrier is low; one event changes the barrier where it floods and leaves
    # plants as they are
    site = barrier_florence_site()
    @test size(site.before) == size(site.after) == size(site.overwash) == (500, 701)
    @test all(site.before[:, 1] .< site.MHW) && maximum(site.before) > 5 &&
          0 < count(site.overwash) < length(site.overwash) && all(0 .<= site.spec1 .<= 1)
    @test (site.Rhigh, site.duration) == (3.32, 83)
    u = tiny()
    u.topo .= repeat(barrier_profile(; offshore=5, beach=20, dune=0.5, dune_width=10,
                                     interior=15, bay=10)', 5)
    u.spec1 .= 0.5
    r = barrier_event!(u, p, 2.0, 2)
    @test any(r.flooded) && any(!=(0), r.change) && all(u.spec1 .== 0.5) &&
          length(r.overwash) == 5
end

@testset "pitch_canker" begin
    E = EcologicalCellularAutomata
    # the gamma function, for the scale of a Weibull flight of given mean
    @test E._pc_gamma(2.0) ≈ 1 && E._pc_gamma(0.5) ≈ sqrt(π) &&
          E._pc_gamma(5 / 3) ≈ 0.9027452929509336
    @test isapprox(E._pc_lgamma(30.0), Float64(log(factorial(big(29)))); atol=1e-10) &&
          isapprox(E._pc_lgamma(7.5), log(E._pc_gamma(7.5)); atol=1e-12)

    # the binomial sampler against the exact distribution, by inversion and by BTRS
    rng = Random.Xoshiro(11)
    for (n, p) in ((20, 0.3), (1000, 0.4), (5000, 0.93))
        N = 200_000
        x = [E._pc_binomial(rng, n, p) for _ in 1:N]
        lp = zeros(n + 1)   # log pmf by recursion from k = 0
        lp[1] = n * log1p(-p)
        for k in 0:(n - 1)
            lp[k + 2] = lp[k + 1] + log((n - k) / (k + 1)) + log(p / (1 - p))
        end
        pmf = exp.(lp)
        μ = n * p
        @test abs(sum(x) / N - μ) <= 4 * sqrt(n * p * (1 - p) / N)
        for k in round(Int, μ - 2sqrt(μ * (1 - p))):round(Int, μ + 2sqrt(μ * (1 - p)))
            f = count(==(k), x) / N
            @test abs(f - pmf[k + 1]) <= 4 * sqrt(pmf[k + 1] * (1 - pmf[k + 1]) / N) + 1e-4
        end
    end

    # flight distances: a negative exponential of mean μ cut at M has mean
    # μ - M e^(-M/μ) / (1 - e^(-M/μ)); a Weibull has the given mean
    d = [E._pc_truncexp(rng, 280.0, 1300.0) for _ in 1:200_000]
    m = 280 - 1300 * exp(-1300 / 280) / (1 - exp(-1300 / 280))
    @test maximum(d) <= 1300 && abs(sum(d) / length(d) - m) <= 4 * 280 / sqrt(length(d))
    base = defaults(build_pitch_canker())
    pw = merge(base, (; insect_distribution=:weibull, insect_shape=1.5, insect_max=Inf))
    d = [E._pc_flight(rng, pw) for _ in 1:200_000]
    @test abs(sum(d) / length(d) - 1000) <= 4 * 700 / sqrt(length(d))
    # the landing kernel of an insect from a cell's midpoint against flown insects
    K, R = E._pc_insect_kernel(base)
    @test R == 5 && isapprox(sum(K), 1; atol=1e-6) &&
          K ≈ reverse(K; dims=1) ≈ permutedims(K)
    hits = zeros(Int, 2R + 1, 2R + 1)
    N = 400_000
    for _ in 1:N
        r = E._pc_truncexp(rng, 1000.0, 4500.0)
        θ = 2π * rand(rng)
        a = floor(Int, (500 + r * cos(θ)) / 1000)
        b = floor(Int, (500 + r * sin(θ)) / 1000)
        hits[a + R + 1, b + R + 1] += 1
    end
    @test all(abs.(hits ./ N .- K) .<= 4 .* sqrt.(K .* (1 .- K) ./ N) .+ 1e-4)

    # one year of the rule (no callback): two spreading years, 20,000 healthy trees in every cell
    # of a 5 × 5 grid with EI = 1, and at the centre 1000 trees infected this year, 500 a year ago
    # and 300 two years ago. The 300 die; 1000 + 500 spread, so with N_max = 50,000 the centre gains
    # 20,000 × 1500 / 50,000 = 600 infections (Eq. 4 with d = 0) and each orthogonal neighbour
    # 600 exp(-1000 / 280) (Eq. 3; its whole part, and one more with probability the rest), the
    # diagonals (1,414 m) none
    model = build_pitch_canker(; spreading_years=2)
    pr = defaults(model)
    function year(seed; t0=0, ei_now=1.0, q=pr)
        z() = zeros(Int32, 5, 5)
        u = (; hosts=fill(Int32(20_000), 5, 5), coverage=fill(0.4, 5, 5),
             ei_now=fill(ei_now, 5, 5), ei_future=ones(5, 5), new=z(), i1=z(), i2=z(),
             dead=z(), infected_year=z())
        u.new[3, 3], u.i1[3, 3], u.i2[3, 3], u.infected_year[3, 3] = 1000, 500, 300, 1
        return solve(CAProblem(model, u, (t0, t0 + 1), q; space=Grid(5, 5), seed)).u[end]
    end
    v = year(1)
    @test v.i1[3, 3] == 1000 && v.i2[3, 3] == 500 && v.dead[3, 3] == 300 &&
          v.new[3, 3] == 600 && v.hosts[3, 3] == 19_400 && v.new[2, 2] == 0 &&
          v.new[1, 3] == 0 && v.infected_year[2, 3] == (v.new[2, 3] > 0 ? 1 : 0)
    λ = 600 * exp(-1000 / 280)
    runs = 4000
    ns = [year(s).new[2, 3] for s in 1:runs]
    @test all(in((floor(Int, λ), floor(Int, λ) + 1)), ns)
    f = λ - floor(λ)
    @test abs(count(==(floor(Int, λ) + 1), ns) / runs - f) <= 4 * sqrt(f * (1 - f) / runs)
    # the index moves linearly from the current climate in the first year to the future one in
    # climate_year: in year 5 of 2005-2015, 0.4 of the way from 0 to 1
    v = year(1; t0=4, ei_now=0.0, q=merge(pr, (; climate_year=2015)))
    @test v.new[3, 3] == 240

    # seedling lots: 6000 lots from the middle of a 101 × 101 grid, up to 40 km, one dispatch each
    # landing within a metre, on hosts with EI = 1 and N_max so large that each infects. Proposals
    # are symmetric, so with coverage 0.2 on the left half and 1 on the right the left gets
    # 0.2 / 1.2 of the lots; on uniform coverage the distance is uniform on [0, 40 km]
    function lots(cov; seed=1)
        z() = zeros(Int32, 101, 101)
        u = (; hosts=fill(Int32(10_000_000), 101, 101), coverage=cov, ei_now=ones(101, 101),
             ei_future=ones(101, 101), new=z(), i1=z(), i2=z(), dead=z(), infected_year=z())
        q = merge(pr,
                  (; trees_per_ha=1e5, transport_max=40_000.0, dispatches=1,
                   spore_mean=1.0, spore_max=1.0,
                   entries=[(; i=51, j=51, first=1, last=1, lots=6000)]))
        return pitch_canker_spread!(u, q, 1, Random.Xoshiro(seed))
    end
    cov = [j < 51 ? 0.2 : j > 51 ? 1.0 : 0.0 for _ in 1:101, j in 1:101]
    u = lots(cov)
    left, right = sum(u.new[:, 1:50]), sum(u.new[:, 52:end])
    @test left + right == 6000 && abs(left / 6000 - 1 / 6) <= 4 * sqrt(5 / 36 / 6000)
    u = lots(ones(101, 101))
    r = [hypot(I[1] - 51, I[2] - 51) * 1000 for I in CartesianIndices(u.new)
         for _ in 1:u.new[I]]
    @test length(r) == 6000 &&
          abs(sum(r) / 6000 - 20_000) <= 4 * 40_000 / sqrt(12 * 6000) + 300

    # insects: 2000 spreading trees at the centre of a 21 × 21 grid where each try infects with
    # probability 0.1 (EI = 1, healthy trees a tenth of N_max, too many to run down), so an insect
    # infects within 10 tries with probability 1 - 0.9^10, at offset o with probability
    # K(o) (1 - 0.9^10); flown one by one and pooled
    function insects(method, seed)
        z() = zeros(Int32, 21, 21)
        u = (; hosts=fill(Int32(500_000), 21, 21), coverage=ones(21, 21),
             ei_now=ones(21, 21),
             ei_future=ones(21, 21), new=z(), i1=z(), i2=z(), dead=z(), infected_year=z())
        u.i1[11, 11] = 2000
        q = merge(pr, (; trees_per_ha=50_000.0, insect_method=method))
        return pitch_canker_spread!(u, q, 1, Random.Xoshiro(seed)).new
    end
    P = 1 - 0.9^10
    sd = sqrt(2000 * P * (1 - P))
    for method in (:each, :pooled)
        news = [insects(method, s) for s in 1:60]
        @test abs(sum(sum, news) / 60 - 2000P) <= 4 * sd / sqrt(60)
        centre = sum(n -> n[11, 11], news) / 60
        pc = K[R + 1, R + 1] * P
        @test abs(centre - 2000pc) <= 4 * sqrt(2000pc / 60) + 1
    end

    # the problem: reproducible, infection arriving from the entry points, trees dying after the
    # spreading years; the infected area is the host area of cells with any infected or dead tree
    prob = problem_pitch_canker(; dims=(60, 30), years=8, seed=3)
    s1 = solve(prob, Synchronous())
    s2 = solve(prob, Synchronous())
    @test s1.u[end] == s2.u[end] && sum(s1.u[end].dead) > 0 &&
          s1.u[2].infected_year[s1.u[2].infected_year .> 0] ==
          fill(1, count(>(0), s1.u[2].infected_year))
    u = (; coverage=[0.5 0.5; 1.0 0.0], new=Int32[1 0; 0 0], i1=Int32[0 0; 0 0],
         dead=Int32[0 0; 2 0])
    @test pitch_canker_infected_area(u) ≈ 100 * 1.5 / 2 &&
          pitch_canker_infected(u) == [1 0; 0 0]
    @test_throws ArgumentError build_pitch_canker(; spreading_years=0)
    @test_throws ArgumentError problem_pitch_canker(; coverage=zeros(5, 5))
    @test_throws ArgumentError problem_pitch_canker(; dims=(30, 10),
                                                    entries=[(; i=40, j=1, first=1, last=2,
                                                              lots=1)])
    @test_throws ArgumentError pitch_canker_spread!(year(1),
                                                    merge(pr, (; insect_method=:bad)),
                                                    1, Random.Xoshiro(1))
    @test_throws ArgumentError pitch_canker_landscape(Random.Xoshiro(1); coast=0.7,
                                                      inland=0.5)
end
