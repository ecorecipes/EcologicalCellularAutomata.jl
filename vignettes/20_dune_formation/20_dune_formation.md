# Coastal dunes
Simon Frost

- [Overview](#overview)
- [The model](#the-model)
  - [Water levels](#water-levels)
  - [Profiles from pieces](#profiles-from-pieces)
- [The same potential transport](#the-same-potential-transport)
- [Fifty years of dune growth](#fifty-years-of-dune-growth)
- [Deviations from the source](#deviations-from-the-source)

## Overview

Coastal foredunes grow where wind blows sand off the beach and dune
plants trap it. Teixeira et al. ([2023](#ref-Teixeira2023)) study this
with DuBeVeg, a cellular automaton in the family of Werner’s sand-slab
model: the surface is a stack of slabs on each cell, slabs are picked up
and hop downwind until they land, avalanches keep slopes below the angle
of repose, the sea resets the beach below the water level, and two kinds
of plant grow or die back with the sand they gain. In the original model
a slab hops one cell; the authors make it hop a saltation length of 5
cells, with slabs five times thinner so that the potential transport is
the same, and let vegetation trap slabs part-way through a hop. On a
wide beach this changes how the dune grows.

This implementation is an **adaptation**: the authors’ code is not
published, the domain is narrower, and several constants come from the
base model’s documentation.

``` julia
ENV["GKSwstype"] = "100"
using EcologicalCellularAutomata
using CellularAutomataCore
using Plots
using Random
using Statistics
```

## The model

Each aeolian iteration (50 a year), every bare, dry cell outside the
shadow of higher ground gives up a slab with probability `(1 - ρ) 0.5`,
where `ρ` is its plant cover. A slab moves landward and, at the end of
each hop of `L` cells, lands with probability `ρ (1 - 0.1) + 0.1`
(always in the shadow); with `L = 5`, it can also land part-way through
a hop on cells whose plant cover is at least 0.3. Every second iteration
the sea, with run-up from waves, moves the cells it reaches towards the
reference beach, removing or adding each slab between them with the same
probability. Once a year the plants grow or die back with the sand each
cell gained, each species spreading from neighbours of its own, and the
year’s sea thins them.

The beach rises from the shoreline to the foot of a foredune 15 m high
and 100 m wide; the source compares a narrow beach of 100 m with a wide
one of 500 m.

``` julia
plot(; xlabel="cross-shore distance (m)", ylabel="elevation (m)", size=(700, 300),
     legend=:topleft)
for (beach, c) in ((100, 1), (500, 2))
    prof = dune_profile(; beach, h=0.02)
    plot!(prof.elevation .* 0.02; lw=2, color=c, label="$beach m beach")
end
hline!([0]; color=:gray, ls=:dash, label="mean sea level")
```

![](20_dune_formation_files/figure-commonmark/cell-3-output-1.svg)

### Water levels

The source forces the sea with the highest water level in each
spring-neap cycle at Scheveningen, 1970 to 2020, which it plots (its
Fig. 5d) but does not tabulate. The measurements are public, from
Rijkswaterstaat, and `scheveningen_water_levels` holds them: the highest
level in each block of 14.6 days, 25 a year. They agree with the figure,
whose median is about 1.5 m, with a tenth of the cycles above about 2 m
and a highest level of about 3.2 m, and they run higher in winter than
in summer. Wave run-up adds a few tenths of a metre more, and the level
the sea reaches most years keeps the beach below it bare. The levels are
against the Dutch datum, NAP, where the source writes mean sea level; no
correction is made. `problem(:dune_formation; levels=:synthetic)` draws
a series with much the same distribution from the problem’s seed, 1.43 m
plus a Gumbel-distributed surge of scale 0.24 m (`dune_water_levels`).

``` julia
levels = scheveningen_water_levels()
plot((1:length(levels)) ./ 25, levels; lw=0.6, color=:black, legend=false,
     xlabel="year", ylabel="water level (m)", size=(700, 240), left_margin=4Plots.mm,
     bottom_margin=4Plots.mm)
```

![](20_dune_formation_files/figure-commonmark/cell-4-output-1.svg)

``` julia
in_year(months...) = [l for (k, l) in enumerate(levels) if mod1(k, 25) in months]
(median=round(median(levels); digits=2), q90=round(quantile(levels, 0.9); digits=2),
 highest=round(maximum(levels); digits=2),
 winter=round(mean(in_year(23, 24, 25, 1, 2, 3)); digits=2),
 summer=round(mean(in_year(9:15...)); digits=2))
```

    (median = 1.52, q90 = 2.04, highest = 3.12, winter = 1.79, summer = 1.38)

### Profiles from pieces

A profile is built by joining pieces from the sea to the land: a Dean
shoreface, a straight beach, the dune ridge and the plateau behind it. A
piece runs from one elevation to another, and `*` joins two pieces only
where the first ends at the elevation where the second starts. Joining
is associative and an empty piece changes nothing, so profiles compose
like paths between elevations (morphisms in a category whose objects are
elevations), and a scenario is a composite:

``` julia
narrow = shoreface_segment(60) * slope_segment(100, 0, 2.5) *
         ridge_segment(100, 2.5, 15, 5) * flat_segment(20, 5; kind=:hinterland)
(narrow, profile_heights(narrow, 280; h=0.02) == dune_profile(; beach=100).elevation)
```

    (ProfileSegment(shoreface * beach * dune * hinterland: -1.8391426377445272 → 5.0 m over 280 m), true)

The reference surface the sea works towards is the same profile with the
dune and the land behind it levelled at the dune’s foot.
`reference_profile` does this piece by piece, so it respects joins: it
maps composites to composites.

Teixeira et al.’s third scenario is a **prograding** beach, which builds
out 4 m a year, from 100 m to 300 m wide over 50 years. In their Fig. 2b
the profile below a berm at 2 m moves seaward and the berm lengthens.
With pieces, that is one more piece: the beach is cut where it reaches 2
m and a level berm of `4t` m is joined in. Level pieces at one elevation
add, a berm of `a` m joined to one of `b` m being a berm of `a + b` m,
so building out is time acting on the profile by joining berms.
`problem(:dune_formation; scenario=:prograding)` puts the narrow beach
on a domain 200 m longer out to sea, anchored at the land, and replaces
the reference surface with the year’s at the end of every year; the sea,
resetting cells towards it, builds the beach out.

``` julia
plot(; xlabel="cross-shore distance (m)", ylabel="elevation (m)", size=(700, 300),
     legend=:topleft, ylims=(-3, 16))
for y in 0:10:50
    ref = reference_profile(dune_landform(; beach=100, offshore=260, berm=4y))
    plot!(profile_heights(ref, 480); lw=1.5, color=cgrad(:copper)[y / 50], label="year $y")
end
plot!(profile_heights(dune_landform(; beach=100, offshore=260), 480); color=:black, lw=2,
      label="the dune at the start")
```

![](20_dune_formation_files/figure-commonmark/cell-7-output-1.svg)

## The same potential transport

The source’s Eq. 3 gives the potential transport
`Q = h L (pe0 / pd0) n`: with slabs of 0.1 m hopping 1 cell or of 0.02 m
hopping 5, 25 m³ per metre of beach a year. On a bare, dry beach the
slabs crossing a line do carry that:

``` julia
function flat_flux(h, L; rows=40, iterations=20)
    p = merge(defaults(build(:dune_formation)),
              (; h, L, rho_star=Inf, hydro_every=10^6, iterations=10^6))
    total = 0
    for t in 1:iterations
        u = (topo=fill(Int32(round(2.9 / h)), rows, 400), pioneer=zeros(rows, 400),
             stabilizer=zeros(rows, 400), balance=zeros(Int32, rows, 400),
             zref=zeros(Int32, rows, 400), sea=zeros(rows, 400))
        u.topo[:, 200:end] .= round(Int32, 3.1 / h)
        total += dune_step!(u, p, t, Random.Xoshiro(t))
    end
    return round(total * h * 50 / (rows * iterations); digits=1)
end
(hop_1=flat_flux(0.1, 1), hop_5=flat_flux(0.02, 5))
```

    (hop_1 = 24.6, hop_5 = 24.9)

## Fifty years of dune growth

One run of 50 years for each beach and representation, on 20 rows of 1 m
(the source runs 400 rows five times). We record the sand crossing the
dune toe, the first cell above 3 m, each year, and the dune’s volume and
width every ten years.

``` julia
function dune_run(scenario, representation; years=50, rows=20, seed=1)
    prob = problem(:dune_formation; scenario, representation, alongshore=rows, years, seed)
    p, u, h = prob.p, deepcopy(prob.u0), prob.p.h
    flux = zeros(years)
    measures = [dune_measures(u; h)]
    profiles = [vec(mean(u.topo; dims=1)) .* h]
    for t in 1:(50years)
        crossed = dune_step!(u, p, t, Random.Xoshiro(hash((prob.seed, t, :dune))))
        flux[cld(t, 50)] += crossed * h / rows
        if t % 500 == 0
            push!(measures, dune_measures(u; h))
            push!(profiles, vec(mean(u.topo; dims=1)) .* h)
        end
    end
    return (; flux, measures, profiles, final=u.topo .* h)
end
cases = [(s, r) for s in (:narrow, :wide, :prograding) for r in (:default, :saltation)]
runs = Dict(zip(cases, fetch.([Threads.@spawn dune_run(c...) for c in cases])))
nothing
```

The same runs drive the integrator: `dune_step!` with the generator
seeded from the problem’s seed and the iteration is what `solve` runs at
each step.

``` julia
plot([begin
          r = runs[(s, rep)]
          plot(r.profiles[1]; color=:black, lw=1.5, label="start",
               title="$s beach, hop $(rep === :saltation ? 5 : 1)", titlefontsize=10,
               xlabel="cross-shore (m)", ylabel="elevation (m)", legend=:topleft,
               xlims=s === :wide ? (400, 680) : s === :prograding ? (100, 480) : (0, 280),
               ylims=(-1, 22))
          plot!(r.profiles[end]; color=2, lw=2, label="after 50 years")
      end
      for (s, rep) in cases]...; layout=(3, 2), size=(850, 800))
```

![](20_dune_formation_files/figure-commonmark/cell-10-output-1.svg)

``` julia
heatmap(hcat(permutedims(runs[(:wide, :default)].final[:, 400:end]),
             permutedims(runs[(:wide, :saltation)].final[:, 400:end])); c=:terrain,
        clims=(-1, 20), xticks=false, ylabel="cross-shore (m, from 400)", size=(600, 360),
        title="wide beach after 50 years: hop 1 (left), hop 5 (right)", titlefontsize=10)
```

![](20_dune_formation_files/figure-commonmark/cell-11-output-1.svg)

``` julia
source_volume = Dict((:narrow, :default) => 1500, (:narrow, :saltation) => 1611,
                     (:wide, :default) => 1563, (:wide, :saltation) => 2234,
                     (:prograding, :default) => 1540, (:prograding, :saltation) => 2130)
source_width = Dict((:narrow, :default) => 145, (:narrow, :saltation) => 150,
                    (:wide, :default) => 220, (:wide, :saltation) => 230,
                    (:prograding, :default) => 304, (:prograding, :saltation) => 218)
[(beach=s, hop=rep === :saltation ? 5 : 1,
  toe_flux=round(mean(runs[(s, rep)].flux[2:end]); digits=1),
  volume_gain=round(Int,
                    runs[(s, rep)].measures[end].volume - runs[(s, rep)].measures[1].volume),
  source_gain=source_volume[(s, rep)] - 1150,
  width=round(Int, runs[(s, rep)].measures[end].width), source_width=source_width[(s, rep)])
 for (s, rep) in cases]
```

    6-element Vector{@NamedTuple{beach::Symbol, hop::Int64, toe_flux::Float64, volume_gain::Int64, source_gain::Int64, width::Int64, source_width::Int64}}:
     (beach = :narrow, hop = 1, toe_flux = 18.2, volume_gain = 371, source_gain = 350, width = 166, source_width = 145)
     (beach = :narrow, hop = 5, toe_flux = 19.4, volume_gain = 572, source_gain = 461, width = 159, source_width = 150)
     (beach = :wide, hop = 1, toe_flux = 14.0, volume_gain = 392, source_gain = 413, width = 235, source_width = 220)
     (beach = :wide, hop = 5, toe_flux = 22.5, volume_gain = 1029, source_gain = 1084, width = 242, source_width = 230)
     (beach = :prograding, hop = 1, toe_flux = 13.9, volume_gain = 376, source_gain = 390, width = 274, source_width = 304)
     (beach = :prograding, hop = 5, toe_flux = 23.3, volume_gain = 1045, source_gain = 980, width = 216, source_width = 218)

The sand crossing the dune toe, in m³ per metre a year, the gain in dune
volume over 50 years and the dune’s width after them, against the
source’s medians (from an initial volume of about 1150 m³/m; the medians
of the narrow and wide beaches are in its text, those of the prograding
beach and the wide beach’s widths with hops of 1 cell read from its Fig.
4). As in the source, hops of 5 cells carry sand across a wide beach to
the dune close to the potential transport, while single-cell hops build
bedforms on the beach that hold sand back. So the dune gains much more
with hops of 5 cells on the wide beach, here as in the source; on the
narrow beach the fluxes differ less.

``` julia
plot(1:50, [runs[c].flux for c in cases]; lw=1.5,
     label=permutedims(["$(c[1]), hop $(c[2] === :saltation ? 5 : 1)" for c in cases]),
     xlabel="year", ylabel="flux at the dune toe (m³/m/yr)", ylims=(0, 30),
     size=(650, 340), legend=:bottomright)
hline!([25]; color=:gray, ls=:dash, label="potential transport")
```

![](20_dune_formation_files/figure-commonmark/cell-13-output-1.svg)

On the prograding beach, with hops of one cell, the beach builds out
faster than the dune can follow and smaller foredunes form on it,
seaward of the main dune, which is wide, as in the source (274 m in this
run, against its 304 m, its widest). With hops of 5 cells the sand
reaches the old dune, which grows a wide new ridge in front of it (216 m
against 218 m), and the flux at the dune toe rises a little from the
first decade (22.4 to about 23.6 m³/m/yr) as the beach widens, as in the
source’s Fig. 5c.

Those are single runs, with means over rows. Over five seeds, with the
source’s statistic, the median over all profiles (and in brackets the
range of the seeds’ medians):

| beach      | hops | width (m)       | source | gain (m³/m) | source |
|------------|------|-----------------|--------|-------------|--------|
| narrow     | 1    | 166 (164-166)   | 145    | 370         | 350    |
| narrow     | 5    | 163 (159-165)   | 150    | 583         | 461    |
| wide       | 1    | 238.5 (228-245) | 220    | 392         | 413    |
| wide       | 5    | 258 (242-259)   | 230    | 1036        | 1084   |
| prograding | 1    | 284 (242-288)   | 304    | 380         | 390    |
| prograding | 5    | 213 (208-216)   | 218    | 1051        | 980    |

The gains are within about 7% of the source’s, except on the narrow
beach with hops of 5 cells, a quarter too large. The widths are off by
17.6 m on average: on the two beaches with a fixed shoreline the dune is
too wide, by 13 to 21 m on the narrow and 18.5 and 28 m on the wide. An
earlier version of this model reset a whole column at once and had a
simpler plant rule. It gave 231 m on the wide beach with hops of 5
cells, but those rules are not the source’s, so that agreement says
little.

The narrow beach’s dune advances at about a metre a year to the end of
the run (127, 141, 151, 159 and 166 m by decade, with hops of 1 cell;
seed 1). Its front is a nearly bare face at the angle of repose, rising
0.5 m per cell from about 1.3 m, with plant cover below a quarter up to
about 7 m, its foot cut by the sea. The sea’s strength does not set the
rate: with no resistance from plants (`Rv = 0`) the width is still 166
m, with `G` of 0.7 or 0.9 it is 166 or 165 m, and with no dissipation
(`Wdiss = 0`) 161 m. Each of these changes it (one run each, against 145
m and +350 m³/m in the source):

| change                                   | width (m) | gain (m³/m) |
|------------------------------------------|-----------|-------------|
| the run above                            | 166       | 369         |
| water levels 0.5 m higher                | 130       | 320         |
| water levels 0.5 m lower                 | 181       | 434         |
| waves at their least strength everywhere | 180.5     | 411         |
| half the erosion probability             | 144       | 234         |
| twice the deposition probability         | 141       | 242         |
| no plants establishing anew              | 152       | 461         |
| no new plants at all                     | 144       | 368         |

The last row is close to the source’s pair, though the source’s plants
do colonise, so it is no explanation. The source does not publish the
probability that plants establish, 0.05 here from DuBeVeg’s code. In a
sweep from 0 to 0.1 over the six cases, run on this model before its
avalanches took their present form, no value fitted both the widths and
the gains: lower values narrowed the dunes on the fixed shorelines and
made the narrow beach gain far too much. Nor does it publish how the
dissipation `D` accumulates. The waves follow DuBeVeg’s code: their
initial strength is the total water level (in slabs) times an energy
factor of 1, which keeps them strong far up the beach and, near the
shoreline, overrules the protection of plants; the source’s Eq. 4 leaves
the factor out. Vignette 28 shows that a wind steered by the topography
also narrows the dunes, though the source’s wind is not steered.

## Deviations from the source

``` julia
foreach(d -> println("- ", d), model_info(:dune_formation).deviations)
```

    - The sea is forced by the water levels measured at Scheveningen (1970-2020), the highest in each spring-neap cycle, as in the source, taken from Rijkswaterstaat rather than from the source, which only plots them. They are against NAP, where the source writes mean sea level, with no correction. A synthetic series with much the same distribution (1.43 m plus a Gumbel-distributed surge of scale 0.24 m) remains available with levels = :synthetic.
    - The domain is 20 to 100 rows of 1 m alongshore (400 in the source) and only the cross-shore profile from 60 m offshore to 20 m behind the dune (1250 cells in the source). The initial profile is a Dean-like shoreface, a linear beach and a smooth dune ridge to 15 m read from the source's figure; the beach rises to 2.5 m at the dune foot, as the figure shows, where the text says 3 m (with 3 m, small bumps on the upper beach count as the dune).
    - Slab transport is synchronous, as in DuBeVeg's code: all pick-ups are drawn from the state at the start of an iteration, with the shadow and the probabilities fixed for its duration.
    - Part-way through a hop (Eq. 5), a slab lands only where vegetation is at least rho_star, and the shadow does not force it to; after a failed landing at a hop's end the next hop begins. Slabs reaching the land side of the domain are lost.
    - The shadow test uses the continuous threshold i tan(15°); on whole slabs it is the same test as DuBeVeg's floored one, except that DuBeVeg's code looks upwind no further than the highest elevation divided by tan(15°), which matters only for cells below zero.
    - Constants not in the source come from DuBeVeg: the growth functions of the pioneer (burial -1.4, 0.2, 0.6, 2.0, 2.2 m/yr, peak 0.2) and stabilizer (-1.4, -0.65, 0, 0.2, 2.8, peak 0.05), lateral expansion 0.2 and establishment 0.05 a year (each species spreading only from neighbours of its own, its yearly gain at most the cover left free, establishment only from 0 to 20 m high, where the source gives every bare cell the same chance, and for new stabilizers only where the sand balance is zero, where DuBeVeg's code asks that no sand moved at all that year), the wave run-up relation, the dissipation by the inverse water depth summed from the sea, the least wave probability 0.1, and the waves' initial strength, the total water level in slabs times Fenergy = 1 (wave_energy), which multiplies WE - Wdiss D in DuBeVeg's code (the source's Eq. 4 leaves it out); near the shoreline it overrules the protection of vegetation, as in DuBeVeg. Plant cover is multiplied each year by one minus the sum of the sea's probabilities over the year (at most 1), as in DuBeVeg's code, and otherwise stays through each reset. The initial split of species is not stated: pioneers seaward of the crest, stabilizers behind.
    - Avalanching follows DuBeVeg's code: each cell passes a slab to its steepest lower neighbour (by height difference) only if that drop exceeds the angle of repose (times the square root of 2 along a diagonal). The source allows only slopes equal to or gentler than the angle of repose; the code leaves an orthogonal slope steeper than that where a lower diagonal neighbour is within its own limit.
    - The dune is wider than in the source after 50 years, medians of five seeds, on the narrow beach (166 and 163 m with hops of 1 and 5 cells, against 145 and 150) and on the wide beach (238.5 and 258 m against 220 and 230), and the narrow beach with hops of 5 cells gains too much (583 m³/m against 461). The cause is not known. On the narrow beach the front is a nearly bare face at the angle of repose, its foot cut by the sea, advancing about a metre a year, and the sea's strength does not set the rate. With no new plants at all the narrow beach with hops of 1 cell gives 144 m and +368 m³/m (the source 145 and +350), but the source's plants do colonise; it does not publish the probability that they establish, nor how the dissipation D accumulates.
    - The prograding beach builds out through the sea alone: the reference surface is replaced at the end of each year by the profile with a berm 4 m longer, the profile below the berm moved seaward within a domain extending 200 m further out (the source's profile reaches a depth of closure of 8 m, beyond the domain here), and the sea's resets move the beach towards it; the shoreline then follows the reference about a year behind.
    - The sensitivity to rho_star is not included; one run per case on 20 rows replaces the source's five runs on 400.
    - Beyond the source: three more landforms built from the same profile pieces (an embryo dune and slack seaward of the foredune on the wide beach, a seawall and promenade at the back of the narrow beach, and a breach through the foredune in a fifth of the rows), and a non-erodible layer (hard) for the seawall; the source runs none of them, so they test the model's behaviour, not its fidelity. A breach fills from its flanks rather than growing into a blowout.
    - Beyond the source: an optional wind steered by the topography (speedup > 0; off by default), from the linear theory of turbulent flow over low hills (Kroy et al. 2002 along the wind, Schwammle and Herrmann 2005 across it), computed on the height above the reference surface with the shadow envelope in the lee. The erosion probability scales with the shear along the wind and flying slabs turn alongshore with the shear across it. With A = 4 it narrows the dunes, a little past the source's width on both beaches with hops of 5 cells (144 and 219 m against 150 and 230), but a breach still fills.

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Teixeira2023" class="csl-entry">

Teixeira, Manuel, Erik M. Horstman, and Kathelijne M. Wijnberg. 2023.
“Conceptualizing Aeolian Sediment Transport in a Cellular Automata Model
to Simulate the Bio-Geomorphological Evolution of Beach–Dune Systems.”
*Journal of Marine Science and Engineering* 11 (7): 1278.
<https://doi.org/10.3390/jmse11071278>.

</div>

</div>
