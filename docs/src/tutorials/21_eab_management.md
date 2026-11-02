# Managing emerald ash borer outliers
Simon Frost

- [Overview](#overview)
- [The treatments](#the-treatments)
- [Treating within 300 m](#treating-within-300-m)
- [The radius of treatment and the year it
  starts](#the-radius-of-treatment-and-the-year-it-starts)
- [How many trees can be treated](#how-many-trees-can-be-treated)
- [The limits on girdling](#the-limits-on-girdling)
- [How stress is counted](#how-stress-is-counted)
- [Deviations from the source](#deviations-from-the-source)

## Overview

The emerald ash borer spreads over long distances when people move
infested wood, founding outlier colonies far ahead of the main front.
Once one is found, what should be done there? Mercader et al.
([2011b](#ref-Mercader2011b)) took their model of local spread (Mercader
et al. ([2011a](#ref-Mercader2011)); see the emerald ash borer vignette)
and added three treatments in a circle around the origin of an outlier:

- **ash reduction**: harvesting the large trees, removing their phloem
  and the larvae in it, once;
- **girdling**: girdling one large tree a cell each year, which draws
  egg-laying females to it, and destroying it with its larvae;
- **insecticide**: protecting the large trees with a systemic
  insecticide, which kills the larvae laid in them and keeps their
  phloem from being eaten.

Large trees (27.9 cm in diameter or more) hold half the phloem by
default. The authors also asked how the results change if females are
drawn to stressed trees.

This implementation is an **adaptation**: the authors’ code was not
available, so the treatments are written from the text, with the gaps
filled as the deviations below list.

``` julia
ENV["GKSwstype"] = "100"
using EcologicalCellularAutomata
using CellularAutomataCore
using Plots
```

## The treatments

The landscape is 101 × 101 cells of 50 m, each with 200 m² of ash
phloem, with 100 adults released at the centre. Treated cells lie within
a radius of the origin. Without treatment the model is the base model
with avoidance: its radius of infestation after 15 years is 1350 m.

``` julia
radii(sol) = [eab_spread(u).radius for u in sol.u]
simulate(; kw...) = solve(problem(:eab_management; kw...); saveat=1)
untreated = simulate()
radii(untreated)'
```

    1×16 adjoint(::Vector{Float64}) with eltype Float64:
     0.0  150.0  150.0  250.0  300.0  400.0  …  1050.0  1150.0  1250.0  1350.0

A girdled tree of 10.32 m² takes a share of the larvae laid in its cell
given by the ratio of the density of larvae on girdled to other trees.
The source fits that ratio as `20.468 x^-0.559` at a density `x` per m²
on other trees, from 16 outlier sites in Michigan; here it is capped at
20.468, its value at 1 larva per m². The ratio falls as density rises:
in a cell of 200 m², a girdled tree takes about half the eggs at the low
densities ahead of a front, and 8% where the phloem is full.

``` julia
x = 0.5:0.5:88
y = min.(20.468 .* x .^ -0.559, 20.468)
share = 10.32 .* y ./ (10.32 .* y .+ (200 - 10.32))
plot(plot(x, y; ylabel="ratio, girdled to other trees"),
     plot(x, 100 .* share; ylabel="share in the girdled tree (%)", ylims=(0, 60));
     xlabel="larvae per m² of phloem", lw=2, legend=false, layout=(1, 2), size=(850, 320),
     left_margin=5Plots.mm, bottom_margin=5Plots.mm)
```

![](21_eab_management_files/figure-commonmark/cell-4-output-1.svg)

Under the **stress bias**, dispersing females favour cells with a
girdled tree (seven times as attractive) or stressed trees, whose
attraction rises linearly from 1 to 7 as the cell loses up to a third of
its phloem.

## Treating within 300 m

The source’s main comparison treats every cell within 300 m of the
origin from the start and follows the infestation for 15 years.

``` julia
treatments = (:none, :girdling, :insecticide, :ash_reduction)
labels = ["standard" "girdled trees" "insecticide" "ash reduction"]
runs = Dict((tr, sb) => simulate(; treatment=tr, stress_bias=sb)
            for tr in treatments for sb in (false, true))
plot([plot(0:15, [radii(runs[(tr, sb)]) for tr in treatments]; label=labels,
           title=sb ? "stress bias" : "no stress bias", titlefontsize=10, xlabel="year",
           ylabel="radius of infestation (m)", lw=2, legend=:topleft, ylims=(0, 1500),
           color=[:black :gray :black :gray], ls=[:solid :solid :dash :dash])
      for sb in (false, true)]...; layout=(1, 2), size=(850, 380), left_margin=5Plots.mm,
     bottom_margin=5Plots.mm)
```

![](21_eab_management_files/figure-commonmark/cell-5-output-1.svg)

``` julia
source = Dict((:none, false) => (1350, 5.45), (:girdling, false) => (1150, 4.45),
              (:insecticide, false) => (950, 3.3), (:ash_reduction, false) => (1350, 5.3),
              (:none, true) => (1100, 2.5), (:girdling, true) => (950, 2.1),
              (:insecticide, true) => (850, 1.8), (:ash_reduction, true) => (1150, 2.8))
[(; treatment=tr, stress_bias=sb, radius=Int(radii(runs[(tr, sb)])[end]),
  source_radius=source[(tr, sb)][1],
  beyond_millions=round(eab_beyond(runs[(tr, sb)].u[end], 300) / 1e6; digits=2),
  source_beyond=source[(tr, sb)][2])
 for sb in (false, true) for tr in treatments]
```

    8-element Vector{@NamedTuple{treatment::Symbol, stress_bias::Bool, radius::Int64, source_radius::Int64, beyond_millions::Float64, source_beyond::Float64}}:
     (treatment = :none, stress_bias = 0, radius = 1350, source_radius = 1350, beyond_millions = 5.46, source_beyond = 5.45)
     (treatment = :girdling, stress_bias = 0, radius = 1100, source_radius = 1150, beyond_millions = 4.25, source_beyond = 4.45)
     (treatment = :insecticide, stress_bias = 0, radius = 1000, source_radius = 950, beyond_millions = 3.71, source_beyond = 3.3)
     (treatment = :ash_reduction, stress_bias = 0, radius = 1400, source_radius = 1350, beyond_millions = 5.67, source_beyond = 5.3)
     (treatment = :none, stress_bias = 1, radius = 1050, source_radius = 1100, beyond_millions = 2.73, source_beyond = 2.5)
     (treatment = :girdling, stress_bias = 1, radius = 900, source_radius = 950, beyond_millions = 2.12, source_beyond = 2.1)
     (treatment = :insecticide, stress_bias = 1, radius = 950, source_radius = 850, beyond_millions = 2.27, source_beyond = 1.8)
     (treatment = :ash_reduction, stress_bias = 1, radius = 1150, source_radius = 1150, beyond_millions = 2.92, source_beyond = 2.8)

The radius after 15 years and the beetles beyond the treated area
(millions), against the source’s Fig. 4 read by eye. As in the source:

- insecticide slows the spread most, then girdling;
- ash reduction does not reduce the beetles beyond the treated area;
- with the stress bias, ash reduction spreads faster than no treatment.

Insecticide is a little weaker here (1000 m against 950 without the
bias) and girdling a little stronger.

## The radius of treatment and the year it starts

The source then treated every cell out to radii of 0 to 1700 m, starting
in the first year or up to four years later, when an outlier might be
found.

``` julia
rs = 0:100:1700
starts = (0, 2, 4)
cases = [(tr, sb, st, r) for tr in (:girdling, :insecticide, :ash_reduction)
         for sb in (false, true) for st in starts for r in rs]
final = Vector{Float64}(undef, length(cases))
Threads.@threads for i in eachindex(cases)
    tr, sb, st, r = cases[i]
    prob = problem(:eab_management; treatment=tr, stress_bias=sb, start=st, radius=r)
    final[i] = eab_spread(solve(prob).u[end]).radius
end
sweep = Dict(cases .=> final)
plot([plot(rs, [[sweep[(tr, sb, st, r)] for r in rs] for st in starts];
           label=permutedims(["year $st" for st in starts]), lw=1.5,
           title="$(replace(String(tr), '_' => ' ')), $(sb ? "" : "no ")stress bias",
           titlefontsize=9, ylims=(600, 1650), xlabel="radius of treatment (m)",
           ylabel="radius after 15 years (m)", legend=:bottomleft)
      for tr in (:girdling, :insecticide, :ash_reduction) for sb in (false, true)]...;
     layout=(3, 2), size=(800, 860))
```

![](21_eab_management_files/figure-commonmark/cell-7-output-1.svg)

The results follow the source’s Fig. 5:

- **Insecticide:** the radius falls to 800 m as the treated radius grows
  (750 m with the stress bias; the source has 800 and 700 m). Later
  starts weaken it.
- **Girdling:** the radius falls to 1100 m once the treated circle
  passes 300–800 m, depending on the start, and stays there.
- **Under the stress bias:** girdling helps most at small radii and less
  at large ones, where the girdled cells draw beetles outwards, as the
  source found. Ash reduction widens spread as the treated circle grows,
  to 1300–1450 m.

## How many trees can be treated

The source varied the share of the phloem in large trees, and so the
share that can be treated, from 0 to 90%. Here the whole lattice is
treated.

``` julia
shares = 0:0.1:0.9
share_cases = [(tr, sb, s) for tr in (:girdling, :insecticide, :ash_reduction)
               for sb in (false, true) for s in shares]
share_final = Vector{Float64}(undef, length(share_cases))
Threads.@threads for i in eachindex(share_cases)
    tr, sb, s = share_cases[i]
    prob = problem(:eab_management; treatment=tr, stress_bias=sb, share=s, radius=1700)
    share_final[i] = eab_spread(solve(prob).u[end]).radius
end
by_share = Dict(share_cases .=> share_final)
plot([plot(100 .* shares,
           [[by_share[(tr, sb, s)] for s in shares]
            for tr in (:girdling, :insecticide, :ash_reduction)];
           label=["girdled trees" "insecticide" "ash reduction"], lw=2,
           title=sb ? "stress bias" : "no stress bias", titlefontsize=10,
           xlabel="phloem available for treatment (%)",
           ylabel="radius after 15 years (m)", ylims=(0, 1700), legend=:bottomleft)
      for sb in (false, true)]...; layout=(1, 2), size=(850, 380), left_margin=5Plots.mm,
     bottom_margin=5Plots.mm)
```

![](21_eab_management_files/figure-commonmark/cell-8-output-1.svg)

Insecticide’s effect grows steadily with the share treated, to 150 m at
90%, as in the source’s Fig. 3 (about 200 m). Girdling helps less,
because each cell girdles only one tree a year and the large trees run
out. Ash removal never slows the spread. Removing 75–90% of the phloem
widens it here, which the source does not show.

## The limits on girdling

Girdling is limited in two ways: only large trees can be girdled, and
girdled trees draw fewer eggs as densities rise. The source’s Fig. 6
lifts each limit in turn, girdling throughout the lattice.

``` julia
scenarios = [("standard", (; treatment=:none)),
             ("all constraints", (; treatment=:girdling)),
             ("all phloem available", (; treatment=:girdling, share=1.0)),
             ("no reduction in attraction",
              (; treatment=:girdling, pmap=(; attraction_decline=false))),
             ("no constraints",
              (; treatment=:girdling, share=1.0, pmap=(; attraction_decline=false)))]
beetles = [[eab_spread(u).beetles / 1e6 for u in simulate(; radius=1700, kw...).u]
           for (_, kw) in scenarios]
plot(0:15, beetles; label=permutedims(first.(scenarios)), lw=2, xlabel="year",
     ylabel="beetles (millions)", legend=:topleft, size=(650, 360))
```

![](21_eab_management_files/figure-commonmark/cell-9-output-1.svg)

Each limit lifted cuts the population further. With neither, the
infestation stays small, as in the source.

## How stress is counted

The source counts a cell as stressed once a third of its “initial
phloem” is gone. Here that means phloem lost by any means. If only
phloem eaten by larvae counted, removing ash would leave the treated
cells unstressed. Under the stress bias it would then barely change the
spread, unlike the source’s Figs. 3 and 5.

``` julia
lost = [eab_spread(solve(problem(:eab_management; treatment=:ash_reduction,
                                 stress_bias=true, radius=r)).u[end]).radius
        for r in (0, 300, 800, 1700)]
eaten = [eab_spread(solve(problem(:eab_management; treatment=:ash_reduction,
                                  stress_bias=true, radius=r,
                                  pmap=(; stress_from=:eaten))).u[end]).radius
         for r in (0, 300, 800, 1700)]
(; radius_of_treatment=[0, 300, 800, 1700], any_loss=lost, eaten_only=eaten)
```

    (radius_of_treatment = [0, 300, 800, 1700], any_loss = [1050.0, 1150.0, 1250.0, 1450.0], eaten_only = [1050.0, 1050.0, 1050.0, 1050.0])

## Deviations from the source

``` julia
foreach(d -> println("- ", d), model_info(:eab_management).deviations)
```

    - The source's code (its Supplemental Appendix S1) was not available; the treatments are written from the text, on the emerald_ash_borer model with avoidance, on 101 x 101 cells (the source's 100 x 100 has no centre cell).
    - The ratio of larval density on girdled to other trees, 20.468 x^-0.559, is taken at the cell's density x = N / A before girdling, and capped at 20.468 (its value at 1 larva per m², near the largest ratio observed); uncapped, it removes nearly all larvae at the front and girdling is far stronger than in the source. Without its decline the ratio is 20.468 throughout.
    - Under the stress bias, a cell's stress is the share of its initial phloem lost by any means (eaten, removed or girdled), so ash reduction stresses the treated cells; counting only eaten phloem (stress_from = :eaten) leaves ash reduction under the stress bias with no effect on the radius, unlike the source.
    - Attraction multiplies each receiving cell's weight in the dispersal kernel before the source-by-source normalisation, as the source says dispersal is proportional. With no stressed cells this is the plain kernel, so the radius after the first year is 150 m with or without the bias; the source's stress-bias runs show 200 m.
    - When a treatment starts in a later year, ash reduction removes the treated share of the larvae in the cell and insecticide kills the larvae in the protected phloem; girdling and insecticide treat min(share x initial phloem, current phloem).
    - With 75-90% of phloem removed, ash reduction widens spread here (1500-1600 m) where the source shows no change; the stress-bias radius at 0 m of treatment is 1000-1050 m against the source's 1050-1100.

```@raw html
<div id="refs" class="references csl-bib-body hanging-indent">
```

```@raw html
<div id="ref-Mercader2011" class="csl-entry">
```

Mercader, Rodrigo J., Nathan W. Siegert, Andrew M. Liebhold, and Deborah
G. McCullough. 2011a. “Influence of Foraging Behavior and Host Spatial
Distribution on the Localized Spread of the Emerald Ash Borer, *Agrilus
Planipennis*.” *Population Ecology* 53 (2): 271–85.
<https://doi.org/10.1007/s10144-010-0233-6>.

```@raw html
</div>
```

```@raw html
<div id="ref-Mercader2011b" class="csl-entry">
```

Mercader, Rodrigo J., Nathan W. Siegert, Andrew M. Liebhold, and Deborah
G. McCullough. 2011b. “Simulating the Effectiveness of Three Potential
Management Options to Slow the Spread of Emerald Ash Borer (*Agrilus
Planipennis*) Populations in Localized Outlier Sites.” *Canadian Journal
of Forest Research* 41 (2): 254–64. <https://doi.org/10.1139/X10-201>.

```@raw html
</div>
```

```@raw html
</div>
```
