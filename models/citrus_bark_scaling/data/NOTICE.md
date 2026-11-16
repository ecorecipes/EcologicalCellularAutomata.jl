# Bahia bark scaling of citrus, Laranjeira et al. (2006)

`trees.csv` and `inspections.csv` hold the survey of 240 grapefruit trees in an experimental grove in
north-eastern Brazil, by Laranjeira, Barbosa, Santos-Filho, Gonçalves and Nickel (2006), Annals of Applied
Biology 148, 187-195 (doi 10.1111/j.1744-7348.2006.00057.x). They were taken from Text S1
(`pcbi.1003753.s002.xls`) of Cunniffe, Laranjeira, Neri, DeSimone and Gilligan (2014), PLoS Computational
Biology 10(8), e1003753 (doi 10.1371/journal.pcbi.1003753), published under the Creative Commons Attribution
4.0 licence (CC BY 4.0).

- `trees.csv`: one row per tree, `row` (1 to 16, row 1 nearest the diseased adult trees), `col` (1 to 15)
  and `first_symptomatic_month`, the month since planting of the inspection at which the tree was first seen
  with symptoms, empty for the three trees never seen with symptoms.
- `inspections.csv`: one row per inspection, its `date`, `month` since planting, and the numbers of trees
  newly and cumulatively seen with symptoms.

Changes made: the spreadsheet's grid labels each tree with an inspection month in reverse rank (its largest
label, 55, marks the first inspection and its 0s the last); the counts of the labels match the newly
symptomatic counts of the spreadsheet's time series exactly only when read that way, and the reading agrees
with Figure 2 of Cunniffe et al. (2014), where the rows nearest the source show symptoms first. Months since
first detection were converted to months since planting by adding 27, since Cunniffe et al. (2014) give the
last survey as 82 months after planting and it was 55 months after the first detection. The data are
otherwise unchanged.

These files are distributed under CC BY 4.0, not under this package's licence.
