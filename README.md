# Partitioning density and size components of biodiversity effects
The R function `densize` (DENSIty and SIZE) partitions the net, complementarity, and selection effects of biodiversity into additive components that reflect diversity-induced changes in plant density and size.

## Installation
```{r}
#library(remotes)
remotes::install_github("communityecologist/densize")
library(densize)
```

## Usage
```{r}
densize(init.dens, final.dens, final.yield, germ.dens)
```
- `init.dens` : A matrix or data frame consisting of initial plant density (the number of plants or seeds sown per area).
- `final.dens` : A matrix or data frame consisting of final plant density (the number of plants survived per area).
- `final.yield` : A matrix or data frame consisting of final yield (yield per area).
- `germ.dens` : Optional. A matrix or data frame consisting of germination density (the number of seedlings emerged per area).

Run `?densize` for detail.

## Datasets
```{r}
BioDivExpt.1
BioDivExpt.2
```
- `BioDivExpt.1` : A dataset of a diversity experiment using annual plants. Five species were sown solely or in mixtures of two, three, or four species in planting pots at a total of 24 seeds.
- `BioDivExpt.2` : A dataset of a diversity experiment using annual plants. Five species were sown solely or in mixtures of two, three, or four species in planting pots at a total of 72 seeds.

Run `?BioDivExpt.1` and `?BioDivExpt.2` for detail.

## Examples
```{r}
dat  <- BioDivExpt.1
res1 <- densize(dat$InitDens, dat$FinalDens, dat$FinalYield)
res2 <- densize(dat$InitDens, dat$FinalDens, dat$FinalYield, dat$GermDens)
```

## Citation
[Tatsumi S & Loreau M (in press) Partitioning the biodiversity effects on productivity into density and size components. *Ecology Letters* doi:10.1111/ele.14300.](https://onlinelibrary.wiley.com/doi/full/10.1111/ele.14300)
