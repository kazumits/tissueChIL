---
title: "Modeling traveling ratio"
output: 
  html_document: 
    toc: yes
    keep_md: yes
---

The read counts of RNAP2 ChIL-seq at TSS- and TES-region were fitted to the following Poisson regression model. The model can evaluate mean, variance and thus confidence intervals of traveling ratio by utilizing all (five) replicates that have different total sequenced reads.

For each gene, assume read count $y_{ij}$ of $i$-th replicate at side $j$ (TSS, TES) follows a Poisson distribution, and the mean $\lambda_{ij}$ satisfies the relation:

$$
\begin{aligned}
\lambda_{ij}/M_i &= \exp(\beta_0+\beta_1 s_{ij})\\
\log\lambda_{ij} &= \log M_i + \beta_0 + \beta_1 s_{ij}
\end{aligned}
$$
, where the offset term $M_i$ is the total reads (in million) of replicate $i$ and $s_{ij}$ is the indicator variable that the read count $y_{ij}$ is either from TSS ($s_{ij}=0$) or TES ($s_{ij}=1$). Since the offsetting $\lambda_{ij}/M_i$ is equivalent to CPM normalization of the mean count, $e^{\beta_0}$ and $e^{\beta_1}$ can be reffered to as the expected CPM at TSS and the travaling ratio (TES/TSS) of the gene, respectively.

*Below is an example of using the Poisson regression model to decompose transcriptional activation and population-size effect in bulk tissue epigenome.*

### Load packages


```r
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
library(purrr)
library(tidyr)
```

### Load TSS- and TES-counts of the example transcripts

The example data of time-course RNAP2 ChIL-seq data of injured muscle.


```r
tbl <- readr::read_tsv("data/sampleCountsTSSTES.txt.gz",col_types="cccccciii")
sample_n(tbl,10)
```

```
## # A tibble: 10 x 9
##    ens_tran       ens_gene       ext_gene cell    day   rep     tss   tes  total
##    <chr>          <chr>          <chr>    <chr>   <chr> <chr> <int> <int>  <int>
##  1 ENSMUST000001… ENSMUSG000000… Tnnt3    SKMusc… Day0  4         7    20 3.69e6
##  2 ENSMUST000001… ENSMUSG000000… Cnp      Immune  Day14 5        16    11 5.19e6
##  3 ENSMUST000000… ENSMUSG000000… Pf4      Immune  Day5  2         3     4 3.07e6
##  4 ENSMUST000001… ENSMUSG000000… Tpm1     Immune  Day3  1        13    14 6.28e6
##  5 ENSMUST000001… ENSMUSG000000… Vim      SKMusc… Day5  5        53    19 2.93e6
##  6 ENSMUST000000… ENSMUSG000000… Cdk1     SKMusc… Day0  5         5     6 4.16e6
##  7 ENSMUST000002… ENSMUSG000000… Egfl7    Immune  Day3  1         8    12 6.28e6
##  8 ENSMUST000002… ENSMUSG000000… Myh11    Immune  Day14 2         1     3 5.07e6
##  9 ENSMUST000001… ENSMUSG000000… Spp1     SKMusc… Day5  2         6     9 3.07e6
## 10 ENSMUST000001… ENSMUSG000000… Myl9     Immune  Day3  1         5     7 4.11e6
```

The table is already annotated as above.

### Fitting to GLM and formatting tables

The offset term is equipped by `offset=log(total/1e6)`.


```r
res <- tbl %>%
  pivot_longer(c(tss,tes),names_to="site",values_to="y") %>%
  mutate(site=factor(site,c("tss","tes"))) %>% 
  arrange(ens_tran) %>% group_by(ens_tran,cell,day) %>% nest() %>%
  mutate(model=map(data,~
    glm(y~site,family=poisson,offset=log(total/1e6),data=.x))
  ) %>% ungroup
```

Formatting tables and the base change.


```r
tbl.pois <- res %>%
  mutate(brm=map(model,broom::tidy,conf.int=TRUE)) %>%
  select(-data,-model) %>% unnest(brm) %>%
  mutate_at(
    vars(estimate,std.error,conf.low,conf.high),
    ~./log(2)
  ) %>% mutate(
    term=sub("\\(Intercept\\)","log2TSS",term),
    term=sub("sitetes","log2TR",term),
    day=factor(day,paste0("Day",c(0,3,5,7,14)))
  )
```


```r
sample_n(tbl.pois,10)
```

```
## # A tibble: 10 x 10
##    ens_tran cell  day   term  estimate std.error statistic  p.value conf.low
##    <chr>    <chr> <fct> <chr>    <dbl>     <dbl>     <dbl>    <dbl>    <dbl>
##  1 ENSMUST… SKMu… Day7  log2…    1.53      0.146    10.5   7.48e-26   1.24  
##  2 ENSMUST… Immu… Day0  log2…    0.222     0.803     0.277 7.82e- 1  -1.37  
##  3 ENSMUST… Immu… Day0  log2…    0.627     0.268     2.34  1.92e- 2   0.0685
##  4 ENSMUST… SKMu… Day0  log2…    2.48      0.141    17.6   1.19e-69   2.20  
##  5 ENSMUST… SKMu… Day0  log2…   -1.28      0.131    -9.77  1.50e-22  -1.55  
##  6 ENSMUST… SKMu… Day5  log2…    0.160     0.213     0.751 4.53e- 1  -0.278 
##  7 ENSMUST… SKMu… Day14 log2…    1.63      0.168     9.74  2.13e-22   1.29  
##  8 ENSMUST… Immu… Day0  log2…   -1.28      0.131    -9.77  1.50e-22  -1.55  
##  9 ENSMUST… Immu… Day5  log2…    0.415     0.275     1.51  1.32e- 1  -0.122 
## 10 ENSMUST… SKMu… Day7  log2…    1.89      0.129    14.7   3.45e-49   1.64  
## # … with 1 more variable: conf.high <dbl>
```


### Calculation of the contrasts: DayX/Day0


```r
tbl.day0 <- tbl.pois %>% filter(day=="Day0") %>%
  select(ens_tran,term,estimate0=estimate,std.error0=std.error)

tbl.conta <- tbl.pois %>% 
  filter(day!="Day0") %>% inner_join(tbl.day0) %>% 
  transmute(
    ens_tran,
    cell,
    day,
    term,
    log2FC = estimate-estimate0,
    se = sqrt(std.error^2+std.error0^2),
    stat = log2FC/se,
    p.value = 2*pnorm(abs(stat),lower.tail=FALSE),
  ) %>% group_by(term) %>%
  mutate(p.adj = p.adjust(p.value,method="BH")) %>%
  ungroup
```

```
## Joining, by = c("ens_tran", "term")
```


```r
sample_n(tbl.conta,10)
```

```
## # A tibble: 10 x 9
##    ens_tran        cell    day   term    log2FC    se     stat  p.value    p.adj
##    <chr>           <chr>   <fct> <chr>    <dbl> <dbl>    <dbl>    <dbl>    <dbl>
##  1 ENSMUST0000022… SKMusc… Day14 log2TR  0.255  0.417   0.611  5.41e- 1 7.02e- 1
##  2 ENSMUST0000016… Immune  Day3  log2TR  0.450  0.610   0.738  4.60e- 1 6.40e- 1
##  3 ENSMUST0000013… Immune  Day5  log2TR  2.13   0.585   3.64   2.77e- 4 2.43e- 3
##  4 ENSMUST0000010… Immune  Day3  log2TR -0.233  0.250  -0.932  3.51e- 1 5.49e- 1
##  5 ENSMUST0000008… Immune  Day5  log2TR -0.0286 0.420  -0.0681 9.46e- 1 9.71e- 1
##  6 ENSMUST0000014… Immune  Day14 log2T… -0.903  0.454  -1.99   4.66e- 2 8.32e- 2
##  7 ENSMUST0000010… Immune  Day3  log2T… -1.86   0.156 -11.9    7.50e-33 3.61e-31
##  8 ENSMUST0000011… SKMusc… Day3  log2TR  1.93   0.615   3.13   1.75e- 3 1.19e- 2
##  9 ENSMUST0000004… Immune  Day7  log2T… -0.723  0.353  -2.05   4.05e- 2 7.37e- 2
## 10 ENSMUST0000010… SKMusc… Day5  log2T…  0.882  0.168   5.25   1.51e- 7 8.58e- 7
```

