## Setup for this tutorial

1. Make sure you have an up-to-date version of [R](https://cloud.r-project.org/). 
   The tutorial should work with R >= 3.3, but it's always a good idea
   to stay up-to-date with R and Bioconductor, as bugs are fixed and
   code is made more efficient.
2. If you haven't already,
   [install Bioconductor](https://bioconductor.org/install). The basic
   installation can be completed with:

```{r}
source("http://bioconductor.org/biocLite.R")
biocLite()
```

3. Install the following packages (if you have some already installed
   you can skip these):

```{r}
pkgs <- c("here","readr","GenomicFeatures","tximport","DESeq2","org.Hs.eg.db","goseq")
find.package(pkgs) # this will give an ERROR if any not found
biocLite(pkgs)
```

If you have any questions about setting up for the tutorial you can
email me at `michaelisaiahlove` at `gmail.com`.

## RNA-seq of colon tissue from patients with Crohn's disease and controls

## References:

*Inflamm Bowel Dis.* 2015 Sep; 21(9): 2178–2187.
**MicroRNAs Classify Different Disease Behavior Phenotypes of Crohn's
Disease and May Have Prognostic Utility.** 
Bailey C. E. Peck, BA, Matthew Weiser, BS, Saangyoung E. Lee, BS,
Gregory R. Gipson, BS, Vishal B. Iyer, Ryan B. Sartor, MD, Hans
H. Herfarth, MD, PhD, Millie D. Long, MD, MPH, Jonathan J. Hansen,
MD, PhD, Kim L. Isaacs, MD, PhD, Dimitri G. Trembath, MD, Reza
Rahbar, MD, Timothy S. Sadiq, MD, Terrence S. Furey,
PhD, Praveen Sethupathy, PhD, and Shehzad Z. Sheikh, MD, PhD

<https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4603665/>

*Gut.* 2016 Oct 14. pii: gutjnl-2016-312518. doi: 10.1136/gutjnl-2016-312518.
**Molecular classification of Crohn's disease reveals two clinically relevant subtypes.**
Weiser M, Simon JM, Kochar B, Tovar A, Israel JW, Robinson
A, Gipson GR, Schaner MS, Herfarth HH, Sartor RB, McGovern DP,
Rahbar R, Sadiq TS, Koruda MJ, Furey TS, Sheikh SZ

<https://www.ncbi.nlm.nih.gov/pubmed/27742763>
