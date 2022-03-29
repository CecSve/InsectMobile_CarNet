# Detecting flying insects using car nets and DNA metabarcoding #

This repository (R project) contains all scripts necessary to run the ecological/statistical analyses from the study [Detecting flying insects using car nets and DNA metabarcoding](https://royalsocietypublishing.org/doi/full/10.1098/rsbl.2020.0833) (Svenningsen et al., 2021).

The data were collected in June 2018 as part of the citizen science project InsectMobile ("Insektmobilen") at the Natural History Museum of Denmark. 

**Data used for analysis**: Svenningsen, Cecilie et al. (2021), Detecting flying insects using car nets and DNA metabarcoding, Dryad, Dataset, https://doi.org/10.5061/dryad.6q573n5z5

## Preparation of ASV tables and sequencing data ##

* Sequencing data was demultiplexed, filtered and cleaned with [DADA2](https://benjjneb.github.io/dada2/tutorial.html). 
* Erronous sequences were detected with [LULU](https://www.nature.com/articles/s41467-017-01312-x) and removed prior to analysis.
* Taxonomy was assigned to .fastas with the [GBIF sequence ID tool](https://www.gbif.org/tools/sequence-id).  

## Description of the sub directories for the sequence data processing ##

* **reports**: a step-wise list of scripts (01_, 02_ etc.) used for processing and analysing the data
