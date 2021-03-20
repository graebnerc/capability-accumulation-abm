# Capability accumulation and product innovation: An agent-based perspective

[Claudius Gräbner](https://claudius-graebner.com/) and [Anna Hornykewycz](https://www.jku.at/en/institute-for-comprehensive-analysis-of-the-economy/about-us/team/anna-hornykewycz/)

This repository contains the code for the model used in the paper "Capability accumulation and product innovation: An agent-based perspective". 
The latter has been published as:

> Gräbner, C. and Hornykewycz, A. (2020): Capability accumulation and product innovation: An agent-based perspective, *ICAE Working Paper*, 
No. 108. URL: [https://www.jku.at/fileadmin/gruppen/108/ICAE_Working_Papers/wp108.pdf](https://www.jku.at/fileadmin/gruppen/108/ICAE_Working_Papers/wp108.pdf)

and is available [here](https://www.jku.at/fileadmin/gruppen/108/ICAE_Working_Papers/wp108.pdf) (OA).

## Folder structure of the repository

`R`: contains all R code, which is used to create the majority of the figures

`python`: contains all Python code; this includes all elements of the actual ABM; the directory `python/parameters` contains all the `.json` files with the parameter settings analyzed in the paper 

`bash`: contains `.sh` scripts that are mainly used for convenience and call a selection of `.py` and `.R` scripts sequentially or run the simulations for a number of different parameter constellations

`output`: contains all the output produced by the source code

`garbage`: contains only a single text file. Used by the script `clean_repo.sh` to clean the repository

## Replication of the paper results

You can run the scripts in the folder `bash` to run both the MCS and to produce
standard visualizations directly afterwards.
Make sure the script contains the right calls, though.

You can run all simulations used for the final paper in parallel using `GNU parallel`
by calling the following from the root repository:

```
zsh bash/run_all.sh
```

But keep in mind that this will consume a bit of computing power and will take quite a while. 
When these simulations are complete you can replicate the figures of the paper using the script `bash/make_paper_vis.sh`. 
Of course, you may also call selected functions to run specific parameter constellations of produce specific figures.

## System configuration used to produce paper results

The Python libraries contained in the virtual environment that has been
used to produce the results in the paper are summarized in the file
`environment.yml`. 
We used the [renv](https://rstudio.github.io/renv/index.html) package 
to keep track of the versions of the R-packages used. They are stored 
in `renv.lock`.
