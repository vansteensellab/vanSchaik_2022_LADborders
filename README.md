#### Computer scripts related to the LAD borders manuscript

This GitHub page contains all scripts that were used to create the figures for
the LAD borders manuscript with AID-mediated depletion of CTCF and proteins
involved in (loop) extrusion, by Tom van Schaik et al., 2021/2022 
(https://doi.org/10.1101/2021.09.13.460079). 

DamID processing was performed with a custom snakemake pipeline 
(damid_processing/bin/snakemake/damid.snake), using a configuration file with
information on the samples (damid_processing/bin/snakemake/config.yaml). The
other data was processed as described in the accompanying manuscript.

Except the initial data processing, all data analysis was performed in R using 
R-markdown scripts. These scripts are also included here and provide 
introductions and conclusions of the described analyses. It is important to 
highlight that these scripts were not filtered to only include the figure panels
that were used in the manuscript. This should provide a more transparant 
overview of our data analysis.