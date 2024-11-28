This is the github repo for "Microbial reaction rate estimation using proteins and proteomes".

If you want to run the analysis, the key script is `running-all-analyses-plots.R`. This 
will take some time, so there is a flag in the script that is designated to reload key 
analyses.

To reload the analyses, you first need to decompress some of the files (they were too 
big to be uploaded raw to github). These are in `data/intermediate_data/`. To unzip 
them, you can use the script `scripts/decompress-intermediate-data.sh` or just manually 
with 7zip.

Before running the scripts, I recommend using 
[renv](https://rstudio.github.io/renv/articles/renv.html) to load the same version of R 
I used, and the same set of packages. For example, you can rerun this on a cluster 
using (within a bash script):

```
renv::restore()
```

And then using a Singularity image to run the analyses, after loading all the libraries 
as above, you can use the following bash script:

```
#!/bin/bash
#SBATCH --job-name=running
#SBATCH --nodes=1
#SBATCH --time=20:00:00
#SBATCH --mem=1gb
#SBATCH --output=running.out
#SBATCH --error=running.err

module load singularity/3.10.4

singularity exec --bind /net:/net docker://rocker/verse:4.2.2 Rscript 
running-all-analyses-plots.R
```
