# sTEP_overlap

This repository contains code for the 2018 global analysis of geographic, genetic, and trait data.  

Making this easily reproducible is a major challenge because the data files are huge >100 GB, and the analysis currently requires somewhere around 70 GB of RAM.  Moreover many of the best software tools for these large datafiles are not platform independent, which adds another layer of difficulty.  Nonetheless, the parts of the code in this repository may help.  

Currently the `Rakefile` downloads the raw data and processes it (rake will likely require some specific software installation depending on your platform; more about rake [here](https://github.com/ruby/rake)).  The the analysis is done via the R scripts in the R folder.  The major memory intensive step is the generalized addative model which is documented [here](https://www.rdocumentation.org/packages/mgcv/versions/1.8-24/topics/bam).  