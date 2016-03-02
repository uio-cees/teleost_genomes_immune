This folder contains all scripts and data needed to re-construct the estimation of MHC gene copy numbers as described in Supplementary Notes 22-28.

For a short video on how to proceed, see https://screencast.uninett.no/relay/ansatt/larssnnmbu.no/2016/25.01/385867/Teleost_copy_number_-_20160125_130342_39/index.html

1. Start R (or RStudio) and set the working directory to this folder.

2. Run the first script

  > source("script1_reference_markers.R")

...or open the script in RStudio, highlight chunks of code, and Run them step by step if you want to study the results and code in detail.

3. Run the second script

  > source("script2_reference_estimates.R")

...or (again) run chunk by chunk and inspect results/code if you like.

4. Run the third script

  > source("script3_MHC_estimates.R")

...or (again) run chunk by chunk and inspect results/code if you like.
