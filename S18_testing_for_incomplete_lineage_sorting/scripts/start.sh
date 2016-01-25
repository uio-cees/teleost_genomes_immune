# Michael Matschiner, 2015-06-24

# Find all indels and the taxa affected by it.
ruby find_indels.rb ../data/02 ../analysis/indels/indels.txt 10 ../data/trees ../analysis/paup

# Run PAUP* to calculate branch lengths and retention indeces for each indel.
bash run_paup.sh

# Produce a table of branch lengths and percentage of indels with retention index 1.
ruby produce_branch_table.rb ../analysis/paup/76g_nucl_conc_fossils.combined.simple.log ../analysis/stats/76g_nucl_conc_fossils.combined.simple.branches.txt

# Produce a plot of log branch duration and proportion of indels with retention index 1.0 (reproduction of Jarvis et al. 2014, Figure 3F).
rscript produce_plot.r ../analysis/stats/76g_nucl_conc_fossils.combined.simple.branches.txt ../analysis/stats/76g_nucl_conc_fossils.combined.simple.branches.pdf