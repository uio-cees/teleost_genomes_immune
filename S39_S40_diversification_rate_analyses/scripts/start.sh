# Michael Matschiner, 2015-09-09.

# Add Lepidogalaxias to the tree, with a divergence age chosen at random between the root and the divergence of Salsal.
bash add_lepidogalaxias_all.sh ../data/trees/76g_nucl_conc_fossils.combined.simple.ladder.tre ../analysis/bamm/replicates 10

# Prepare BAMM input.
ruby prepare_bamm_input.rb ../data/tables/species.txt ../data/tables/clades.txt ../data/templates/divcontrol.txt ../analysis/bamm/replicates 10 resources/bamm

# Run BAMM.
bash run_bamm.sh ../analysis/bamm/replicates

# Analyze BAMM output.
bash produce_bamm_plots_all.sh ../analysis/bamm/replicates

# Compare per branch net diversification rate and mean trait value.
ruby compare_per_branch_rates_and_traits.rb ../analysis/bamm/replicates/r001/branches.txt ../data/outfiles/11001011.3d.txt ../../info/traits_vs_net_diversification.txt ../../info/traits_vs_net_diversification.svg

# Produce a phylomorphospace plot for net diversification and mean trait value, per terminal branch.
ruby prepare_phylomorphospace_plot_input.rb ../data/trees/76g_nucl_conc_fossils.combined.pruned.tre ../analysis/bamm/replicates/r001/branches.txt ../data/outfiles/11001011.txt ../analysis/phytools/input/76g_nucl_conc_fossils.combined.pruned.tre ../analysis/phytools/input/tips.txt
bash produce_phylomorphospace_plot.sh ../analysis/phytools/input/76g_nucl_conc_fossils.combined.pruned.tre ../analysis/phytools/input/tips.txt ../analysis/phytools/output/phylomorphospace.pdf

# Prune the tree to a single tip per clade.
bash prune_tree.sh ../data/tables/species.txt ../data/tables/order.txt ../data/trees/76g_nucl_conc_fossils.combined.simple.ladder.tre ../analysis/diversitree/input/76g_nucl_conc_fossils.clades.tre

# Prepare a table of species richness values and trait values for each clade.
ruby prepare_diversitree_table.rb ../analysis/diversitree/input/76g_nucl_conc_fossils.clades.tre ../data/tables/clades.txt ../data/tables/species.txt ../data/outfiles/11001011.txt ../analysis/diversitree/input/table.txt

# Prepare diversitree input with different thresholds for high/low trait values.
ruby prepare_diversitree_input.rb ../analysis/diversitree/input/76g_nucl_conc_fossils.clades.tre ../analysis/diversitree/input/table.txt ../analysis/diversitree ../data/tables/species.txt ../analysis/bamm/replicates/r001/branches2.txt run_diversitree.r

# Run diversitree analyses.
bash run_diversitree.sh ../analysis/diversitree

# Summarize results of diversitree analyses.
ruby summarize_diversitree.rb ../analysis/diversitree ../analysis/diversitree/summary

# Plot diversitree results.
rscript plot_diversitree.r ../analysis/diversitree/summary
