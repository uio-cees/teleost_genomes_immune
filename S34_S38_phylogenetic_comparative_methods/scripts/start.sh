# Michael Matschiner, 2015-07-09

# Prepare a table with the sum of the mean copy numbers (U+Z lineage).
ruby prepare_copy_number_table.rb "../data/tables/MHCU_boot.txt" "../data/tables/MHCZ_boot.txt" "../data/tables/species.txt" "../analysis/pcm/input/log_copy_numbers.txt" "../analysis/pcm/input/taxa.txt"

# Prune the tree to remove taxa without copy number information.
rscript prune_tree.r "../data/trees/76g_nucl_conc_fossils.combined.tre" "../analysis/pcm/input/76g_nucl_conc_fossils.combined.pruned.tre" "../analysis/pcm/input/76g_nucl_conc_fossils.combined.pruned.nex" "../analysis/pcm/input/taxa.txt"

# Produce phenograms for log copy numbers.
bash produce_phenogram.sh "../analysis/pcm/input/76g_nucl_conc_fossils.combined.pruned.tre" "../analysis/pcm/input/log_copy_numbers.txt" "../analysis/phytools/output/phenogram.pdf"

# Run trait evolution analyses.
bash analyse_trait_evolution.sh "../analysis/pcm/input/76g_nucl_conc_fossils.combined.pruned.tre" "../analysis/pcm/input/log_copy_numbers.txt" "../analysis/pcm/output/plots.pdf" "../analysis/pcm/output/results.txt"

# Prepare slouch input files.
ruby prepare_slouch_input.rb "../analysis/pcm/input/76g_nucl_conc_fossils.combined.pruned.nex" "../analysis/pcm/input/log_copy_numbers.txt" "../analysis/slouch/input"

# Run slouch.
bash run_slouch.sh "../analysis/slouch/input" "../analysis/slouch/output"

# Produce a 3D figure based on the model inferred by slouch.
bash calculate_ancestral_states.sh "../analysis/pcm/input/76g_nucl_conc_fossils.combined.pruned.tre" "../data/tables/order.txt" "../analysis/slouch/input" "../analysis/slouch/output"

# Make summary table of Slouch runs.
ruby make_summary_table.rb "../analysis/slouch/output" "../../info/trait_evolution_models.txt"