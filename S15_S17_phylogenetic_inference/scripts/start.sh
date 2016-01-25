# 1) Analyses of the concatenated dataset.

# Prepare BEAST input files for unconstrained (except the root) analyses for the concatenated nuclear dataset.
ruby prepare_concatenated_BEAST_analyses.rb ../data/nuclear ../analysis/nuclear_concatenated/beast2/unconstrained ../data/constraints/root_constraints.xml

# Run unconstrained (except the root) BEAST analyses for the concatenated nuclear dataset.
ruby run_BEAST_analyses.rb ../analysis/nuclear_concatenated/beast2/unconstrained

# Combine results of unconstrained (except the root) BEAST analyses for the concatenated nuclear dataset.
bash combine_BEAST_trees.sh ../analysis/nuclear_concatenated/beast2/unconstrained/replicates ../analysis/nuclear_concatenated/beast2/unconstrained/combined

# Prepare BEAST input files for fossil-constrained analyses for the concatenated nuclear dataset.
ruby prepare_concatenated_BEAST_analyses.rb ../data/nuclear ../analysis/nuclear_concatenated/beast2/fossil_constrained ../data/constraints/fossil_constraints.xml ../analysis/nuclear_concatenated/beast2/unconstrained/combined

# Produce five replicates of a fossil-constrained BEAST phylogeny for the concatenated nuclear dataset.
ruby run_BEAST_analyses.rb ../analysis/nuclear_concatenated/beast2/fossil_constrained

# Combine results of fossil-constrained BEAST analyses for the concatenated nuclear dataset.
bash combine_BEAST_trees.sh ../analysis/nuclear_concatenated/beast2/fossil_constrained/replicates ../analysis/nuclear_concatenated/beast2/fossil_constrained/combined


# 2) Identification of alignment bins, and analyses for binned datasets.

# Prepare statistical binning.
ruby prepare_statistical_binning.rb ../data/nuclear ../analysis/nuclear_binned/statistical_binning/01

# Run statistical binning.
bash run_statistical_binning.sh ../analysis/nuclear_binned/statistical_binning/01 ../analysis/nuclear_binned/statistical_binning/02 75

# Prepare BEAST input files for fossil-constrained analyses for each of the bins.
ruby prepare_binned_BEAST_analyses.rb ../data/nuclear ../analysis/nuclear_binned/statistical_binning/02 ../analysis/nuclear_binned/beast2/fossil_constrained_separate ../data/constraints/fossil_constraints.xml ../analysis/nuclear_concatenated/beast2/fossil_constrained/combined

# Run BEAST analyses for each bin (performed on abel).
bash run_binned_BEAST_analyses.sh ../analysis/nuclear_binned/beast2/fossil_constrained_separate

# Combine results from BEAST analyses for each bin.
bash combine_binned_BEAST_analyses.sh ../analysis/nuclear_binned/beast2/fossil_constrained_separate

# Prepare RAxML input files for each of the bins.
ruby prepare_binned_RAxML_analyses.rb ../data/nuclear ../analysis/nuclear_binned/statistical_binning/02 ../analysis/nuclear_binned/raxml

# Run RAxML analyses for the concatenated sequences of each bin.
bash run_binned_RAxML_analyses.rb ../analysis/nuclear_binned/raxml

# Run Astral analyses with the RAxML and BEAST trees for each bin.
bash run_astral_analyses_with_binned_trees.sh ../analysis/nuclear_binned/astral/input ../analysis/nuclear_binned/astral/output ../analysis/nuclear_binned/raxml ../analysis/nuclear_binned/beast2/fossil_constrained_separate


# 3) Gene tree analyses of unbinned alignments for individual genes.

# Run RAxML for each single gene alignment.
bash run_single_RAXML_analyses.sh ../data/nuclear ../analysis/nuclear_single/raxml

# Run Astral analyses.
bash run_astral_analyses_with_single_trees.sh ../analysis/nuclear_single/astral/input ../analysis/nuclear_single/astral/output ../data/nuclear ../../ortholog_identification/analysis/alignments/nuclear/09 ../analysis/nuclear_single/raxml

#  Root and ladderize all Astral trees.
rscript root_astral_trees.r ../analysis/nuclear_binned/astral/output/beast2.bins.tre ../analysis/nuclear_binned/astral/output/beast2.bins.rooted.tre ../analysis/nuclear_binned/astral/output/raxml.bins.tre ../analysis/nuclear_binned/astral/output/raxml.bins.rooted.tre ../analysis/nuclear_single/astral/output/beast2.single.tre ../analysis/nuclear_single/astral/output/beast2.single.rooted.tre ../analysis/nuclear_single/astral/output/raxml.single.tre ../analysis/nuclear_single/astral/output/raxml.single.rooted.tre
