# Michael Matschiner, 2015-07-07

replicates_dir=$1
for rep in ${replicates_dir}/r???
do
	tree_file_name="${rep}/*.tre"
	mcmc_out_file_name="${rep}/mcmc_out.txt"
	event_data_file_name="${rep}/event_data.txt"
	prior_probs_file_name="${rep}/prior_probs.txt"
	rates_file_name="${rep}/rates.txt"
	post_probs_file_name="${rep}/post_probs.txt"
	netdiv_best_config_file_name="${rep}/best_netdiv.pdf"
	lambda_best_config_file_name="${rep}/best_lambda.pdf"
	mu_best_config_file_name="${rep}/best_mu.pdf"
	set_of_4_file_name="${rep}/set4.pdf"
	set_of_9_file_name="${rep}/set9.pdf"
	bayes_factors_file_name="${rep}/bayes_factors.pdf"
	cohorts_file_name="${rep}/cohorts.pdf"
	net_div_tree_file_name="${rep}/net_div_tree.pdf"
	branch_file_name="${rep}/branches.txt"
	branch2_file_name="${rep}/branches2.txt"
	rscript produce_bamm_plots.r ${tree_file_name} ${mcmc_out_file_name} ${event_data_file_name} ${prior_probs_file_name} ${rates_file_name} ${post_probs_file_name} ${netdiv_best_config_file_name} ${lambda_best_config_file_name} ${mu_best_config_file_name} ${set_of_4_file_name} ${set_of_9_file_name} ${bayes_factors_file_name} ${cohorts_file_name} ${net_div_tree_file_name} ${branch_file_name} ${branch2_file_name}
	exit
done
