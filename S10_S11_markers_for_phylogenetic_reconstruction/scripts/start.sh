# This pipeline requires the following tools (tested versions in parentheses):

# Uncompress the file of ENSEMBL gene trees.
gunzip ../data/Compara.78.protein.nhx.emf.gz 

exit
# Prepare the zebrafish query sequences for blast searches.
ruby prepare_queries.rb
