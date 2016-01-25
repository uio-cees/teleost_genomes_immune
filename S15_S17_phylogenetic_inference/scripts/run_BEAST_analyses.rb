# Michael Matschiner, 2015-04-08.

# Get the command line arguments.
replicates_directory = ARGV[0].chomp("/") + "/replicates"

# Collect names of nexus files in the input directory.
dir_entries_in = Dir.entries(replicates_directory)
replicate_directories = []
dir_entries_in.each {|e| replicate_directories << e if e.match(/^r\d+/)}

# Get the current directory.
home = Dir.pwd

# For each replicate directory.
replicate_directories.each do |r|

	# Change directory into the first replicate directory.
	Dir.chdir("#{replicates_directory}/#{r}")

	# Get the start command for the BEAST run.
	command = `cat start.sh | grep "java -jar"`
	system(command)

	# Change directory back home.
	Dir.chdir(home)

end