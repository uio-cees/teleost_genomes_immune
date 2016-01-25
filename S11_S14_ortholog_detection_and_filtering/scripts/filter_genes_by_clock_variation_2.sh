for i in ${1}ENSDARG*
do
    cp resources/beast.jar $i
done

for i in ${1}ENSDARG000000*
do
    cd $i
    trees_file=`basename ${i}.trees`
    if [ -f $trees_file ]
    then
        echo "Skipping directory ${i} as file ${trees_file} exists."
    else
        java -jar -Xmx4g beast.jar -threads 1 -seed $RANDOM -beagle `basename ${i}.xml`
    fi
    cd ../../../../scripts
done
