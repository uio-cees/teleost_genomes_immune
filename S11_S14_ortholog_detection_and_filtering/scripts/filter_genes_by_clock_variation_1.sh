for i in ${1}*
do
    ruby resources/beauti.rb -id `basename $i` -n $i -u -s -o $i -c resources/root_constraints.xml
done
