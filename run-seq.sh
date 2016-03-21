#!/bin/bash
#
el_type=0
dom_size=3
#
#
path_to_example="examples/meshgenerator/cube_elasticity_fixed_bottom.txt"
#
searchedString0="nan"
#searchedString1="iter     |r|         r          e     time[s]"
searchedString1="Projector_l"
#
N=10000
cnt1=0
file_updated_each_sol="numberOfIter1.txt"
rm $file_updated_each_sol
#
for ((i = 1; i <= N; i++)); do
#
  d=$dom_size
   x=1; y=1; z=1
   X=1; Y=1; Z=1
   log_file=LOG-$X:$Y:$Z-$x:$y:$z-$d:$d:$d.log
   date | tee $log_file
   echo "Config: dom_size=  $d | cluster_size = $x:$y:$z | clusters = $X:$Y:$Z  " 
   date | tee -a $log_file
#
   mpirun  -n $(( X * Y * Z ))  ./espreso  $path_to_example \
          ${el_type[0]} ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d} | tee -a $log_file
#
  if  grep -qi $searchedString0 $log_file ;
    then
      echo "last i=$i"
    break 
  fi
#
  ii=$(grep "Projector_l" $log_file)
  arrj=(`echo ${ii}`);
  for j in $( seq 0 $((${#arrj[@]}-1)) )
  do
    if [ "${arrj[j]}" = "count:" ] ; then
      numberOfIter[cnt1]=${arrj[j+1]}
      echo "$cnt1/$N, it: ${arrj[j+1]} " >> $file_updated_each_sol 
      cnt1=$(($cnt1+1))
    fi
    continue
  done
#
done

echo "${numberOfIter[@]}" >  numberOfIter.txt

#for i in ${numberOfIter[@]}
#do
#  echo $i
#done
