N=1
core_num=4
wd=./workdir
echo Merging events.
for i in `seq 1 $N`
do
  echo $i
  cat ${wd}/${i}_event_be.txt ${wd}/${i}_event_fe.txt > ${wd}/temp
  ./sort_base_on_first_number.py ${wd}/temp ${wd}/${i}.raw
done

for i in `seq 1 $N`
do
  echo $i
  ./output_blame.py $core_num ${wd}/${i}.raw ${wd}/${i}.out
done
