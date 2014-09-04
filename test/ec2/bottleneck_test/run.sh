N=100
core_num=4
wd=./working_dir_scale_X
# merge/sort worker raw output. 
for i in `seq 1 $N`
do
  echo $i
  cat ${wd}/${i}_event_be.txt ${wd}/${i}_event_fe.txt > ${wd}/temp
  ./sort_base_on_first_number.py ${wd}/temp ${wd}/${i}.raw
done
# merge/sort scheduler raw output.
./sort_base_on_first_number.py ${wd}/job_assigner_log ${wd}/before_set.raw

for i in `seq 1 $N`
do
  echo $i
  ./output_blame.py $core_num ${wd}/${i}.raw ${wd}/${i}.step_1
done

for i in `seq 1 $N`
do
  echo $i
  ./sort_base_on_first_number.py ${wd}/${i}.step_1 ${wd}/${i}.step_2
done
./prune_blame.py ${N} ${wd}/%d_event_be.txt ${wd}/%d.step_2 ${wd}/before_set.raw ${wd}/%d.out
