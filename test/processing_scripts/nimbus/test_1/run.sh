N=8
core_num=8
WD=../adjust_data
START_TIME=-3000
END_TIME=-2775
# Output the blocking timestamp of each job.
./output_blocking_timestamp.py $N ${WD}/{}_event_be.txt ${WD}/log_before_set ${WD}/blocking_time.txt
for i in `seq 1 $N`
do
  echo merging $i
  cat ${WD}/${i}_event_be.txt ${WD}/${i}_event_fe.txt > ${WD}/temp
  ./sort_base_on_first_number.py ${WD}/temp ${WD}/${i}.raw
done
./output_state_change.py $N ${WD}/{}.raw ${WD}/{}.state ${WD}/blocking_time.txt
for i in `seq 1 $N`
do
  ./sort_base_on_first_number.py ${WD}/${i}.state ${WD}/${i}.sort_state
done
./prepare_final.py ${N} ${WD}/{}.sort_state ${WD}/{}.prep_data
./make_all_figure.py ${core_num} ${N} ${WD}/{}.prep_data $START_TIME $END_TIME
./print_summary.py ${core_num} ${N} ${WD}/{}.prep_data
