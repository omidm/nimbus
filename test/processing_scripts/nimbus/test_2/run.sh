N=8
core_num=8
WD=../output-sep15
#Output the blocking timestamp of each job.
# ./output_blocking_timestamp.py $N ${WD}/{}_event_be.txt ${WD}/log_before_set ${WD}/blocking_time_2.txt
./prepare_final.py $N ${WD}/{}_event_fe.txt ${WD}/{}_event_be.txt ${WD}/blocking_time_2.txt
