N=8
core_num=8
WD=../adjust_data
#Output the blocking timestamp of each job.
./output_blocking_timestamp.py $N ${WD}/{}_event_be.txt ${WD}/log_before_set ${WD}/blocking_time.txt
./prepare_final.py $N ${WD}/{}_event_fe.txt ${WD}/{}_event_be.txt ${WD}/blocking_time.txt
