N=100
ORIG_D=_save_output_21
WD=result
mkdir -p ${WD}
cp ${ORIG_D}/log_before_set ${WD}/
for i in `seq 1 $N`
do
	./trim.py ${ORIG_D}/${i}_event_fe.txt ${WD}/${i}_event_fe.txt
	./trim.py ${ORIG_D}/${i}_event_be.txt ${WD}/${i}_event_be.txt
done
