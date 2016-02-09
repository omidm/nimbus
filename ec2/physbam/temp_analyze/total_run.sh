for i in 256 344 400 440 480 512
do
    echo $i
    ./analyze_summary.py 64 ../output-scale${i}/mpi{}.log ../output-scale${i}/loop > summary_${i}.txt
done
