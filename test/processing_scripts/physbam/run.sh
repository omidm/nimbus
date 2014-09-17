mkdir -p adjust_data/
./adjust_timedrifting.py 8 8 raw_data/mpi{}.log adjust_data/mpi{}.log
./make_all_figure.py 64 3 adjust_data/mpi{}.log
./print_summary.py 64 adjust_data/mpi{}.log
