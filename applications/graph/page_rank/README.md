
# PageRank

The PageRank algorithm implemented for edge-cut partitioned graphs.
You can build the application library by issuing `make` in the current
directory. The application can be run by launching a controller and worker and
passing the page rank library to the worker, as follows

    $ cd <nimbus-root>
    $ ./scripts/start-controller.sh
    $ ./scripts/start-workers.sh 1 --fg -l applications/graph/page_rank/libpagerank.so



You can pass options to the application, as follows:

    $ ./scripts/start-workers.sh 1 --fg -l applications/graph/page_rank/libpagerank.so
              -h [ --help ]                         produce help message
              -g [ --input_dir ] arg (=../../applications/graph/page_rank/input-sample)
                                                    input directory that holds partitioned 
                                                    graph (absolute path or relative to 
                                                    nimbus_worker executable located at 
                                                    <nimbus-root>/nodes/nimbus_worker)
              -o [ --output_dir ] arg (=output)     output directory that saves the ranks 
                                                    (absolute path or relative to 
                                                    nimbus_worker executable located at 
                                                    <nimbus-root>/nodes/nimbus_worker)
              -i [ --iteration ] arg (=10)          number of iterations

The input files to page rank include the graph partitions, each with sets of
nodes and edges grouped as a single logical object. The application is written
to match the output files that the `graph_partitioner` generates. There are two
sample input files in this directory. Both partition the wikispeedia graph in
`graph_samples` folder into 10 pieces:

1. `input-sample': random partitioning (the default input)
2. `input-sample-metis': partitioned by METIS

For more information on the input format of the page rank algorithm and how to
generate the graph partitions, please refer to the README file in the
`graph_partitioner` folder.


