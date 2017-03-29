
# Sample Graphs

## wikipedia

This folder has the scripts to generate the graph meta data from wikidump. You
van download the wikidump online. The pipeline line is to first trim the
wikidump with the `trim_wikidump.py` script:


    $ ./trim_wikidump.py <pages-file> <links-file>

This will create a folder names `result/trim` with two files, `pages` and
`linkes`, in it.  Next step is to feed these files into the
`generate_wikigraph.py` script, as follows:

    $ ./generate_wikigraph.py -p <trimmed-page-file>
                              -l <trimmed-link-file>
                              -n <output-node-file-name>
                              -e <output-edges-file-name>
                              -j <output-node-number-file-name>
                              -k <output-edge-number-file-name>

The final results would be a four files as follows:

    1. nodes: each line is an integer.
    2. edges: each line is two integers specifying the source and destination of
              the edge.
    3. num_nodes: an integer specifying the number of nodes.
    4. num_edges: an integer specifying the number of edges.

Since the files are big, they are not included in the repository. But, you can
find a compressed version of the files (1.5GB) in the nimbus-data repo.


## wikispeedia

Smaller version of wikipedia with two files:

1. articles.tsv: each line is a name of a article.
2. links.tsv: each line is to article names specifying the source and
   destination of the edge. 

