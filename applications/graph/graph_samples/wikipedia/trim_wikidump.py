#!/usr/bin/env bash

# This script trims wiki dumps into only relevant information.
# The output is:
#   - results/pages_trimmed containing page id, namespace, title
#   - results/links_trimmed containing pade id (from), namespace, title (to)
# This should be further parsed by parse_and_save to construct edges as
# (src id, dst id) which can be parsed and partitioned into partitions and
# disjoint logical data objects.
# Author: chinmayee Shah

echo "Pages file is $1"
echo "Links file is $2"

RDIR=result/trim

mkdir -p ${RDIR}

echo ""
echo "Trimming pages ..."
grep INSERT $1 > ${RDIR}/pages
echo "Done trimming pages"

echo ""
echo "Trimming links ..."
grep INSERT $2 > ${RDIR}/links
echo "Done trimming links"
