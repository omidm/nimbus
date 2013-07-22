#!/bin/bash

scp -P 1361 -l 35000 grid.gz rfedkiw@bioxcluster:~/input

for i in `seq $1 $2`; 
do 

scp -P 1361 -l 35000 levelset_1.$i.* rfedkiw@bioxcluster:~/input;
scp -P 1361 -l 35000 levelset_2.$i.* rfedkiw@bioxcluster:~/input;
scp -P 1361 -l 35000 levelset_3.$i.* rfedkiw@bioxcluster:~/input;
scp -P 1361 -l 35000 levelset_4.$i.* rfedkiw@bioxcluster:~/input;

scp -P 1361  -l 35000 object_*.$i.* rfedkiw@bioxcluster:~/input;
scp -P 1361  -l 35000 rigid_*.$i.* rfedkiw@bioxcluster:~/input;

done

