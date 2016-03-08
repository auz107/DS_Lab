#!/bin/bash -l
#$-l h_rt=24:00:00

cd /usr2/postdoc/alizom/work/models/ModelSEED 

#$ -j y

#$ -o job_read_tsv_files.out

source activate /projectnb/bioinfor/SEGRE/alizom

python -c "import time;print '\n**Job started at ',time.strftime('%c'),'\n'" 

python  read_tsv_files.py 

python -c "import time;print '\n**Job ended at ',time.strftime('%c'),'\n'" 


