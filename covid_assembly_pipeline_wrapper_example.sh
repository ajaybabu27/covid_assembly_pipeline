cd <REPO DIR>

source activate covid
module load python/3.8.2

export sample_folder=$1
export sample_name=`basename $sample_folder`
export run_ID=`basename $(dirname "$sample_folder")`

snakemake --cores 12


