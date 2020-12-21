#set input vars

export MSSQL_driver=/sc/arion/work/kumara22/covid_assembly_pipeline/drivers/libmsodbcsql-17.2.so.0.1
export nextstrain_dir=/sc/arion/projects/PVI/analyses/nextstrain_ncov2
export nextstrain_meta=/sc/arion/projects/PVI/analyses/nextstrain_ncov2/data/metadata_2020-11-16_11-06.tsv
export PVI_fasta=/sc/arion/projects/PVI/analyses/nextstrain_ncov2/data/sinai_seq.fasta
nextstrain_gisaid_fasta=/sc/arion/projects/PVI/analyses/nextstrain_ncov2/data/sequences_2020-11-16_07-46.fasta
nextstrain_output_fasta=/sc/arion/projects/PVI/analyses/nextstrain_ncov2/data/sequences.fasta
export admit_status_report=/sc/arion/projects/PVI/analyses/nextstrain_ncov2/data

module load mysql/5.7.19
module load R

Rscript scripts/nextstrain_tree_input.R

cat $nextstrain_gisaid_fasta $PVI_fasta > $nextstrain_output_fasta



