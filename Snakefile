rule all:
	input:
		pdb_upload_check=os.environ.get("sample_folder")+'/pdb_upload_complete.txt',
		gt_file=os.environ.get("sample_folder")+'/genotype/'+'genotype_report.csv'

rule assemble:
	message: "Assembling SARS-CoV-2 genome"
	input:
		sample_folder=os.environ.get("sample_folder")
	params:		
		sample_name=os.environ.get("sample_name"),
		repo_dir=os.environ.get("repo_dir")		
	output:		
		consensus_fasta_file=os.environ.get("sample_folder")+'/pipeline/'+os.environ.get("sample_name")+'.fasta'
	shell:		
		"""
		if [ -d {input.sample_folder}/bams ]; then
			python scripts/run_pipeline.py -rd {params.repo_dir} -b {input.sample_folder}
		else
			python scripts/run_pipeline.py -rd {params.repo_dir} -i {input.sample_folder} -s {params.sample_name} -r1 _R1_001.fastq.gz -r2 _R2_001.fastq.gz			
		fi
		"""

rule variant_analysis:
	message: "Perform intra host variant analysis on SARS-CoV-2 library"
	input:
	    sample_folder=os.environ.get("sample_folder"),
		consensus_fasta_file=os.environ.get("sample_folder")+'/pipeline/'+os.environ.get("sample_name")+'.fasta'
	params:		
		repo_dir=os.environ.get("repo_dir")
	output:
		pileup_file=os.environ.get("sample_folder")+'/variants/'+'pileup'
	shell:
		"""
		if [ -d {input.sample_folder}/bams ]; then
			python scripts/variant_analysis.py -rd {params.repo_dir} -b {input.sample_folder}
		else
			python scripts/variant_analysis.py -rd {params.repo_dir} -i {input.sample_folder}
		fi
		"""

rule QC_analysis:
	message: "Perform QC analysis of the SARS-CoV-2 library"
	input:
	    sample_folder=os.environ.get("sample_folder"),
		pileup_file=os.environ.get("sample_folder")+'/variants/'+'pileup'
	params:		
		sample_name=os.environ.get("sample_name"),
		repo_dir=os.environ.get("repo_dir"),
		run_ID=os.environ.get("run_ID")
	output:
		qc_file=os.environ.get("sample_folder")+'/QC/'+'quality_control.pdf'
	shell:
		"""
		if [ -d {input.sample_folder}/bams ]; then
			python scripts/run_QC.py -rd {params.repo_dir} -b {input.sample_folder} -kdb /sc/arion/projects/PVI/db/minikraken_8GB_20200312
		else
			python scripts/run_QC.py -rd {params.repo_dir} -i {input.sample_folder} -kdb /sc/arion/projects/PVI/db/minikraken_8GB_20200312
		fi
		module load R
		Rscript scripts/plot-coverage-report.R -i {input.sample_folder}/variants/variable_bases.tsv -o {input.sample_folder}/variants/{params.sample_name}"_"var

		pdfunite {input.sample_folder}/QC/quality_control.pdf {input.sample_folder}/variants/{params.sample_name}"_"var.pdf {input.sample_folder}/QC/quality_control2.pdf
		mv {input.sample_folder}/QC/quality_control2.pdf {input.sample_folder}/QC/quality_control.pdf
		rm {input.sample_folder}/variants/{params.sample_name}"_"var.pdf
		"""

rule vadr_analysis:
	message: "Perform VADR analysis on the SARS-CoV-2 genome"
	input:
	    consensus_fasta_file=os.environ.get("sample_folder")+'/pipeline/'+os.environ.get("sample_name")+'.fasta'
	params:
		sample_folder=os.environ.get("sample_folder"),	
		sample_name=os.environ.get("sample_name")		
	output:
		vadr_error_file=os.environ.get("sample_folder")+'/pipeline/VADR/VADR.vadr.fail.tbl'
	shell:
		"python scripts/vadr_run.py {input.consensus_fasta_file} {params.sample_folder}/pipeline/VADR {params.sample_folder}/pipeline/VADR/VADR.gff /sc/arion/projects/PVI/db/vadr-models-corona-1.1-1"
		#python scripts/vadr_run.py os.environ.get("sample_folder")+"/pipeline/"+os.environ.get("sample_name")+".fasta" os.environ.get("sample_folder")+"/pipeline/VADR" os.environ.get("sample_folder")+"/pipeline/VADR/VADR.gff" /sc/arion/projects/PVI/db/vadr-models-corona-1.1-1

rule genotyping:
	message: "Perform Nextclade analysis to genotype the SARS-CoV-2 genome"
	input: 
	    consensus_fasta_file=os.environ.get("sample_folder")+'/pipeline/'+os.environ.get("sample_name")+'.fasta'
	params:
		sample_folder=os.environ.get("sample_folder"),	
		sample_name=os.environ.get("sample_name")
	output:
	    gt_file=os.environ.get("sample_folder")+'/genotype/'+'genotype_report.csv'
	shell:
	    """
	    #bash scripts/run-nextclade.sh {input.consensus_fasta_file} {params.sample_folder}/genotype/nextclade.tsv
	    nextclade -i {input.consensus_fasta_file} -t {params.sample_folder}/genotype/nextclade.tsv
	    python scripts/nextclade_parse.py {params.sample_folder}/genotype/nextclade.tsv {params.sample_folder}/genotype/genotype_report.csv	    	    
	    """

rule push_data_pathogendb:
	message: "Push genome assembly data to pathogenDB"
	input:
		qc_file=os.environ.get("sample_folder")+'/QC/'+'quality_control.pdf'
	params:
		sample_folder=os.environ.get("sample_folder"),
		sample_name=os.environ.get("sample_name"),
		run_ID=os.environ.get("run_ID").split("_")[0]
	output:
		pdb_upload_check=os.environ.get("sample_folder")+'/pdb_upload_complete.txt'
	shell:
		"""	
		module purge
		module load python/2.7.16
		python scripts/push_pathogendb.py {params.sample_name} {params.run_ID}
		touch {params.sample_folder}/pdb_upload_complete.txt
		"""
