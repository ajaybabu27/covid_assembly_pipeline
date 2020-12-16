rule all:
	input:
		vadr_error_file=os.environ.get("sample_folder")+'/pipeline/VADR/VADR.vadr.fail.tbl'

rule assemble:
	message: "Assembling SARS-CoV-2 genome"
	input:
		sample_folder=os.environ.get("sample_folder")
	params:		
		sample_name=os.environ.get("sample_name"),
		run_ID=os.environ.get("run_ID")
	output:
		#consensus_fasta_file="{input.sample_folder}/pipeline/{params.sample_name}.fasta"
		consensus_fasta_file=os.environ.get("sample_folder")+'/pipeline/'+os.environ.get("sample_name")+'.fasta'
	shell:		
		"""
		if [ -d {input.sample_folder}/bams ]; then
			python scripts/run_pipeline.py -b {input.sample_folder}
		else
			python scripts/run_pipeline.py -i {input.sample_folder} -s {params.sample_name} -r1 _R1_001.fastq.gz -r2 _R2_001.fastq.gz			
		fi
		"""

rule variant_analysis:
	message: "Perform intra host variant analysis on SARS-CoV-2 library"
	input:
	    sample_folder=os.environ.get("sample_folder"),
		consensus_fasta_file=os.environ.get("sample_folder")+'/pipeline/'+os.environ.get("sample_name")+'.fasta'
	output:
		pileup_file=os.environ.get("sample_folder")+'/variants/'+'pileup'
	shell:
		"""
		if [ -d {input.sample_folder}/bams ]; then
			python scripts/variant_analysis.py -b {input.sample_folder}
		else
			python scripts/variant_analysis.py -i {input.sample_folder}
		fi
		"""

rule QC_analysis:
	message: "Perform QC analysis of the SARS-CoV-2 library"
	input:
	    sample_folder=os.environ.get("sample_folder"),
		pileup_file=os.environ.get("sample_folder")+'/variants/'+'pileup'
	params:		
		sample_name=os.environ.get("sample_name"),
		run_ID=os.environ.get("run_ID")
	output:
		qc_file=os.environ.get("sample_folder")+'/QC/'+'quality_control.pdf'
	shell:
		"""
		if [ -d {input.sample_folder}/bams ]; then
			python scripts/run_QC.py -b {input.sample_folder} -kdb /sc/arion/projects/PVI/db/minikraken_8GB_20200312
		else
			python scripts/run_QC.py -i {input.sample_folder} -kdb /sc/arion/projects/PVI/db/minikraken_8GB_20200312
		fi
		module load R
		Rscript plot-coverage-report.R -i {input.sample_folder}/variants/variable_bases.tsv -o {input.sample_folder}/variants/{params.sample_name}"_"var

		pdfunite {input.sample_folder}/QC/quality_control.pdf {input.sample_folder}/variants/{params.sample_name}"_"var.pdf {input.sample_folder}/QC/quality_control2.pdf
		mv {input.sample_folder}/QC/quality_control2.pdf {input.sample_folder}/QC/quality_control.pdf
		rm {input.sample_folder}/variants/{params.sample_name}"_"var.pdf

		"""
