#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
@Created: 07-21-2020
@Description: This script updates assembly stats in the assembly table and the genome, cds and protein sequences in the assembly_sequences table in pathogendb
"""

from __future__ import division
from Bio import SeqIO
import mysql.connector
from collections import defaultdict
import sys
from itertools import chain
import os
import os.path
import shutil
from datetime import datetime

path = os.path.expanduser('~') + '/.my.cnf.pdbrw'
with open(path) as cnf_file:
	for line in cnf_file:
		if line.startswith('user='):
			user = line.rstrip()[5:]
		if line.startswith('password='):
			pw = line.rstrip()[9:]
		if line.startswith('host='):
			host = line.rstrip()[5:]
		if line.startswith('database='):
			database = line.rstrip()[9:]

db = mysql.connector.connect(host=host,user=user,passwd=pw,db=database)

old_run_list=['TD01619_TD01622_TD01623','TD01627','TD01628_TD01629','TD01630_TD01631','TD01635','TD01637_TD01638',
'TD01640','TD01640_TD01644_TD01654','TD01640_TD01644_TD01654_TD01753','TD01644','TD01646','TD01654','TD01655','TD01662','TD01686',
'TD01694','TD01697','TD01704','TD01719','TD01724','TD01735','TD01748','TD01753','TD01778','TD01785','TD01809','TD01833','TD01836','TD01836_rerun','TD01836_TD01836_rerun',
'TD01866','TD01866_rerun','TD01866_TD01866_rerun','TD01889','TD01916','TD01924']

sample=sys.argv[1]
systematic_id=sample
#systematic_id=sample.replace('-','_')
extract_id=sample.split('_')[1]
colloborator_id=sample.split('_')[0]
if sample.split('_')[0]=='Pool':
	sys.exit(0)
if extract_id=='neg':
	sys.exit(0)
run=sys.argv[2]
cur=db.cursor()	
assembly_check_cmd="Select assembly_id from tCEIRS_assemblies where Extract_ID='"+str(extract_id)+"' and assembly_run='"+str(run)+"'"
cur.execute(assembly_check_cmd)
data=cur.fetchall()
if len(data)>0:
	for row in data:		
		assembly_id=row[0]
	update='yes'
	#sys.exit(0)
else:
	update='no'
if colloborator_id=='SS':
    base_dir='/sc/arion/projects/CHARM/genomes/assembly'
else:
    base_dir='/sc/arion/projects/PVI/genomes/assembly'
try:	
	sample_path_qc=base_dir+'/'+run+'_pipeline/'+sample+'/QC'	
	sample_path_variants=base_dir+'/'+run+'_pipeline/'+sample+'/variants'	
	bacteria_reads=0
	eukaryote_reads=0
	coronavirus_reads=0
	iav_reads=0
	ibv_reads=0
	var15_count=0
	with open("%s/refbam.flagstat" % sample_path_qc) as f:		
		total_reads = int(f.readline().split()[0])
		f.readline()
		f.readline()
		f.readline()
		mapped_reads = int(f.readline().split()[0])
	with open("%s/kraken_report.out" % sample_path_qc) as f:
		for line in f:
			read_count_clade = int(line.split()[1])
			classification = line.split()[5]
			if classification=='Bacteria':
				bacteria_reads=read_count_clade
			if classification=='Eukaryota':
				eukaryote_reads=read_count_clade
			if classification=='Coronaviridae':
				coronavirus_reads=read_count_clade
			if classification=='Influenza A virus':
				iav_reads=read_count_clade
			if classification=='Influenza B virus':
				ibv_reads=read_count_clade
	with open("%s/variable_bases.tsv" % sample_path_variants) as f:
		f.readline()
		for line in f:						
			if line.split('\t')[2]=='FLAGGED':
				
				var15_count+=1
	if coronavirus_reads==0:
		print sample,'no corona det'
	coronavirus_reads=coronavirus_reads+mapped_reads
	
	unique_mapped_per=round((mapped_reads/total_reads)*100,2)
	short_unmapped_per=round((100-unique_mapped_per),2)
	bact_per=round((bacteria_reads/total_reads)*100,2)
	euk_per=round((eukaryote_reads/total_reads)*100,2)
	cor_per=round((coronavirus_reads/total_reads)*100,2)
	iav_per=round((iav_reads/total_reads)*100,2)
	ibv_per=round((ibv_reads/total_reads)*100,2)
	
	fasta_path=base_dir+'/'+run+'_pipeline/'+sample+'/pipeline/'+sample+".fasta"
	fastaSequence=SeqIO.read(fasta_path, "fasta")
	nCount = fastaSequence.seq.lower().count('n') 
	length = len(fastaSequence.seq)
	non_amb_length=length-nCount
	
	completeness=round((non_amb_length/29870.0),2)
	
	if completeness>1:
		completeness=1
        
	seq=fastaSequence.seq.upper()
	
	if completeness>=0.95:
		assembly_status='Complete'        
	else:
		assembly_status='Partial'
		

	if cor_per<4:        
		assembly_quality='Failed'
	else:
		if run in old_run_list:
			if var15_count>4:
				assembly_quality='Check'
			else:
				assembly_quality='Passed'
		else:
			if var15_count>9:
				assembly_quality='Check'
			else:
				assembly_quality='Passed'
            
	
	
	if update=='no':		
		cur=db.cursor()
		cur.execute("INSERT INTO `tCEIRS_assemblies` (`Extract_ID`, `assembly_run`, `assembly_status`, `Total_reads`,`Uniq_mapped_read_percent`,`Short_unmapped_read_percent`,`completeness`,`IAV_percent`,`IBV_percent`,`coronavirus_percent`,`Eukaryota_percent`,`Bacteria_percent`,`Assembly_quality`,`Variant_pos_sum_15pct`) VALUES ('%s','%s','%s','%i','%f','%f','%f','%f','%f','%f','%f','%f','%s','%i')" %(extract_id,run,assembly_status,total_reads,unique_mapped_per,short_unmapped_per,completeness,iav_per,ibv_per,cor_per,euk_per,bact_per,assembly_quality,var15_count))
		assembly_id=cur.lastrowid
		cur=db.cursor()
		cur.execute("INSERT INTO `tCEIRS_assembly_sequences`(`AssemblyID`, `Sequence_type`, `Sequence_name`, `Sequence`,`Sequence_length`,`Sequence_quality`,`Variant_pos_sum`) VALUES ('%s', 'Genome', 'SARS-CoV-2', '%s','%i','%s','%i')" % (assembly_id,seq,length,assembly_quality,var15_count))
		
		

	else:
		
		cur=db.cursor()
		cur.execute("UPDATE `tCEIRS_assemblies` SET assembly_status ='%s', completeness=%f  where assembly_ID='%i'" % (assembly_status,completeness,assembly_id))
		#cur=db.cursor()
		#cur.execute("UPDATE `tCEIRS_assembly_sequences` SET Sequence_quality ='%s', Variant_pos_sum=%i where AssemblyID='%i'" % (assembly_quality,var15_count,assembly_id))

	output_dir='/sc/arion/projects/CRIPT/crip-surveillance/www/crip-final-assemblies/'+str(assembly_id)
	os.system('rm -r '+output_dir)
	os.system('mkdir -p '+output_dir+';cp '+fasta_path+' '+output_dir+'/aid_'+str(assembly_id)+'.fa')
	
	os.system('cp '+sample_path_qc+'/quality_control.pdf'+' '+output_dir+'/'+str(systematic_id)+'_final.report.pdf')
	os.system('cp '+sample_path_variants+'/variable_bases.tsv'+' '+output_dir+'/'+str(systematic_id)+'_final.variants.calls.txt')
	
	vadr_path=base_dir+'/'+run+'/'+sample+'/pipeline/VADR'
	os.system('cp '+vadr_path+'/'+'VADR.vadr.pass.tbl'+' '+output_dir+'/'+str(systematic_id)+'_final.features_table.txt')
	os.system('cp '+vadr_path+'/'+'VADR.vadr.fail.tbl'+' '+output_dir+'/'+str(systematic_id)+'_final.features_table_errors.txt')



except IOError:
	print sample,'missing file',run

db.commit()
		
