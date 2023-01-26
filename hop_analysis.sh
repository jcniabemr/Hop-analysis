#!/usr/bin/env bash
#SBATCH -p medium
#SABTCH -J hop
#SBATCH --mem=2G
#SBATCH --cpus-per-task=4

####Script to assemble genomes and conduct analysis

####Run FASTQC
# for RawData in /home/jconnell/hop/A111055/*.fq.gz; do
#     OutDir=/home/jconnell/hop/A111055_fastqc
#     mkdir -p $OutDir
#     ProgDir=/home/jconnell/hop/scripts
#     sbatch $ProgDir/fastqc_ud.sh $RawData $OutDir
# done
  
####Collect QC data
# ls -l /home/jconnell/niab/andrea_rna_seq | awk '{print $9}' | grep -v "scripts\|result_table_untrim.txt" > filenames
# ls -l /home/jconnell/niab/andrea_rna_seq | grep -v "scripts\|result_table_untrim.txt" | awk '{print $9 "-F"}' > newnames_f
# ls -l /home/jconnell/niab/andrea_rna_seq | grep -v "scripts\|result_table_untrim.txt" | awk '{print $9 "-R"}' > newnames_r
# sed -i 1d newnames_f
# sed -i 1d newnames_r
# #Get read counts 
# for name in $(cat filenames); do 
# 	for file in /home/jconnell/niab/andrea_rna_seq/$name/*_1_fastqc; do 
# 		cat $file/fastqc_data.txt | sed -n 7p | awk '{print $3}' >>  f_read_counts_pre
# 		cat $file/fastqc_data.txt | sed -n 10p | awk '{print $2}' >>  f_GC_counts_pre
# 		cat $file/summary.txt | cut -f1 | paste -s >> result_f_pre
# 	done 
# done 
# for name in $(cat filenames); do 
# 	for file in /home/jconnell/niab/andrea_rna_seq/$name/*_2_fastqc; do 
# 		cat $file/fastqc_data.txt | sed -n 7p | awk '{print $3}' >>  r_read_counts_pre
# 		cat $file/fastqc_data.txt | sed -n 10p | awk '{print $2}' >>  r_GC_counts_pre
# 		cat $file/summary.txt | cut -f1 | paste -s >> result_r_pre
# 	done 
# done  
# paste newnames_f f_read_counts_pre f_GC_counts_pre result_f_pre >> final_result_f_pre
# paste newnames_r r_read_counts_pre r_GC_counts_pre result_r_pre >> final_result_r_pre
# cat final_result_r_pre >> final_result_f_pre
# echo -e "Read name"'\t'"Read count"'\t'"GC content" > headder1
# cat /home/jconnell/niab/andrea_rna_seq/$name/*_2_*/summary.txt | cut -f2 | paste -s > headder2
# paste headder1	headder2 > headders
# sort -V final_result_f_pre > final_result
# cat final_result >> headders	
# mv headders	/home/jconnell/niab/andrea_rna_seq/result_table_untrim.txt
# rm newnames_f filenames newnames_r f_read_counts_pre r_read_counts_pre f_GC_counts_pre r_GC_counts_pre result_f_pre	result_r_pre final_result_f_pre	final_result_r_pre headder1 headder2 headders final_result	

####Trim_reads
#  IlluminaAdapters=/mnt/shared/scratch/agomez/apps/git_repos/bioinformatics_tools/SEQdata_qc/illumina_full_adapters.fa
#   for StrainPath in /home/jconnell/hop/A111055; do
#    ReadsF=$(ls $StrainPath/*.fq.gz | grep "_1.fq.gz")
#    ReadsR=$(ls $StrainPath/*.fq.gz | grep "_2.fq.gz")
#    F_outdir=/home/jconnell/hop/A111055_trimmed_reads/F
#    R_outdir=/home/jconnell/hop/A111055_trimmed_reads/R
#    mkdir -p $F_outdir $R_outdir
#    ProgDir=/home/jconnell/niab/andrea_rna_seq/scripts
#    sbatch $ProgDir/fastq_mcf.sh $ReadsF $ReadsR $IlluminaAdapters $F_outdir $R_outdir
#   done   

####Run fastQC on trimmed reads
# v=$(ls -l /mnt/shared/projects/niab/pseudomonas_RNAseq/bacteria/*/01.RawData | awk '{print $9}')
# for f in $(echo $v); do
#   for RawData in $(ls /home/jconnell/projects/niab/andrea_rna_seq/trimmed_reads/${f}/*/*.fq.gz); do
#       OutDir=/home/jconnell/niab/andrea_rna_seq/trimmed_readsQC/${f}
#       mkdir -p $OutDir
#       ProgDir=/home/jconnell/projects/niab/andrea_rna_seq/scripts
#       sbatch $ProgDir/fastqc_ud.sh $RawData $OutDir
#   done
# done  

####Collect QC data past trim 
# ls -l /home/jconnell/niab/andrea_rna_seq/trimmed_readsQC | awk '{print $9}' > filenames
# ls -l /home/jconnell/projects/niab/andrea_rna_seq/trimmed_reads | awk '{print $9 "-F"}' > newnames_f
# ls -l /home/jconnell/projects/niab/andrea_rna_seq/trimmed_reads | awk '{print $9 "-R"}' > newnames_r
# sed -i 1d newnames_f
# sed -i 1d newnames_r
# #Get read counts 
# for name in $(cat filenames); do 
# 	for file in /home/jconnell/niab/andrea_rna_seq/trimmed_readsQC/$name/*_1_trim_fastqc; do 
# 		cat $file/fastqc_data.txt | sed -n 7p | awk '{print $3}' >>  f_read_counts_pre
# 		cat $file/fastqc_data.txt | sed -n 10p | awk '{print $2}' >>  f_GC_counts_pre
# 		cat $file/summary.txt | cut -f1 | paste -s >> result_f_pre
# 	done 
# done 
# for name in $(cat filenames); do 
# 	for file in /home/jconnell/niab/andrea_rna_seq/trimmed_readsQC/$name/*_2_trim_fastqc; do 
# 		cat $file/fastqc_data.txt | sed -n 7p | awk '{print $3}' >>  r_read_counts_pre
# 		cat $file/fastqc_data.txt | sed -n 10p | awk '{print $2}' >>  r_GC_counts_pre
# 		cat $file/summary.txt | cut -f1 | paste -s >> result_r_pre
# 	done 
# done  
# paste newnames_f f_read_counts_pre f_GC_counts_pre result_f_pre >> final_result_f_pre
# paste newnames_r r_read_counts_pre r_GC_counts_pre result_r_pre >> final_result_r_pre
# cat final_result_r_pre >> final_result_f_pre
# echo -e "Read name"'\t'"Read count"'\t'"GC content" > headder1
# cat /home/jconnell/niab/andrea_rna_seq/trimmed_readsQC/$name/*_2_*/summary.txt | cut -f2 | paste -s > headder2
# paste headder1	headder2 > headders
# sort -V final_result_f_pre > final_result
# cat final_result >> headders	
# mv headders	/home/jconnell/niab/andrea_rna_seq/result_table_trim_QC.txt
# rm newnames_f filenames newnames_r f_read_counts_pre r_read_counts_pre f_GC_counts_pre r_GC_counts_pre result_f_pre	result_r_pre final_result_f_pre	final_result_r_pre headder1 headder2 headders final_result	

####Assemble genomes 
# for StrainPath in /home/jconnell/hop/A111055_trimmed_reads; do
#     F_Read=$(ls $StrainPath/F/*_1_*)
#     R_Read=$(ls $StrainPath/R/*_2_*)
#     Organism=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
#     OutDir=/home/jconnell/hop/spades
#     mkdir -p $OutDir
#     ProgDir=/home/jconnell/hop/scripts
#     sbatch $ProgDir/slurm_spades_30cpu.sh $F_Read $R_Read $OutDir correct 10
# done

####Run Busco
# for Assembly in /home/jconnell/hop/genomes/A111055/A111055_contigs_unmasked.fasta; do
#     ProgDir=/home/jconnell/git_repos/emr_repos/Fv_C-variants/Busco
#     b=$(basename $Assembly .fasta)
#     BuscoDB=/home/jconnell/hop/busco_database/hypocreales_odb10
#     OutDir=/home/jconnell/hop/A111055/busco/
#     mkdir -p $OutDir
#     sbatch $ProgDir/busco.sh $Assembly $BuscoDB $OutDir $strain
# done

####Run checkM
# data=/home/jconnell/hop/genomes/A111055/A111055_contigs_unmasked.fasta
# temp_files=/home/jconnell/hop/temp_checkm
# mkdir -p $temp_files	
# outdir=/home/jconnell/hop/A111055/checkm
# cp --symbolic ${data} $temp_files
# scriptdir=/home/jconnell/git_repos/tools/analysis_tools
# sbatch $scriptdir/checkM.sh $temp_files $outdir

####Collect busco data
# data=/home/jconnell/hop/A111055/busco/A111055_contigs_unmasked
#   cat ${data}/*.txt | sed -n 9p | awk '{print $1}' | cut -c 3-8 | sed 's/\[//g' >> completeness 
#   cat ${data}/*.txt | sed -n 10p | awk '{print $1}' >> complete_buscos 
#   cat ${data}/*.txt | sed -n 11p | awk '{print $1}' >> complete_single_copy_buscos
#   cat ${data}/*.txt | sed -n 12p | awk '{print $1}' >> complete_and_duplicated_buscos 
#   cat ${data}/*.txt | sed -n 13p | awk '{print $1}' >> fragmented_buscos 
#   cat ${data}/*.txt | sed -n 14p | awk '{print $1}' >> missing_buscos 
#   cat ${data}/*.txt | sed -n 15p | awk '{print $1}' >> total_buscos_searched 
#   echo $x >> file_names
#   done 
# done 
# echo -e "Genome""\t""Busco % complete""\t""Complete BUSCOs""\t""Complete and single-copy BUSCOs""\t""Complete and duplicated BUSCOs""\t""Fragmented BUSCOs""\t""Missing BUSCOs""\t""Total BUSCO groups searched" > headder
# paste file_names completeness complete_buscos complete_single_copy_buscos complete_and_duplicated_buscos fragmented_buscos missing_buscos total_buscos_searched > res1
# cat res1 | sort -k2 -n >> headder 
# mv headder A111055_busco_results.txt
# rm completeness complete_buscos complete_single_copy_buscos complete_and_duplicated_buscos fragmented_buscos missing_buscos total_buscos_searched file_names res1 
