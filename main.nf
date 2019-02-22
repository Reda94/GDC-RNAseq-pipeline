#!/usr/bin/env nextflow

input_dir = params.in //input directory containing all fastq files
input_data = "${input_dir}/*.gz"
fastq_files = Channel.fromPath(input_data).map{file -> tuple(file.simpleName, file)}

output_dir = params.out //output directory containing results for all files
threads = params.threads

gdc_ref_gen_idx_dir = params.ref_gen_idx_dir
gdc_ref_gen = params.ref_gen
gdc_gtf = params.gtf
rnaseqc_gtf = params.rnaseqc_gtf

path_to_data = "\$HOME/.nextflow/assets/Reda94/GDC-RNAseq-pipeline/data/"

//First alignment pass:

process alignmentFirstPass {

  publishDir "${output_dir}/${sample_name}"

  input:
  set sample_name, file(fq) from fastq_files

  output:
  set sample_name, file(fq), file("${sample_name}_first_pass/*SJ.out.tab") into firstpass_results

  """
  module load STAR/2.4.2a-foss-2018b
  mkdir ${sample_name}_first_pass

  STAR \\
  --genomeDir $gdc_ref_gen_idx_dir \\
  --readFilesIn $fq \\
  --outFileNamePrefix ./${sample_name}_first_pass/${sample_name} \\
  --runThreadN $threads \\
  --outFilterMultimapScoreRange 1 \\
  --outFilterMultimapNmax 20 \\
  --outFilterMismatchNmax 10 \\
  --alignIntronMax 500000 \\
  --alignMatesGapMax 1000000 \\
  --sjdbScore 2 \\
  --alignSJDBoverhangMin 1 \\
  --genomeLoad NoSharedMemory \\
  --readFilesCommand zcat \\
  --outFilterMatchNminOverLread 0.33 \\
  --outFilterScoreMinOverLread 0.33 \\
  --sjdbOverhang 100 \\
  --outSAMstrandField intronMotif \\
  --outSAMtype None \\
  --outSAMmode None
  """
}

//Intermediate indexing:

process intermediateIdx {

  publishDir "${output_dir}/${sample_name}"

  input:
  set sample_name, file(fq_sec_pass), file(sj) from firstpass_results

  output:
  set sample_name, file(fq_sec_pass), file(out_inter_idx) into inter_idx_results

  script:
  out_inter_idx = "${sample_name}_inter_idx"
  gd = "$out_inter_idx/"
  """
  module load STAR/2.4.2a-foss-2018b
  mkdir $out_inter_idx

  STAR \\
  --runMode genomeGenerate \\
  --genomeDir $gd \\
  --genomeFastaFiles $gdc_ref_gen \\
  --sjdbOverhang 100 \\
  --runThreadN $threads \\
  --sjdbFileChrStartEnd $sj
  """
}

//Second alignment pass:

process alignmentSecondPass {

  publishDir "${output_dir}/${sample_name}"

  input:
  set sample_name, file(fq_sp), file(idx) from inter_idx_results

  output:
  set sample_name, file("${sample_name}_second_pass/*Aligned.sortedByCoord.out.bam") into final_bam
  set sample_name, file("${sample_name}_second_pass/*Aligned.sortedByCoord.out.bam.bai") into final_bam_idx
  set sample_name, file("${sample_name}_second_pass/*Log.final.out") into fb_log_final_out
  set sample_name, file("${sample_name}_second_pass/*.out") into fb_out
  set sample_name, file("${sample_name}_second_pass/*SJ.out.tab") into fb_sj_out_tab
  set sample_name, file("${sample_name}_second_pass/*Log.out") into fb_log_out
  set sample_name, file("${sample_name}_second_pass/*Aligned.sortedByCoord.out.bam"), file("${sample_name}_second_pass/*Aligned.sortedByCoord.out.bam.bai") into rseqc_input

  """
  module load STAR/2.4.2a-foss-2018b
  module load SAMtools/1.1-foss-2018b
  mkdir ${sample_name}_second_pass

  STAR \\
  --genomeDir $idx \\
  --readFilesIn $fq_sp \\
  --outFileNamePrefix ./${sample_name}_second_pass/${sample_name} \\
  --runThreadN $threads \\
  --outFilterMultimapScoreRange 1 \\
  --outFilterMultimapNmax 20 \\
  --outFilterMismatchNmax 10 \\
  --alignIntronMax 500000 \\
  --alignMatesGapMax 1000000 \\
  --sjdbScore 2 \\
  --alignSJDBoverhangMin 1 \\
  --genomeLoad NoSharedMemory \\
  --limitBAMsortRAM 0 \\
  --readFilesCommand zcat \\
  --outFilterMatchNminOverLread 0.33 \\
  --outFilterScoreMinOverLread 0.33 \\
  --sjdbOverhang 100 \\
  --outSAMstrandField intronMotif \\
  --outSAMattributes NH HI NM MD AS XS \\
  --outSAMunmapped Within \\
  --outSAMtype BAM SortedByCoordinate \\
  --outSAMheaderHD @HD VN:1.4 \\
  --outSAMattrRGline ID:${sample_name} SM:${sample_name}

  samtools index ./${sample_name}_second_pass/${sample_name}Aligned.sortedByCoord.out.bam
  """
}

//Raw read counting:

process rawReadCount {

  publishDir "${output_dir}/${sample_name}"

  input:
  set sample_name, file(bam) from final_bam

  output:
  set sample_name, file("${sample_name}_raw_read_counts/*") into raw_counts_results

  """
  module load SAMtools/1.1-foss-2018b
  module load HTSeq/0.6.1p1-foss-2016b-Python-2.7.12
  mkdir ${sample_name}_raw_read_counts

  samtools view -F 4 $bam | \\
  htseq-count \\
  -m intersection-nonempty \\
  -i gene_id \\
  -r pos \\
  -s no \\
  - $gdc_gtf \\
  | grep 'ENS' > ./${sample_name}_raw_read_counts/${sample_name}raw_counts.txt
  """
}

process FPKM_TPM {

  publishDir "${output_dir}/${sample_name}"

  input:
  set sample_name, file(raw_counts) from raw_counts_results

  output:
  set sample_name, file("${sample_name}_FPKM_TPM/*.FPKM.txt"), file("${sample_name}_FPKM_TPM/*.TPM.txt") into FPKM_TPM_results

  """
  module load Python/3.6.6-foss-2018b
  mkdir ${sample_name}_FPKM_TPM

  FPKM_script.py $raw_counts ./${sample_name}_FPKM_TPM/${sample_name} $path_to_data
  """
}

process RNASeQC {

  publishDir "${output_dir}/${sample_name}"

  input:
  set sample_name, file(bam), file(bai) from rseqc_input

  output:
  set sample_name, file("${sample_name}_QC/*") into qc_results

  """
  module load RNA-SeQC/1.1.8-Java-1.7.0_80
  mkdir ${sample_name}_QC

  java -jar \${EBROOTRNAMINSEQC}/RNA-SeQC_v1.1.8.jar -o ${sample_name}_QC -r $gdc_ref_gen -s \"${sample_name}|$bam|notes\" -t $rnaseqc_gtf
  """
}
