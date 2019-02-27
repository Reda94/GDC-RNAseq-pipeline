#!/usr/bin/env nextflow

input_dir = params.in //input directory containing all fastq files
input_data = "${input_dir}/*.htseq.counts"
htseq_files = Channel.fromPath(input_data).map{file -> tuple(file.simpleName, file)}

output_dir = params.out //output directory containing results for all files
path_to_data = "\$HOME/.nextflow/assets/Reda94/GDC-RNAseq-pipeline/data/"

process FPKM_TPM {

  publishDir "${output_dir}/${sample_name}"

  input:
  set sample_name, file(raw_counts) from htseq_files

  output:
  set sample_name, file("${sample_name}_FPKM_TPM/*.FPKM.txt"), file("${sample_name}_FPKM_TPM/*.TPM.txt") into FPKM_TPM_results

  """
  module load Python/3.6.6-foss-2018b
  mkdir ${sample_name}_FPKM_TPM

  FPKM_script.py $raw_counts ./${sample_name}_FPKM_TPM/${sample_name} $path_to_data
  """
}
