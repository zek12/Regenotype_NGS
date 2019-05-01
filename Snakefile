configfile: "config.yaml"

from os.path import *

rule all:
	input:
		[
		"DRG_037/DRG_037-F01/realigned_bams/DRG_037-F01_LIBx1_FINAL_b37_PE.aln.RG.dupes_marked.flagstat_cram"
		]



# rule sort_bam:
# 	input: config["root_dir"] + "/{run}/{sample}/{basename}.bam"
# 	output: "{run}/{sample}/bams_pre/{basename}.sorted.bam"
# 	# shell: "samtools sort -o {output} -n {input}"
# 	shell: "java -Xmx2g -jar " + config["path_picard"] + " SortSam SO=queryname I={input} O={output}"

# rule fixmate:
# 	input: "{run}/{sample}/bams_pre/{basename}.sorted.bam"
# 	output: "{run}/{sample}/bams_pre/{basename}.fixed_mate.bam"
# 	# shell: "samtools fixmate {input} {output}"
# 	shell: "java -Xmx2g -jar " + config["path_picard"] + " FixMateInformation I={input} O={output}"


# another option to the following step is samtools collate + samtools fastq
rule picard_SamToFastq:
	input: config["root_dir"] + "/{run}/{sample}/{basename}.bam"
	output:
		fq1="{run}/{sample}/regenerated_fastqs/{basename}.R1.fq.gz",
		fq2="{run}/{sample}/regenerated_fastqs/{basename}.R2.fq.gz",
		fq0="{run}/{sample}/regenerated_fastqs/{basename}.unpaired.fq.gz"
	shell: "java -Xmx2g -jar " + config["path_picard"] + " SamToFastq I={input} FASTQ={output.fq1} SECOND_END_FASTQ={output.fq2} UNPAIRED_FASTQ={output.fq0}"


rule realign:
	input:
		fq1="{run}/{sample}/regenerated_fastqs/{basename}.R1.fq.gz",
		fq2="{run}/{sample}/regenerated_fastqs/{basename}.R2.fq.gz",
		ref=config["ref_fa"]
	output:
		"{run}/{sample}/realigned_bams/{basename}.realigned.sam"
	shell:
		"bwa mem {input.ref} {input.fq1} {input.fq2} > {output}"


rule cram:
	input:
		sam="{run}/{sample}/realigned_bams/{basename}.realigned.sam",
		ref=config["ref_fa"]
	output:
		"{run}/{sample}/realigned_bams/{basename}.realigned.cram"
	shell:
		config["path_to_scramble"] + " -I sam -O cram -7 -V 3.0 -t 32 -B -r {input.ref} {input.sam} {output}"


rule flagstats:
	input:
		sam="{run}/{sample}/realigned_bams/{basename}.realigned.sam",
		cram="{run}/{sample}/realigned_bams/{basename}.realigned.cram"
	output:
		flagstat_sam="{run}/{sample}/realigned_bams/{basename}.flagstat_sam",
		flagstat_cram="{run}/{sample}/realigned_bams/{basename}.flagstat_cram"
	shell:
		"samtools flagstat {input.sam} > {output.flagstat_sam} && "
		"samtools flagstat {input.cram} > {output.flagstat_cram}"





