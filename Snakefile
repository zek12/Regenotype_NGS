configfile: "config.yaml"

from os.path import *
import glob


FILES = glob.glob(config["input_dir"] + "DRG_037/**/DRG_037-F01_LIBx1_FINAL_b37_PE.aln.RG.dupes_marked.bam", recursive = True)

d = {}
for f in FILES:
	d[basename(f)[:-4]] = f

rule all:
	input:
		[
			# config["dest_dir_FASTQs"] + s + ".R1.fq.gz" for s in d.keys()
			config["dest_dir_BAMs"] + s + ".flagstat_cram" for s in d.keys()
		]



# rule sort_bam:
# 	input: config["input_dir"] + "{run}/{sample}/{basename}.bam"
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
	input:
		lambda wildcards: d[wildcards.sample]
	output:
		fq1 = config["dest_dir_FASTQs"] + "{sample}.R1.fq.gz",
		fq2 = config["dest_dir_FASTQs"] + "{sample}.R2.fq.gz",
		fq0 = config["dest_dir_FASTQs"] + "{sample}.unpaired.fq.gz"
	shell: "java -Xmx2g -jar " + config["path_picard"] + " SamToFastq I={input} FASTQ={output.fq1} SECOND_END_FASTQ={output.fq2} UNPAIRED_FASTQ={output.fq0}"


rule realign:
	input:
		fq1 = config["dest_dir_FASTQs"] + "{sample}.R1.fq.gz",
		fq2 = config["dest_dir_FASTQs"] + "{sample}.R2.fq.gz",
		ref = config["ref_fa"]
	output:
		config["dest_dir_BAMs"] + "{sample}.sam"
	shell:
		"bwa mem {input.ref} {input.fq1} {input.fq2} > {output}"


rule cram:
	input:
		sam = config["dest_dir_BAMs"] + "{sample}.sam",
		ref = config["ref_fa"]
	output:
		config["dest_dir_BAMs"] + "{sample}.cram"
	shell:
		config["path_to_scramble"] + " -I sam -O cram -7 -V 3.0 -t 32 -B -r {input.ref} {input.sam} {output}"


rule flagstats:
	input:
		sam = config["dest_dir_BAMs"] + "{sample}.sam",
		cram = config["dest_dir_BAMs"] + "{sample}.cram"
	output:
		flagstat_sam = config["dest_dir_BAMs"] + "{sample}.flagstat_sam",
		flagstat_cram = config["dest_dir_BAMs"] + "{sample}.flagstat_cram"
	shell:
		"samtools flagstat {input.sam} > {output.flagstat_sam} && "
		"samtools flagstat {input.cram} > {output.flagstat_cram}"





