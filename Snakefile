configfile: "config.yaml"


rule regenerate_fastqs:
	input:
		"{sample}.bam"
	output:
		fq1="{sample}.R1.fq",
		fq2="{sample}.R2.fq"
	shell:
		"samtools fastq -1 {output.fq1} -2 {output.fq2} {input}"


rule realign:
	input:
		fq1="{sample}.R1.fq",
		fq2="{sample}.R2.fq",
		ref=config["ref_fa"]
		# ref=REF
	output:
		sam="{sample}.realigned.sam"
	shell:
		"bwa mem {input.ref} {input.fq1} {input.fq2} > {output.sam}"


rule cram:
	input:
		sam="{sample}.realigned.sam",
		ref=config["ref_fa"]
	output:
		"{sample}.realigned.cram"
	shell:
		config["path_to_scramble"] + " -I sam -O cram -7 -V 3.0 -t 32 -B -r {input.ref} {input.sam} {output}"


rule flagstats:
	input:
		sam="{sample}.realigned.sam",
		cram="{sample}.realigned.cram"
	output:
		flagstat_sam="{sample}.flagstat_sam.txt",
		flagstat_cram="{sample}.flagstat_cram.txt"
	shell:
		"samtools flagstat {input.sam} > {output.flagstat_sam} && "
		"samtools flagstat {input.cram} > {output.flagstat_cram}"




