IDS, = glob_wildcards("isolates/raw_reads/{id}_R1.fastq.gz")
configfile: "config.yaml"

rule all:
	input:
		expand("kraken_out/{sample}_kraken2_report.txt", sample = IDS),
		expand("quast_out/{sample}", sample = IDS),
		expand("abricate_out/vfdb/{sample}_vfdb.tsv", sample = IDS),
		expand("amrfinder_out/{sample}.tsv", sample = IDS),
		expand("mlst/{sample}.tsv", sample = IDS),
		expand("prokka_out/{sample}", sample = IDS),
                expand("kraken_unmapped_out/{sample}_kraken2_report.txt", sample = IDS),
		"coverage.tsv",
		"summary.tsv"

rule fastp:
	input:
		fw = "isolates/raw_reads/{sample}_R1.fastq.gz",
		rv = "isolates/raw_reads/{sample}_R2.fastq.gz"
	output:
		fw = "isolates/trimmed_reads/{sample}_R1.fastq.gz",
		rv = "isolates/trimmed_reads/{sample}_R2.fastq.gz",
		html = "fastp_out/html/{sample}_fastp.html",
		json = "fastp_out/json/{sample}_fastp.json"
	conda:
		"envs/fastp.yaml"
	params:
		compression_level = config["fastp"]["compression_level"],
		general = config["fastp"]["general"]
	log:
		"logs/fastp/{sample}.log"
	threads: 8
	shell:
		"""
		fastp -w {threads} -z {params.compression_level} -i {input.fw} -o {output.fw} -I {input.rv} -O {output.rv} {params.general} --html {output.html} --json {output.json} 2>&1>{log}
		"""

rule kraken2:
	input:
		fw = "isolates/trimmed_reads/{sample}_R1.fastq.gz",
		rv = "isolates/trimmed_reads/{sample}_R2.fastq.gz"
	output:
		report = "kraken_out/{sample}_kraken2_report.txt"
	conda:
		"envs/kraken.yaml"
	params:
		general = config["kraken"]["general"],
		db = config["kraken"]["db"]
	log:
		"logs/kraken2/{sample}.log"
	threads: 8
	shell:
		"""
		kraken2 --db {params.db} {params.general} --threads {threads} --report {output.report} {input.fw} {input.rv} 2>&1>{log}
		"""

rule shovill:
	input:
		fw = "isolates/trimmed_reads/{sample}_R1.fastq.gz",
		rv = "isolates/trimmed_reads/{sample}_R2.fastq.gz"
	output:
		assembly = "genomes/{sample}.fasta",
		shovill = directory("shovill_out/{sample}")
	params:
		minlen = config["shovill"]["minlen"],
		ram = config["shovill"]["ram"],
		depth = config["shovill"]["depth"],
		assembler = config["shovill"]["assembler"],
		tmpdir = config["shovill"]["tmpdir"]
	conda:
		"envs/shovill.yaml"
	log:
		"logs/shovill/{sample}.log"
	threads: 16
	shell:
		"""
		shovill --assembler {params.assembler} --outdir {output.shovill} --tmp {params.tmpdir} --depth {params.depth} --cpus {threads} --ram {params.ram} --minlen {params.minlen} --R1 {input.fw} --R2 {input.rv} 2>&1>{log}
		cp {output.shovill}/contigs.fa {output.assembly}
		"""

rule quast:
	input:
		"genomes/{sample}.fasta"
	output:
		directory("quast_out/{sample}")
	conda:
		"envs/quast.yaml"
	log:
		"logs/quast/{sample}.log"
	threads: 8
	shell:
		"""
		quast --threads {threads} -o {output} {input} 2>&1>{log}
		"""

rule abricate:
	input:
		assembly = "genomes/{sample}.fasta"
	output:	
		vfdb = "abricate_out/vfdb/{sample}_vfdb.tsv"
	conda:
		"envs/abricate.yaml"
	params:
		minid = config["abricate"]["minid"],
		mincov = config["abricate"]["mincov"],
		vfdb = config["abricate"]["vfdb"]
	log:
		"logs/abricate/{sample}.log"
	threads: 8
	shell:
		"""
		abricate --db {params.vfdb} --nopath --minid {params.minid} --mincov {params.mincov} {input.assembly} > {output.vfdb} 2>>{log}
		"""

rule amrfinder:
	input:
		"genomes/{sample}.fasta"
	output:
		"amrfinder_out/{sample}.tsv"
	conda:
		"envs/amrfinder.yaml"
	log:
		"logs/amrfinder/{sample}.log"
	threads: 16
	shell:
		"""
		amrfinder -u
		amrfinder --threads {threads} --nucleotide {input} --output {output} 2>&1>{log}
		"""

rule mlst:
	input:
		"genomes/{sample}.fasta"
	output:
		"mlst/{sample}.tsv"
	conda:
		"envs/mlst.yaml"
	log:
		"logs/mlst/{sample}.log"
	shell:
		"""
		mlst {input} > {output} 2>>{log}
		"""

rule prokka:
	input:
		"genomes/{sample}.fasta"
	output:
		directory("prokka_out/{sample}")
	conda:
		"envs/prokka.yaml"
	params:
		kingdom = config["prokka"]["kingdom"],
		prefix = "{sample}"
	log:
		"logs/prokka/{sample}.log"
	threads: 8
	shell:
		"""
		prokka --force --outdir {output} --kingdom {params.kingdom} --cpus {threads} --prefix {params.prefix} {input} 2>&1>{log}
		if [ -f {output}/*.gff ]; then echo "{output} exists"; else exit 1; fi
		"""

rule snippy:
	input:
		reads = "metagenomics/{sample}_clean.fq.gz",
		ref = "genomes/{sample}.fasta"
	output:
		directory("snippy_out/{sample}")
	conda:
		"envs/snippy.yaml"
	params:
		general = config["snippy"]["general"]
	log:
		"logs/snippy/{sample}.log"
	threads: 8
	shell:
		"""
		snippy {params.general} --cpus {threads} --outdir {output} --ref {input.ref} --peil {input.reads} 2>{log}
		"""

rule kraken2_unmapped:
	input:
		snippy = "snippy_out/{sample}"
	output:
		report = "kraken_unmapped_out/{sample}_kraken2_report.txt"
	conda:
		"envs/kraken.yaml"
	params:
		general = config["kraken"]["general"],
		db = config["kraken"]["db"]
	log:
		"logs/kraken2_unmapped/{sample}.log"
	threads: 8
	shell:
		"""
		kraken2 --db {params.db} {params.general} --threads {threads} --report {output.report} {input.snippy}/snps.unmapped_R1.fq.gz {input.snippy}/snps.unmapped_R2.fq.gz 2>&1>{log}
		"""

rule bamcov:
	input:
		"snippy_out/{sample}"
	output:
		"bamcov_out/{sample}.tsv"
	log:
		"logs/bamcov/{sample}.log"
	shell:
		"""
		bamcov {input}/snps.bam > {output} 2>{log}
		"""

rule print_coverage:
	input:
		expand("bamcov_out/{sample}.tsv", sample=IDS)
	output:
		"coverage.tsv"
	log:
		"logs/print_coverage.log"
	shell:
		"""
		bash scripts/print_coverage.sh {input} > {output} 2>{log}
		"""

rule summary:
	input:
		"coverage.tsv",
		expand("mlst/{sample}.tsv", sample = IDS),
		expand("quast_out/{sample}", sample = IDS),
		expand("kraken_out/{sample}_kraken2_report.txt", sample = IDS)
	output:
		"summary.tsv"
	log:
		"logs/summary.log"
	shell:
		"""
		bash scripts/summary.sh > {output} 2>{log}
		"""
