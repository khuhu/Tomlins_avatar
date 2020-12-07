import pandas as pd
data = pd.read_csv("/home/kevhu/data/20180220mouseBams.txt")
SAMPLES = data['SampleName'].tolist()
#PATH = data['path'].tolist()

rule tvc_caller:
	input:
		bed="/srv/tvc-5.0.2-Ubuntu_12.04_x86_64-binary/ref/ptrim_bed/IAD124056_167_Designed.bed",
		ref="/srv/tvc-5.0.2-Ubuntu_12.04_x86_64-binary/ref/mm10/mm10.fasta",
		JSON="/srv/tvc-5.0.2-Ubuntu_12.04_x86_64-binary/ref/json/somatic_low_stringency_proton.json",
		bin="/srv/tvc-5.0.2-Ubuntu_12.04_x86_64-binary/bin/",
		trimBed="/srv/tvc-5.0.2-Ubuntu_12.04_x86_64-binary/ref/ptrim_bed/IAD124056_167_Designed.bed",
		hotspot="/mnt/DATA4/kevhu/choLab/vcfs/test/allMouseVarsCombined.hotspot.vcf",
		out="/mnt/DATA4/kevhu/choLab/vcfs/",
		bam=expand("test/{sampleName}.vcf.gz", sampleName=SAMPLE)
	output:
		"/mnt/DATA4/kevhu/choLab/vcfs/{sampleName}.vcf"
	shell:
		"""
		/mnt/DATA4/kevhu/scripts/copyOfTVC.py -b {input.bed} -i {input.bam}  -r {input.ref} -p {input.JSON} -B {input.bin}  -o {input.out} --primer-trim-bed {input.trimBed} -s {input.hotspot} --file-name {sampleName}
		"""

#rule rename_vcf:
#	input:
#		"/mnt/DATA4/kevhu/choLab/vcfs/TSVC_variants.vcf",
#		names=expand("{sampName}", sampName=SAMPLES)
#	output:
#		expand("/mnt/DATA4/kevhu/choLab/vcfs/{sample}.vcf", sample=SAMPLES)
#	shell:
#		"rename 's/TSVC_variants/{input.names}/ *.vcf*"

