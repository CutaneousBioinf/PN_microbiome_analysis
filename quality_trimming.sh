pre="conda activate fastqc-0.11.9"
eval pre$

files_tab=$(ls $1|grep "fastq.gz$"|cut -d"_" -f1,2|uniq)

for i in $files_tab;do
	f1=$folder$i"_R1_001.fastq.gz"
	f2=$folder$i"_R2_001.fastq.gz"
	CMD="trim_galore -q 25 --phred33 --stringency 3 -o $2 --paired $f1 $f2 --fastqc_args \"--outdir $3\""
	eval $CMD	
done

qc_report="multiqc $3"
eval qc_report

files=$(ls $2|grep "fq.gz$"|cut -d"_" -f1,2|uniq)

for i in $files_tab;do
	in1=$2$i"_R1_001_val_1.fq.gz"
	out1=$2$i"_R1_001.fq.gz"
	CMD="mv $in1 $out1"
	eval $CMD
done

for i in $files_tab;do
	in2=$2$i"_R2_001_val_2.fq.gz"
	out2=$2$i"_R2_001.fq.gz"
	CMD="mv $in2 $out2"
	eval $CMD
done