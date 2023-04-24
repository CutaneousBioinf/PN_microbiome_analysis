system(echo 'Processing quality control on MOTHUR')

make.file(inputdir=., type=gz, prefix=galderma)

system(echo 'Reducing sequencing and PCR errors')

make.contigs(file=galderma.files)
screen.seqs(fasta=current, count=current, maxambig=0, maxhomop=8, minlength=450)
unique.seqs(fasta=current, count=current)
summary.seqs(fasta=current, count=current)

system(echo 'Processing improved sequence')

align.seqs(fasta=current, reference=silva.v1v3.align)
summary.seqs(fasta=current, count=current)
screen.seqs(fasta=current, count=current, start=1, end=12818, minlength=450)
filter.seqs(fasta=current, vertical=T, trump=.)
unique.seqs(fasta=current, count=current)
summary.seqs(fasta=current, count=current)

pre.cluster(fasta=current, count=current, diffs=5)
summary.seqs(fasta=current, count=current)

chimera.vsearch(fasta=current, count=current, dereplicate=t)
remove.seqs(fasta=current, accnos=current)
summary.seqs(fasta=current, count=current)

classify.seqs(fasta=current, count=current, reference=trainset18_062020.pds.fasta, taxonomy=trainset18_062020.pds.tax)
remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)
summary.seqs(fasta=current, count=current)

rename.file(fasta=current, count=current, taxonomy=current, prefix=galderma)

quit()