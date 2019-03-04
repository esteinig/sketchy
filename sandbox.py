from sketchy.sketchy import Sketchy

sketchy = Sketchy()

sketchy.sort_fastq(
    file='st243_r9.4_reads.fq',
    fastq='st243_r9.4_reads_sorted.fq',
    shuffle=False
)

sketchy.predict_nanopore(
    sketch='saureus-v1.msh',
    fastq='st243_r9.4_reads_sorted.fq',
    tmp='./tmp',
    cores=16,
    extension='.fq',
    header=True
)


sketchy.predict_assemblies(
    assemblies='test/*.fasta',
    sketch='saureus-v1.mah',
    glob='*.fasta',
)
