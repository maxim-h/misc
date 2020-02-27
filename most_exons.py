#! /usr/bin/python3

import os
import click
import subprocess
import operator


@click.command()
@click.option('--ann', help='Gencode annotation gtf. Can be gzipped.', type=str, default="/home/max/Work/Bioinformatics/refs/gencode.v33lift37.annotation.gtf.gz")
@click.option('--genes', help='Gene list file', type=str, default="/home/max/Work/Bioinformatics/refs/Human_RTK_names")
@click.option('--ref', help='Genom reference indexed fasta', type=str, default="/home/max/Work/Bioinformatics/refs/GRCh37.primary_assembly.genome.fa")
def main(ann, genes, ref):

    def longest_transcript(gene):
        transcripts = []
        with open(f"genes/{gene}.enst", "r") as ensfile:
            for line in ensfile:
                transcripts.append(line.strip())

        #TODO: Super slow! parallelise with tee
        #os.system(f"zgrep '\"{gene}\"' {ann} | awk '$3 == \"exon\" {{print $0}}' > {gene}.gtf")

        for tr in transcripts:
            os.system(f"grep {tr} {gene}.gtf | gtf2bed > {tr}.bed")


        os.system("bedops --everything " + " ".join([f"{tr}.bed" for tr in transcripts]) + f" | awk '!a[$2*$3]++' > {gene}.bed")
        if not os.path.exists("./seq"):
            os.makedirs("./seq")
        os.system(f"bedtools getfasta -s -fi {ref} -bed {gene}.bed > ./seq/{gene}_exons.fa")

        for tr in transcripts:
            os.system(f"rm {tr}.bed")

        os.system(f"rm {gene}.bed {gene}.gtf")
        return(0)

    os.nice(20)
    gene_list = []

    with open(f"{genes}", "r") as kfile:
        for line in kfile:
            gene_list.append(line.strip())


    for gene in gene_list:
        os.system(f"mkfifo {gene}")

    if not os.path.exists("./genes"):
        os.makedirs("./genes")


    start = f"zcat {ann} | awk '$3 == \"exon\" {{print $0}}' | "
    tee = "tee " + " ".join(gene_list) + " > /dev/null & "
    process = " & ".join([f"grep '\"{gene}\"' {gene} | tee {gene}.gtf | awk '{{print $12}}' | sed -E 's/\"(ENST[0-9]*\.[0-9]*)_[0-9]*\";/\\1/g' | uniq > genes/{gene}.enst" for gene in gene_list])

    cmd = start + tee + process

    os.system(cmd)

    for gene in gene_list:
        os.system(f"rm {gene}")

    with open(f"{genes}", "r") as inp:
        res = [longest_transcript(rtk.strip()) for rtk in inp]

    if sum(res) == 0:
        print("Finished without errors")


if __name__ == '__main__':
    main()
