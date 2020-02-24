#! /usr/bin/python3

import os
import click
import subprocess
import operator


@click.command()
@click.option('--ann', help='Gencode annotation gtf. Can be gzipped.', type=str, default="/home/max/Work/Bioinformatics/refs/gencode.v33lift37.annotation.gtf.gz")
@click.option('--genes', help='Gene list file', type=str, default="/home/max/Work/Bioinformatics/refs/Human_RTK_names")
def main(ann, genes):

    def longest_transcript(gene):
        transcripts = {}
        with open(f"kinases/{gene}.enst", "r") as ensfile:
            for line in ensfile:
                transcripts[line.strip()] = 0

        os.system(f"zgrep '\"{gene}\"' {ann} | awk '$3 == \"CDS\" {{print $0}}' > {gene}.gtf")

        for tr in transcripts:
            transcripts[tr] = int(subprocess.run([f"grep {tr} {gene}.gtf | awk '{{print $5 - $4}}' | paste -sd+ - | bc"], stdout=subprocess.PIPE, shell = True, text = True).stdout)

        os.system(f"rm {gene}.gtf")

        return(max(transcripts.items(), key=operator.itemgetter(1))[0])

    os.nice(20)

    kinases = []

    with open(f"{genes}", "r") as kfile:
        for line in kfile:
            kinases.append(line.strip())


    for gene in kinases:
        os.system(f"mkfifo {gene}")

    if not os.path.exists("./kinases"):
        os.makedirs("./kinases")


    start = f"zcat {ann} | awk '$3 == \"CDS\" {{print $0}}' | "
    tee = "tee " + " ".join(kinases) + " > /dev/null & "
    process = " & ".join([f"grep '\"{gene}\"' {gene} | awk '{{print $12}}' | sed -E 's/\"(ENST[0-9]*\.[0-9]*)_[0-9]*\";/\\1/g' | uniq > kinases/{gene}.enst" for gene in kinases])

    cmd = start + tee + process

    os.system(cmd)

    for gene in kinases:
        os.system(f"rm {gene}")

    with open("Human_RTK_representative_transcripts_CDS", "w") as out:
        with open(f"{genes}", "r") as inp:
            out.write("\n".join([f"{rtk.strip()}\t{longest_transcript(rtk.strip())}" for rtk in inp]))


if __name__ == '__main__':
    main()
