#! /usr/bin/python3


from __future__ import annotations
import click
import subprocess
from typing import List, Dict, NamedTuple
from collections import namedtuple
from dataclasses import dataclass


@dataclass
class Exon:
    start: int
    end: int

    def __lt__(self, other: Exon):
        """Make Exon sortable"""
        return(self.start < other.start)

    def __len__(self):
        return(abs(self.end - self.start))

    def reverse(self) -> Exon:
        return(Exon(start=self.end, end=self.start+1))


class CCDS:
    def __init__(self, chrom: str, gene: str, ccds_id: str, strand: str, exons: str):
        self.chrom: str = chrom
        self.gene: str = gene
        self.ccds_id: str = ccds_id
        self.strand: str = strand
        # Saves a function call. Not evaluated until you access exons field
        self.exons: List[Exon] = self._parse_cds_locations(exons)

    def _parse_cds_locations(self, cds_str: str) -> List[Exon]:
        exl = []
        for exon in cds_str.strip('[]').split(", "):
            for ex in exon.strip().split('='):
                e = ex.split('-')
                exl.append(Exon(start=int(e[0]), end=int(e[1])))
        if self.strand == "+":
            return(sorted(exl))
        elif self.strand == "-":
            m = map(lambda x: x.reverse(), reversed(exl))
            return(sorted(m))
        else:
            raise ValueError("strand attribute must be either '+' or '-'")


class Gene:
    def __init__(self, name: str):
        self.name = name
        self.ccdss = {}

    def add_ccds(self, chrom: str, gene: str, ccds_id: str, strand: str, exons: str):
        self.ccdss[ccds_id] = CCDS(chrom, gene, ccds_id, strand, exons)

    def union(self) -> List[Exon]:

        if list(self.ccdss.values()):
            strand = list(self.ccdss.values())[0].strand
        else:
            raise ValueError(f"No CDS for this {self.name} gene")

        un = {}

        for ex in [ex for cds in self.ccdss.values() for ex in cds.exons]:
            if ex.start in un:
                if len(ex) < len(un[ex.start]):
                    un[ex.start] = ex
            else:
                un[ex.start] = ex

        if strand == "+":
            return(sorted(un.values()))
        elif strand == "-":
            return(list(reversed(sorted(un.values()))))
        else:
            raise ValueError(f"Strand for gene {self.name} must be either '+' or '-'")


@click.command()
@click.option('--ccds', help='CCDS current release for GRCh38', type=str, default="/home/max/Work/Bioinformatics/refs/CCDS.current.txt")
@click.option('--genes', help='Gene list file', type=str, default="/home/max/Work/Bioinformatics/refs/Human_RTK_names")
@click.option('--ref', help='Genom reference indexed fasta', type=str, default="/home/max/Work/Bioinformatics/refs/GRCh38.p13.genome.fa")
def main(ccds, genes, ref):


    with open(f"{genes}", "r") as gl:
        gene_list = [Gene(line.strip()) for line in gl]

    #gene = Gene("EGFR")

    with open(f"{ccds}", "r") as cdsf:
        for line in cdsf:
            l = line.strip().split("\t")
            for gene in gene_list:
                if l[2] == gene.name and l[5] == "Public":
                    gene.add_ccds(l[0], l[2], l[4], l[6], l[9])

    for gene in gene_list:
        print(gene.union())


if __name__ == '__main__':
    main()
