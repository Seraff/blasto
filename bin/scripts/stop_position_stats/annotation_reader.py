from BCBio import GFF
from Bio import SeqIO
from copy import deepcopy

class PredictedGene():
    def __init__(self, seq_feature, sequence):
        self.seq_feature = contig
        self.sequence = sequence

    def seq(self):
        return self.seq_feature.extract(self.sequence)

    def n_right_nucleotides(self, number):
        new_feature = deepcopy(self.seq_feature)
        pass

    def n_left_nucleotides(self, number):
        pass



class AnnotationReader():
    def __init__(self, genome_path, annotation_path):
        with open(genome_path) as f:
            self.seq_dict = SeqIO.to_dict(SeqIO.parse(f, "fasta"))

        with open(annotation_path) as f:
            self.recs = list(GFF.parse(f, base_dict=seq_dict))

        self.genes = []

        for rec in self.recs:
            for feature in rec.features:
                self.genes.append(PredictedGene(feature, rec.seq))

if __name__ == '__main__':
    gff = '/Users/serafim/bio/blasto/data/new_assembly/results/annotation.gff'
    genome = '/Users/serafim/bio/blasto/data/new_assembly/p57_genome_ra_polished.fa'

    print('hello')
