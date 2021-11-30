from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis as pt

fasta_sequence = SeqIO.parse(
    '/Users/morganbairamkulow/Downloads/genome_assemblies_genome_gb/ncbi-genomes-2021-11-30/GCF_004115585.1_ASM411558v1_genomic.gbff',
    'genbank')


class Bioin: # Путь до файла genebank
    def __init__(self, path_fasta_genbank):
        self.fasta_sequence = path_fasta_genbank

    def gene_id(self) -> list: 
        gene_res = []
        for fasta in self.fasta_sequence:  # Принимает genbank
            for feature in fasta.features:
                if feature.type == 'CDS':  # Выведем экзоны
                    for gene_id in feature.qualifiers['db_xref']:
                        print(f'Идентификатор гена: {gene_id}')
                        gene_res.append(gene_id)
        print(f'Список идентификаторов:\n{gene_res}')
        return gene_res

    def gene_location(self) -> dict:
        gene_loc = {}
        for fasta in self.fasta_sequence:
            for feature in fasta.features:
                if feature.type == 'CDS':
                    for gene_id in feature.qualifiers['db_xref']:
                        gene_loc[gene_id] = feature.location

        print(gene_loc)
        return gene_loc

    def calculate_len_exons(self) -> int:
        exon_count = 0
        for fasta in self.fasta_sequence:
            for feature in fasta.features:
                if feature.type == 'CDS':
                    for gene_id in feature.qualifiers['db_xref']:
                        exon_count += len(feature.extract(fasta.seq))
                        print(f'Ген:{gene_id}\tДлина: {len(feature.extract(fasta.seq))}')
        print(f'Суммарная длина экзонов {exon_count}')
        return exon_count

    def calculate_len_seq(self) -> int:
        sum_seq = 0
        for k in self.fasta_sequence:
            sum_seq += len(k.seq)
            print(k.name, len(k.seq))
        print(sum_seq)
        return sum_seq

    # def calculate_percent(self):
    #   ex = calculate_len_exons(self)
    #  seq_len = calculate_len_seq(self)
    # print(ex)
    # print(seq_len)
    # print(f'Длина последовательности:{seq_len}\tКоличесво экзонов:{ex}\tЧисло интронов: {seq_len - ex}\tОтношение экзон /интрон:{ex/(seq_len-ex)}')

    def amino_attitude(self):
        for fasta in self.fasta_sequence:
            for feature in fasta.features:
                if feature.type == 'CDS':
                    try:
                        for prot in feature.qualifiers['translation']:
                            pr = pt(prot)
                            print(pr.get_amino_acids_percent())

                    except:
                        print('Аминокислотная последовательность отсутсвует')

    def secondary_structure(self):
        for fasta in self.fasta_sequence:
            for feature in fasta.features:
                if feature.type == 'CDS':
                    try:
                        print('Спираль, поворот', 'b-стурктура')
                        for prot in feature.qualifiers['translation']:
                            pr = pt(prot)
                            print(pr.secondary_structure_fraction())
                    except:
                        print('Аминокислотная последовательность отсутсвует')


bb = Bioin(fasta_sequence)
