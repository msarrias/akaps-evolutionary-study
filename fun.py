import subprocess, re, alv
from Bio import SeqIO
from Bio import AlignIO
from io import StringIO

def seq_domain_alignment(msa_seqs, start, end):
    temp = {}
    max_len = max(map(len, msa_seqs.keys()))
    for specie, msa_seq in msa_seqs.items():
        temp_species = specie + ' ' *  (max_len - len(specie))
        temp[temp_species] = msa_seq[start:end]
    count = len(temp)
    length = max(map(len, temp.values()))
    msa = f" {count} {length}\n"
    msa += '\n'.join(f"{prot_id[0:10]} {sequence}" for prot_id,
                     sequence in temp.items())
    aln = AlignIO.read(StringIO(msa), 'phylip')
    alv.view(aln)