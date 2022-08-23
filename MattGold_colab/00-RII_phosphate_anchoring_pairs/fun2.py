from fun import *


def summary(feature, hmm_hits, hmm_obj):
    prot_with_hits = [key for key, value in hmm_hits.items() if value != {}]
    print(f'The {feature} was found in:')
    seqs = {key: hmm_obj.parse_fasta_file('fasta_files/' + key +'_oma.fa') for key in prot_with_hits}
    print([(prot, str(len(hmm_hits[prot])) + '/' + str(len(seqs[prot])))for prot in prot_with_hits])
    print('Missing in:')
    print([key for key, value in hmm_hits.items() if value == {}])
    return [prot for prot in prot_with_hits]


def align_ref_region(hmm_hits, prots, hmm_obj, feature):
    ref_dict = {prot: [(key, value) for key,
                       value in hmm_hits[prot].items() if 'HUMAN' in key][0] 
                for prot in prots}
    regions_dict = {}
    aln = {}
    missing_seqs = {}
    for prot_ID in prots:
        msa = hmm_obj.parse_fasta_file('fasta_files/msa/'+ prot_ID +'_oma_msa.fa')
        seqs = hmm_obj.parse_fasta_file('fasta_files/'+ prot_ID +'_oma.fa')
        base_string = ref_dict[prot_ID][1][1]
        ref = ref_dict[prot_ID][0].rsplit('-')[0]
        start, end = hmm_obj.find_binding_region(base_string, ref, msa)
        filename = 'fasta_files/binding_regions/'+ prot_ID + feature + '.fa'
        missing_seqs[prot_ID] = [id_ for id_ in seqs.keys() if id_ not in
           [i.rsplit('-')[0] for i in hmm_hits[prot_ID].keys()]]
        regions_dict[prot_ID], aln[prot_ID] = seq_domain_alignment(msa, seqs, start, end, filename,
                                                           binding_partner = feature)
    return missing_seqs, regions_dict, aln


def visualize_aln(aln, missing_seqs):
    for prot_ID, value in aln.items():
        print('--------')
        print(prot_ID)
        print('--------')
        alv.view(value)
        print('--------')
        print('species with no hits:', missing_seqs[prot_ID])
        print('--------')
        print()
        