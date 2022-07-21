import subprocess, re, alv
from Bio import SeqIO,AlignIO
from io import StringIO
from dna_features_viewer import GraphicFeature, GraphicRecord
import matplotlib.pyplot as plt

def seq_domain_alignment(msa_seqs,
                         seqs,
                         start,
                         end,
                         binding_partner,
                         binding_part_dict = {}):
    temp_msa = {}
    if not binding_part_dict:
        binding_part_dict = {specie:[] for specie in msa_seqs.keys()}
    max_len = max(map(len, msa_seqs.keys()))
    for specie, seq in msa_seqs.items():
        temp_species = specie + ' ' *  (max_len - len(specie))
        temp_msa[temp_species] = seq[start:end]
        binding_part_dict[specie].append((binding_partner,
                                          re.search(temp_msa[temp_species].replace('-',''),
                                                    seqs[specie]).span()))
        count = len(temp_msa)
        length = max(map(len, temp_msa.values()))
        msa = f" {count} {length}\n"
        msa += '\n'.join(f"{prot_id[0:10]} {sequence}" for prot_id,
                         sequence in temp_msa.items())
        aln = AlignIO.read(StringIO(msa), 'phylip')
    return binding_part_dict, aln


def visualiza_structure(host_guest_dict,
                        colors_,
                        regions_dict,
                        title,
                        seqs_dict,
                        figsize=(4,20)
                       ):
    fig, ax = plt.subplots(len(regions_dict), figsize=figsize)
    ax[0].set_title(title, loc='center', weight = 'bold')
    for idx, (specie_id, regions_coord) in enumerate(regions_dict.items()):
        name = specie_id.rsplit('_akap')[0]
        features = [GraphicFeature(start=region[1][0], end=region[1][1],
                                  color=colors_[region[0]],
                                   label=region[0]) for region in regions_coord if region[1] != (0,0)]
        record = GraphicRecord(sequence_length=len(seqs_dict[specie_id]), features=features)
        ax[idx].text(-10, 0, 
                     host_guest_dict[name].replace('_', ' ').capitalize(),
                     ha="right", va="center",
                     fontsize=11)
        record.plot(ax=ax[idx])
    
    
def colors():
    return {'WSK': '#ffcccc',
            'WAK':"#cffccc",
            'PIIT': "#ccccff", 
            'RII_binding_1':"#ffd700",
            'PKC':"#00FFFF",
            'PKA':"#069AF3"}


def akap5_species_dic():
    species_dna_ID = {'ENSG00000179841': 'Homo_sapiens',        
                  'ENSMUSG00000021057': 'Mus_musculus',     
                  'ENSSSCG00000002268': 'Sus_scrofa',       
                  'ENSBTAG00000009642': 'Bos_taurus',       
                  'ENSCAFG00000016026': 'Canis_lupus',       
                  'ENSMODG00000009460': 'Monodelphis_domestica',   
                  'ENSOANG00000014857': 'Ornithorhynchus_anatinus',  
                  'ENSGALG00000011721': 'Gallus_gallus',   
                  'XP_008110209': 'Anolis_carolinensis',           
                  'XP_030130993': 'Taeniopygia_guttata',     
                  'XP_004917258': 'Xenopus_tropicalis',    
                  }
    species_prot_ID = {'ENSP00000378207': 'Homo_sapiens',        
                  'ENSMUSP00000114495': 'Mus_musculus',     
                  'ENSSSCP00000002463': 'Sus_scrofa',       
                  'ENSBTAP00000012705': 'Bos_taurus',       
                  'ENSCAFP00000023607': 'Canis_lupus',       
                  'ENSMODP00000011817': 'Monodelphis_domestica',   
                  'ENSOANP00000023393': 'Ornithorhynchus_anatinus',  
                  'ENSGALP00000019127': 'Gallus_gallus',   
                  'ENSTRUG00000018243' : 'Takifugu_rubripes_akap12b',
                  'XP_008110209': 'Anolis_carolinensis',           
                  'XP_030130993': 'Taeniopygia_guttata',     
                  'XP_004917258': 'Xenopus_tropicalis'
                  }

    species_guest_to_host_map_dic = {'human': 'Homo_sapiens',
                                 'mouse': 'Mus_musculus',
                                 'pig': 'Sus_scrofa',
                                 'cow': 'Bos_taurus',
                                 'dog': 'Canis_lupus',
                                 'opossum': 'Monodelphis_domestica',
                                 'platypus': 'Ornithorhynchus_anatinus',
                                 'chicken': 'Gallus_gallus',
                                 'anole':  'Anolis_carolinensis',
                                 'zebra_finch': 'Taeniopygia_guttata',
                                 'zebra_fish' : 'Danio_rerio',
                                 'clawed_frog': 'Xenopus_tropicalis',
                                 'takifugu' : 'Takifugu_rubripes'}
    species_host_to_guest_map_dic = {'Homo_sapiens':'human',
                                     'Mus_musculus': 'mouse',
                                     'Sus_scrofa':'pig',
                                     'Bos_taurus': 'cow',
                                     'Canis_lupus':'dog',
                                     'Monodelphis_domestica': 'opossum',
                                     'Ornithorhynchus_anatinus':'platypus',
                                     'Gallus_gallus':'chicken',
                                     'Anolis_carolinensis':'anole lizard',
                                     'Taeniopygia_guttata': 'zebra_finch',
                                     'Danio_rerio': 'zebra_fish',
                                     'Xenopus_tropicalis': 'clawed_frog',
                                     'Takifugu_rubripes': 'takifugu' }
    
    return(species_dna_ID, species_prot_ID, species_guest_to_host_map_dic,species_host_to_guest_map_dic)


def akap12_species_dic():
    species_prot_ID = {'ENSG00000131016': 'Homo_sapiens',        
                  'ENSMUSG00000038587': 'Mus_musculus',     
                  'ENSSSCG00000031733': 'Sus_scrofa',       
                  'ENSBTAP00000057364': 'Bos_taurus',       
                  'ENSCAFG00000000419': 'Canis_lupus',       
                  'ENSMODG00000046586': 'Monodelphis_domestica',   
                  'ENSOANG00000039295': 'Ornithorhynchus_anatinus',  
                  'ENSGALG00000036240': 'Gallus_gallus',   
                  'ENSACAG00000006324': 'Anolis_carolinensis',           
                  'ENSTGUG00000029620': 'Taeniopygia_guttata', 
                  'ENSEBUG00000011356': 'Eptatretus_burgeri',
                  'ENSXETG00000037915': 'Xenopus_tropicalis',
                  'ENSTRUG00000018243' : 'Takifugu_rubripes',
                  'ENSCMIG00000006916': 'Callorhinchus_milii',
                  'ENSDARG00000091792' : 'Danio_rerio_akap12a',
                  'ENSDARG00000055678' : 'Danio_rerio_akap12b'
                  }

    species_guest_to_host_map_dic = {'human': 'Homo_sapiens',
                                 'mouse': 'Mus_musculus',
                                 'pig': 'Sus_scrofa',
                                 'cow': 'Bos_taurus',
                                 'dog': 'Canis_lupus',
                                 'opossum': 'Monodelphis_domestica',
                                 'platypus': 'Ornithorhynchus_anatinus',
                                 'chicken': 'Gallus_gallus',
                                 'anole':  'Anolis_carolinensis',
                                 'zebra_finch': 'Taeniopygia_guttata',
                                 'zebra_fish' : 'Danio_rerio',
                                 'inshore hagfish':'Eptatretus burgeri',
                                 'clawed_frog': 'Xenopus_tropicalis',
                                 'elephant_shark':'Callorhinchus_milii'}
    return(species_prot_ID, species_guest_to_host_map_dic)
