import subprocess, re, alv, os, warnings
import random
import numpy as np
from Bio.Blast import NCBIXML, NCBIWWW
from Bio import SeqIO,AlignIO
from io import StringIO
from dna_features_viewer import GraphicFeature, GraphicRecord
import matplotlib.pyplot as plt
from ete3 import PhyloTree, Tree, faces, TreeStyle


def seq_domain_alignment(msa_seqs,
                         seqs,
                         start,
                         end,
                         file_name,
                         motif_coords=[],
                         binding_partner = '',
                         binding_part_dict = {}):
    temp_msa = {}
    if not binding_part_dict:
        binding_part_dict = {specie:{} for specie in msa_seqs.keys()}
    max_len = max(map(len, msa_seqs.keys()))
    with open(file_name, "w") as ofile_dna:
        for specie, seq in msa_seqs.items():
            temp_species = specie + ' ' *  (max_len - len(specie))
            temp_msa[temp_species] = seq[start:end]
            temp_seq = temp_msa[temp_species].replace('-','')
            if motif_coords:
                temp_binding_partner = ''.join([seq[i] for i in motif_coords])
                if '-' in temp_binding_partner:
                    pass
                else:
                    binding_partner = temp_binding_partner
            if binding_partner not in binding_part_dict[specie]:
                binding_part_dict[specie][binding_partner] = []
            binding_part_dict[specie][binding_partner].append(re.search(temp_seq,
                                                                        seqs[specie]).span())
            ofile_dna.write(">" + specie.rsplit('_akap')[0] +  "\n" + temp_msa[temp_species] + "\n")
        count = len(temp_msa)
        length = max(map(len, temp_msa.values()))
        msa = f" {count} {length}\n"
        msa += '\n'.join(f"{prot_id[0:10]} {sequence}" for prot_id,
                         sequence in temp_msa.items())
        aln = AlignIO.read(StringIO(msa), 'phylip')
    return binding_part_dict, aln


def visualize_structure(host_guest_dict,
                        colors_,
                        regions_dict,
                        title,
                        seqs_dict,
                        figsize=(5,20),
                        save_figs = False
                       ):
    if not save_figs:
        fig, ax = plt.subplots(len(regions_dict), figsize=figsize)
        ax[0].set_title(title, loc='center', weight = 'bold')
        for idx, (specie_id, regions_coord) in enumerate(regions_dict.items()):
            name = specie_id.rsplit('_akap')[0]
            features = []
            for feature, coordinates in regions_coord.items():
                features += [GraphicFeature(start=coord[0], end=coord[1],
                                          color=colors_[feature],
                                           label=feature) for coord in coordinates if coord != (0,0)]

            record = GraphicRecord(sequence_length=len(seqs_dict[specie_id]), features=features)
            ax[idx].text(-10, 0, 
                         host_guest_dict[name].replace('_', ' ').capitalize(),
                         ha="right", va="center",
                         fontsize=11)
            record.plot(ax=ax[idx],figure_width=3, figure_height=1.3)
    else:
        for idx, (specie_id, regions_coord) in enumerate(regions_dict.items()):
            name = specie_id.rsplit('_akap')[0]
            features = [GraphicFeature(start=region[1][0], end=region[1][1],
                                      color=colors_[region[0]],
                                       label=region[0]) for region in regions_coord if region[1] != (0,0)]
            record = GraphicRecord(sequence_length=len(seqs_dict[specie_id]), features=features)

            record.plot(figure_width=4, figure_height=1.3) 
            plt.savefig('figs/' + title+ '-' + str(name), bbox_inches='tight') 
            

def layout_tol(node):
    img_path = 'phylo_figs/'
    host_guest_dict = species_host_to_guest_map_dic()
    faces_dict = {x.rsplit('.')[0] : faces.ImgFace(img_path+x) 
              for x in [i for i in os.listdir(img_path) if '.jpeg' in i]}
    if node.is_leaf():
        temp_node_name = host_guest_dict[node.name].replace('_', ' ').capitalize()
        node.add_face(faces.TextFace(temp_node_name,
                                     fsize = 20),
                      column=0)
        if node.name in faces_dict:
            descFace = faces_dict[node.name]
            descFace.border.margin = 1
            faces.add_face_to_node(descFace, node, column=1, aligned=True)

            
def layout_feat(node):
    img_path = 'figs/'
    host_guest_dict = species_host_to_guest_map_dic()
    faces_dict = {x.rsplit('-')[1].rsplit('.')[0] : faces.ImgFace(img_path+x,100,250 ) 
              for x in [i for i in os.listdir(img_path) if '.png' in i]}
    if node.is_leaf():
        temp_node_name = host_guest_dict[node.name].replace('_', ' ').capitalize()
        node.add_face(faces.TextFace(temp_node_name,
                                     fsize = 20),
                      column=0)
        if node.name in faces_dict:
            descFace = faces_dict[node.name]
            descFace.border.margin = 1
            faces.add_face_to_node(descFace, node, column=1, aligned=True)

            
def render_phylogeny(nw_tree_dir, layout):  
    t = Tree(nw_tree_dir)
    ts = TreeStyle()
    ts.layout_fn = layout
    ts.show_leaf_name = False
    ts.scale =  1
    ts.show_scale = False
    return ts, t


def render_msa_phylo(nwk_direct, align, host_guest_dict):
    t = PhyloTree(nwk_direct, format = 1, alignment=align, alg_format="fasta")
    for leaf in t.iter_leaves():
        leaf.name = host_guest_dict[leaf.name].capitalize().replace('_', ' ')
    return t


def read_fasta(filename):
    with open(filename, 'r') as f:
        fasta = f.read()
    return fasta


def read_blast_into_dict(blast_xml_dir, key='5'):
    blast_dict = {}
    with open(blast_xml_dir,"r") as blast_temp:
        records= NCBIXML.parse(blast_temp)
        i = 0
        for q, item in enumerate(records):
            for alignment in item.alignments:
                for hsp in alignment.hsps:
                    hit_dic = {}
                    temp_acc_list = [blast_dict[j]["accession"] for j in list(blast_dict.keys())]
                    temp_acc = (alignment.title.split("|"))[1].split("|")[0]
                    if temp_acc not in temp_acc_list:
                        name = (alignment.title.split("| "))[1].split(" [")[0]
                        if key in name:
                            hit_dic['accession'] = temp_acc
                            hit_dic['specie'] = (alignment.title.split("["))[1].split("]")[0]
                            hit_dic['name'] = name
                            hit_dic['sequence'] = alignment.title
                            hit_dic['length'] = alignment.length
                            hit_dic['e-value'] = hsp.expect
                            hit_dic['sbjct'] = hsp.sbjct.replace('-','') #GenBank sequence
                            blast_dict[i] = hit_dic
                            i += 1
    return blast_dict


def NCBI_order_ids():
    return {'rodent':'9989',
            'Artiodactyla':'91561',
            'Carnivora':'33554',
            'Didelphimorphia':'38605',
            'Monotremata':'9255',
            'Galliformes':'8976',
            'Squamata':'8509',
            'Passeriformes':'9126',
            'Cypriniformes':'7952',
            'Anura':'8342',
            'Tetraodontiformes':'31022'}


def regions_colors(regions_dict):
    regions = set(sum([list(species_dict.keys()) for species_dict in regions_dict.values()], []))
    binding_regions_colors = {'WSK': '#ffcccc',
                              'WAK':"#cffccc",
                              'PIIT': "#ccccff", 
                              'RII_binding_1':"#ffd700",
                              'PKC':"#00FFFF",
                              'PKA':"#069AF3"}
    for region in regions:
        if region not in binding_regions_colors:
            binding_regions_colors[region] = "#"+''.join([random.choice('0123456789ABCDEF') for i in range(6)])
    return binding_regions_colors


def species_host_to_guest_map_dic():
    return {'Homo_sapiens':'human',
            'Mus_musculus': 'mouse',
            'Sus_scrofa':'wild boar',
            'Bos_taurus': 'cow',
            'Canis_lupus':'dog',
            'Monodelphis_domestica': 'opossum',
            'Ornithorhynchus_anatinus':'platypus',
            'Gallus_gallus':'chicken',
            'Anolis_carolinensis':'anole lizard',
            'Taeniopygia_guttata': 'zebra_finch',
            'Danio_rerio': 'zebra_fish',
            'Xenopus_tropicalis': 'clawed_frog',
            'Takifugu_rubripes': 'takifugu',
            'Danio_rerio_akap12b': 'zebra_fish_b',
            'Danio_rerio_akap12a': 'zebra_fish_a',
            'Callorhinchus_milii': 'elephant_shark',
            'Eptatretus_burgeri': 'hagfish'}


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
                                 'wild_boar': 'Sus_scrofa',
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
    
    return(species_dna_ID, species_prot_ID, species_guest_to_host_map_dic)


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
