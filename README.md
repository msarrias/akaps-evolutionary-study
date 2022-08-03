#### We will study the evolutionary history of the AKAP79 protein, first, by looking at this protein in different species.

For our analysis we selected the next model organism, each chosen to provide coverage of the tree of life:


<img width="444" alt="image" src="https://user-images.githubusercontent.com/16377368/182665357-7f4e8649-141c-4eb0-904e-5aa2f5e8f6dd.png">

### AKAP79
* We use as reference the AKAP79 sequence in human
- **Gene:** AKAP5
- **Protein name:** AKAP79, AKAP5 ,…
- **Binding partners:** PKC,PKA, CaM, PP2B,…

#### Questions we want to address:

- Is the AKAP79 protein architecture the same in all species?
- How can we find the binding partners in the sequences? 
    - As a first step, by assesing the conservation of the msa, i.e., looking that the aligned sequences match with the reference sequence region. 
- Can we find repeated instances?
    - We can use a profile hmm for identifying binding regions and possible repetitions.
    
Disadvantages:
- One cannot conclude that a binding region is present in a species solely by looking at a msa
- the profile hmm might fail to identify binding regions, specially when few seed sequences are used for creating the model.


### Results
- All AKAP79 protein sequences in the selected model species are  human orthologs, with fish being the only exception.

Since AKAP79 is known to be important for human synapses, and fish has synapses, one would expect to find it also in fish.
- What could this mean?
    * That some functions of AKAP79 are not important for fish and they can survive without it.
    * That there is another protein performing the same function as AKAP79 in humans.
    
For that reason, we searched for the AKAP79 protein in the genome of 80 fish species (85 including subspecies) on the Ensembl database.
Among the species, we included:  jaw-less, cartilaginous and bony fish.
Our findings show that the AKAP79 protein did not appear in any of the fish species before mentioned, which could indicate that AKAP79 
was not part of first vertebrates. Nevertheless, unlike AKAP79, Gravin (AKAP12) appeared in almost all fish species 
*Gasterosteus aculeatus*, which raises the questions of;
- Which protein serves the same function as AKAP79 in fish? and, 
- How is Gravin related to AKAP79? (considering that they are not homologous)

Note that, even if  Gravin sequences are homologous, it would be difficult to compare them with AKAP79 sequences as they have low similarities,
which could indicate that AKAP79 did not originate from Bravin. Hence, one can say that their similarities are restricted to the WSK motif.
In terms of the WSK motif architecture, the Pfam sequences database exhibited 4 different architectures (that include more than one sequence),
of which only one of them includes an interaction with another domain. Hence, besides showing different repetition patterns,
the WSK motif only interacts with the RII$\_\text{binding}\_1$ domain on Gravin protein sequences.

### Gravin

* We use as reference the AKAP12 sequence in human
- **Gene:** AKAP12
- **Protein name:** AKAP12, Gravin,…
- **Binding partners:** PKC,PKA, CaM,…




