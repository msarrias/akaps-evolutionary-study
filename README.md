## Investigating the evolutionary history of the AKAP79 protein

### We will study the evolutionary history of the AKAP79 protein.

We start by looking at this protein in different species. For our analysis we selected the next model organism, each chosen to provide coverage 
 of the tree of life:


<img width="609" alt="image" src="https://user-images.githubusercontent.com/16377368/182667771-85023bb1-8c58-44dc-8a03-035febca9eef.png">


### AKAP79
We use as reference the AKAP79 sequence in human
- **Gene:** AKAP5
- **Protein name:** AKAP79, AKAP5 ,…
- **Binding partners:** PKC,PKA, CaM, PP2B,…

#### Questions we want to address:

- Is the AKAP79 protein architecture the same in all species?
- How can we find the binding partners in the sequences? 
    - As a first step, by assesing the conservation of the msa, i.e., looking that the aligned sequences match with the reference sequence region. 
- Can we find repeated instances?
    - We can use a profile hmm for identifying binding regions and possible repetitions.
    
Note:
- One cannot conclude that a binding region is present in a species solely by looking at a msa
- the profile hmm might fail to identify binding regions, specially when few seed sequences are used for creating the model.


### key points:
- All AKAP79 protein sequences in the selected model species are  human orthologs, with fish being the only exception.
- Since AKAP79 is known to be important for human synapses, and fish has synapses, one would expect to find it also in fish.
- What could this mean?
    * That some functions of AKAP79 are not important for fish and they can survive without it.
    * That there is another protein performing the same function as AKAP79 in humans.
    
For that reason, we searched for the AKAP79 protein in the genome of 80 fish species (85 including subspecies) on the Ensembl database.
- Our findings show that the AKAP79 protein did not appear in the genome of 80 fish species (including jaw-less, cartilaginous and bony fish), 
which could indicate that AKAP79 was not part of first vertebrates.

## AKAP79 architecture

<img width="463" alt="image" src="https://user-images.githubusercontent.com/16377368/182708369-df5f62af-62c1-47c0-93eb-5161715cd11e.png">


## [WSK domain organisation](http://pfam.xfam.org/family/PF03832#tabview=tab1)
As for the WSK motif, the Pfam database shows 4 different architectures (including more than 1 sequence) with different repeat patterns.
As can be seen in the image below these patterns are only found in the AKAP79 and Gravin proteins. 
Which prompts the question: what is the relationship between AKAP79 and Gravin? - Can we find Gravin in fish species?

<img width="1087" alt="image" src="https://user-images.githubusercontent.com/16377368/185901669-980263bc-f594-4d87-a91d-7d7164e9a8e5.png">

### key points:
- Unlike AKAP79, Gravin (AKAP12) appeared in almost all fish species 
- Similarites between AKAP79 and Gravin are low. This could indicate that AKAP79 did not originate from Gravin. 
- One can say that similarites between AKAP79 and Gravin are restricted to the WSK motif.


### Gravin

<img width="751" alt="image" src="https://user-images.githubusercontent.com/16377368/182708630-99115926-e8c0-453a-9f0d-b643e6948646.png">

* We use as reference the AKAP12 sequence in human
- **Gene:** AKAP12
- **Protein name:** AKAP12, Gravin,…
- **Binding partners:** PKC,PKA, CaM,…




