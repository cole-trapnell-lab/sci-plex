# Sci-Plex


[Sci-Plex] (https://science.sciencemag.org/content/early/2019/12/04/science.aax6234
) is a labeling sampling that enables sample labeling within sci-RNA-seq. This github repository documents processing pipelines, primary and secondary analyses, and other files used in the generation, analysis and interpretation of the data. 

### Contents
##### Primary Analysis and data
Our data derive primarily from 4 sequencing experiments. Which are all available on GEO both in their raw and processed forms.

| Name        | Experiment           |GEO Accession  | Analysis Scripts |
| :-------------: |:-----------:| :----:| :---:|
| sciPlex1      | Barnyard  Experiment| [GSM4150376](sciPlex1_HEK293T_NIH3T3)| barnyard_experiment|
| sciPlex2| Small Screen      |  [GSM4150377](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4150377) | small_screen|
| sciPlex3 | Large Screen      | [GSM4150378](sciPlex3_A549_MCF7_K562_screen) | large_screen |
| sciPlex4 | HDACi Phenocopy Rescue      | [GSM4150379](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4150379)| large_screen |

##### Processing the raw data 
The scripts and code used to process the data can be accessed within the [process\_from_raw](https://github.com/cole-trapnell-lab/sci-plex/tree/master/process_from_raw) folder. This code was used to process the data from a demulitplexed sequencing run to a gene expression matrix. 

##### Ordering Reagents for Hashing

As highlighted by the paper, oligo-hashing is an simple procedure employing off the shelf reagents. To make this more accessible we have included  96-well format IDT order forms for the hash-oligo sequences used in this study. These excel files can be found in the [reagents] (https://github.com/cole-trapnell-lab/sci-plex/tree/master/reagents/IDT_hash_oligo_order_forms) folder in this github repository. A single plate of 40nmols of hash-oligos from IDT should be sufficient for up to 800 labeling reactions per well. 

We are currently working on a protocols.io page for both [sci-RNA-seq](https://science.sciencemag.org/content/357/6352/661) and [sci-plex](https://science.sciencemag.org/content/early/2019/12/04/science.aax6234). As soon as this is ready we will update with a link here and advertise it more broadly. 
