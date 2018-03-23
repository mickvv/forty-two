*Forty-Two*, a next-generation tool for enriching alignments of orthologous genes
and gene families.

## Abstract

Problématique ~ reformuler premier paragraphe de l'intro.

To this end, we have developed a novel tool, "42", along the lines of *HaMStR*
[Ebersberger et al. (2009) *BMC Evol Biol* 9:157] and the newly released
*Orthograph* [M. Petersen et al. (2017), *BMC Bioinformatics], whose aim is to
add (and optionally align) sequences to thousands of preexisting multiple
sequence alignments (MSA) while controlling for orthology relationships and
potentially contaminating sequences. "42" uses advanced heuristics based on a
multiple Best Reciprocal Hit (multi-BRH) strategy against reference proteomes to
distinguish orthologous from paralogous sequences among homologues. It is fully
functional and has already been used in several phylogenomic manuscripts dealing
with various parts of the tree of life.

Precision-recall analyses of eukaryotic orthologous groups enrichment with
complete proteome and genome data have shown its heuristics to perform as
expected. By force of circumstances we designed comparison
tests with *Orthograph*. As long as we deal with transcriptomic data, both
sofwares play on an even field.  Where Forty-Two denotes is when genomic data is
used; it seems that its ability to deal with long sequence by trimming the
coding sub-sequence assessed as ortholog allows Forty-Two to detect sequences
from long genomic sequences.

Forty-Two has proven to be an all terrain tool, able to handle heterogenous data
types while performing as expected. Forty-Two can be directly installed from
CPAN. It is provided with configuration helper tools that makes it easy to use.
Nonetheless the experimented user will be satisfied with the generous set of
fine-tunable parameters. The current version is freely available at:
https://bitbucket.org/dbaurain/42/.

## Introduction

Identifying orthology relationships among sequences is fundamental in
phylogenomics; indeed, those are essential to understand evolution, diversity of
life and ancestry among organisms. To build alignments of orthologous sequences,
phylogenomic pipelines often start with a step of all-vs-all similarity search
followed by a clustering with an algorithm such as *OrthoFinder* [Emms and Kelly
(2015) *Genome Biol* 16:157]. For it to be as accurate as possible, proteomes of
good quality are needed. Since their availability is limited to a small subset
of the living beings, large-scale taxonomic phylogenomic analyses imply the
enrichment of preexisting orthologous groups with transcriptomic or genomic
data, hence the need for robust tools for identifying orthologues from
heterogeneous sequence data. However, obtaining good quality transcriptomic data
is far from trivial, that's why Forty-Two is provided with contamination
detection features that can prevent enrichment of orthologous MSAs with
contaminating sequences.

Besides, Forty-Two can be used as a stand alone contamination detection tool ,
hence without MSAs enrichment, in order to spot contaminants in transcriptomic
or genomic data. To this end Forty-Two comes along with a set of ribosomal MSAs
manually currated, maintained and continuously enriched with new species by
Hervé Philippe. Indeed, it is crucial that the taxonomic richness within these
MSAs is maximised for an ever more precise detection.

Here, we present the principles and algorithms underlying Forty-Two as well as
the results of an extensive test suite of its features. We designed a series of
tests in an ideal case scenario to prove Forty-Two heuristics that we further
transposed to a "real-world" case scenario using a home-made dataset of 480
orthologous groups (OGs). Then, we compared performances with a software,
Orthograph, with both transcriptomic and genomic data. In addition, we tested
its contamination-proof enrichment features with a dataset of OGs built for
Hymenoptera. Finally, using a dataset of transcriptomes from MMESTP, we
highlighted its stand-alone contamination detection feature.

## Methods
### F1-score

The F1 score is an average of precision and recall that does not take into 
account true negatives. Indeed, our choice is based on the fact that is it 
complex to define what true negatives are. It could be all the sequences 
constituting a genome but it would tend to overestimate specificity (TN/TN+FP) 
and accuracy (TP+TN/TP+TN+FP+FN) because of the huge difference in numbers with 
TP, FP and FN. An alternative would be to limit true negatives to the pool of 
homologous but being dependent on *BLAST* parameters it is likely subjective.

### Implementation
The Forty-Two software package is available on CPAN and is divided into three
main tools that handles (i) configuration setup (yaml-generator-42.pl), (ii)
MSAs enrichment by orthology assessment based on a multi-BRH approach
(Forty-Two) and (iii) debriefing tools including summarizers (debrief-42.pl and
debrief-taxR-42.pl) and enrichment control (one-on-one.pl). While the
user-friendly assistant configurator helps you reproduce previous runs alike,
*TODO???* __it can also spare computation time by re-using previous runs data__.

*Forty-Two* requires only one pre-processing step that is formatting
databases i.e. genome/transcriptome set to mine for candidate orthologs and the
reference proteome set used to assess orthology. Formatting implies two step:
(i) renaming deflines as what we call "must ids" (e.g. "genus species") followed
by (ii) turning files into blast databases with the makeblastdb command. All is
left is to follow the instructions given by the assistant configurator
(yaml-generator.pl).

*Forty-Two*'s heuristic is based on a multiple best reciprocal hit criterion; it
can be divided into three main steps:
1. Best-hit checklist building
    - Orthologous queries (sequences from the MSA to enrich, user defined) are
      used to perform a BLAST against the reference proteome. Each best-hit(s)
      in each proteome for each query are collected thus constituting a
      checklist of orthologous sequences in the reference proteome set. The idea
      behind this checklist is that any other sequence that is ortholgous to the
      concerned gene should point to the same best-hit(s) in each and every
      reference proteome.
    - The checklist is further consolidated by performing a BLAST of each
      best-hit against the remaining reference proteomes. If a collected
      best-hit does not point to the same best-hit for each and every reference
      proteomes it is discarded from the list. This is to robustly ensure
      orthology between best-hits.
2. Candidate orthologs mining
    - The same orthologous queries are used to mine the candidate
      genomes/transcriptomes for homologs (candidates orthologs) with BLAST.
3. Orthology test
    - Each homolog is used to perform a BLASTX versus each proteome of the
      reference proteome set and their best-hit are collected. If a candidate
      orhtolog point to the same best-hit as in the checklist for each and every
      reference proteome, it is then considered as an ortholog.

Figure 1. Forty-Two's heuristic illustrated. (1) Alignement to enrich. (a) Query
extraction step. (b) Blastp of queries VS a reference proteome set (2) and
collection of best-hits in order to build (c) a best-hit checklist. (d) tblastn
of queries VS genome set (3) to mine for orthologs in order to build (e) a list
of homologs (candidate orthologs). (f) blastx of homologs VS the same reference
proteome set and collection of best-hits. (g) Orthology affiliation by
confirming that a homologs’ best-hits list is strictly included in the best-hit
checklist (c).


# Results and discussion
## Proving heuristics
In order to prove *Forty-Two*'s heuristics, we desgined an ideal test case
scenario. Starting from our own OGs set that we obtained by running OrthoFinder
on 57 broadly sampled Eukaryotes and distributed in 8 taxonomic groups (14
Heterokonts, 7 Alveolates, 8 Viridiplantae, 4 Rhodophytes, 3 Kinetoplastids, 3
Amoebozoans, 12 Fungi + Fonticulate and 6 Metazoans + Ichthyosporea). We then
filtered the XXXXXXX groups according to their complexity:
    - Taxonomic richness
        - Universal (every species are present)
        - Widespread (approximately 66% of species in each group)
        - Spread (approximately 33% of species in each group)
    - Average copy number
        - 1.1
        - 1.2
        - 2.0
We thus obtained 480 OGs distributed as follow:

            Max copy mean
           | 1.1      | 1.2      | 2.0
Universal  | 17       | 33       | 36
Widespread | 76       | 130      | 59
Sparse     | 32       | 60       | 37

### OGs reconstruction test

The main idea was to prune targeted sequences in all OGs and try to recover
the their original composition by running *Forty-Two* with a mining set
constituted of the proteomes of the species that were pruned. As for the
reference proteome set, we used all species but the ones concerned by the
pruning to avoid that the same sequences are used for the orhtology assessment
step.

We executed this protocol 8 times in parallel by pruning a different taxonomic
group each time and tried to recover the original group by running *Forty-Two*
with a mining set constituted of the (i) proteomes (ideal case scenario) and 
(ii) genomes (real-world case scenario) of the species that were pruned.

It was also the oportunity to test the influence of a parameter called
'ref_org_mul' which is a multiplier that allows to select best suited reference 
organisms for each OGs by ranking them by the average bit score of the obtained
best-hits when performing the BLAST from the OG's queries. We actually tested 10 
values of the 'ref_org_mul' ranging from 0.1 to 1. In the end, we performed 
38,400 pruning/enrichment experiment (480 OGs x 8 taxonomic groups x 10 
'ref_org_mul' values).

#### Proteomic data
When considering proteome as the mining set, a list comparison of protein 
identifiers was enough to assess if an OG had been restored as the original.


#### Genomic data

OGs reconstitution assessment was more challenging when considering genomes as the 
mining set for identifiers comparison was not possible. Instead, we had to 
build phylogenetic trees and develop a script for parsing bipartitions.

We computed precision and recall curves for 

## Assessing performance - F1-score

