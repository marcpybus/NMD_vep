# NMD_rules

This is a plugin for Ensembl Variant Effect Predictor (VEP) software that predicts whether a stop gained variant triggers nonsense-mediated 
mRNA decay (putative_NMD_triggering) or not (canonical_NMD_escaping, noncanonical_NMD_escaping) based on the following rules (1):

1. If the variant is located in an intronless transcript, meaning only one exon exist in the transcript (canonical_NMD_escaping;intronless). 
2. The variant location  falls in the last exon of the transcript (canonical_NMD_escaping;last_exon).
3. The variant location falls  50 bases upstream of the penultimate exon (canonical_NMD_escaping;50bp_penult_exon).
4. The variant falls in the first 150 coding bases in the transcript (noncanonical_NMD_escaping;first_150bp). 
5. If the variant is in an exon bigger than 407 bp (noncanonical_NMD_escaping;lt_407bp_exon). 

The plugin also shows the 2 codons/aminoacids before the novel stop codons, as well as the next nucleotide for futher analysis, as they may influence the NMD efficiency (2):
    i.e. GCC(Ala)CTG(Leu)TGA(Stop)C

Finally, it computes the distance in aminoacids to the next Metionine (start codon) after the stop codon to leverage an hypothetical translation reinitiation (3). 

The plugin analyzes variants that contain a defined Ter in their HGVSp notation (stop_gained, stop_loss, frameshift). Exceptions are synonymous variants at the stop codon (i.e. p.Ter811=) or when a stop codon is not inferrer downstream of a frameshift variant (i.e. p.Ter257GlufsTer?). Similarly, splicing variants are no considered, as VEP does not infer the consequence on the protein.

```
NMrules:                   Considering 5 rules for NMD prediction of stop gain variants (Ter in HGVSp notation) & their surronding genomic context
Annotation format:         NMD_prediction:rule_used:-2codon(-2aa)-1codon(-1aa)stop_codon(Stop)fourth_letter:distance_in_aa_to_next_Met
```

This plugin has been run on the gnomad.exomes.r2.1.1.sites.vcf file without throwing any errors (https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz). Beware that more complex variants may trigger novel bugs.

REFERENCES :
(1) The impact of nonsense-mediated mRNA decay on genetic disease, gene editing and cancer immunotherapy (Lindeboomet al. 2019, Nature Genetics)
(2) Systematic analysis of nonsense variants uncovers peptide release rate as a novel modifier of nonsense-mediated mRNA decay efficiency (Kolokada et al. 2024, bioRxiv)
(3) Advanced variant classification framework reduces the false positive rate of predicted loss-of-function variants in population sequencing data (Singer-Berk et al. 2023, Am J Hum Genet)


