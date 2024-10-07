## NMDrules

### Description

This is a plugin for Ensembl Variant Effect Predictor (VEP) software that predicts whether a stop gained variant triggers nonsense-mediated 
mRNA decay (putative_NMD_triggering) or not (canonical_NMD_escaping, noncanonical_NMD_escaping) based on the following rules (1):

1. The variant is located in an intronless transcript, i.e. there is only one exon exist in the transcript (canonical_NMD_escaping:intronless). 
2. The variant is located in the last exon of the transcript (canonical_NMD_escaping:last_exon).
3. The variant is located within teh last 50 bases of the penultimate exon (canonical_NMD_escaping:50bp_penult_exon).
4. The variant is located within the first 150 coding bases in the transcript (noncanonical_NMD_escaping:first_150bp). 
5. The variant is located in an exon larger than 407 bp (noncanonical_NMD_escaping:lt_407bp_exon). 

The plugin also shows the 2 codons/amino acids before the novel stop codon, plus the next nucleotide, for a detailed analysis, as they may influence the NMD efficiency (2):
    i.e. GCC(Ala)CTG(Leu)TGA(Stop)C

Finally, it calculates the distance in amino acids to the next metionine (Met, start codon) after the stop codon to leverage an hypothetical translation reinitiation (3). 

The plugin analyses variants that contain a stop codon Ter in their HGVSp notation (stop_gained, stop_loss, frameshift). Exceptions are synonymous variants at the stop codon (e.g. p.Ter811=) or when a stop codon is not inferrer downstream of a frameshift variant (e.g. p.Ter257GlufsTer?). Similarly, splicing variants are not considered, as VEP does not infer the consequence on the protein.

Annotation format:
> NMD_prediction:rule_used:-2codon(-2aa)-1codon(-1aa)stop_codon(Stop)fourth_letter:distance_in_aa_to_next_Met

This plugin has been run on the gnomad.exomes.r2.1.1.sites.vcf file without throwing any errors (https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz). Be aware that more complex variants may cause new bugs.

REFERENCES :
(1) The impact of nonsense-mediated mRNA decay on genetic disease, gene editing and cancer immunotherapy (Lindeboomet al. 2019, Nature Genetics)
(2) Systematic analysis of nonsense variants uncovers peptide release rate as a novel modifier of nonsense-mediated mRNA decay efficiency (Kolokada et al. 2024, bioRxiv)
(3) Advanced variant classification framework reduces the false positive rate of predicted loss-of-function variants in population sequencing data (Singer-Berk et al. 2023, Am J Hum Genet)

### Usage

```
mv NMDrules.pm  ~/.vep/Plugins
./vep -i variations.vcf --plugin NMDrules
```

@marcpybus - https://github.com/marcpybus/

