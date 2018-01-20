## Introduction
This is a basic vcf annotation tool. The future version will be updated if possible. 
Currently the tool only support the annotation based on http://exac.hms.harvard.edu RESTful API. The deleterious level of each variant is based on https://www.targetvalidation.org/variants.

## Usage:
### Annotate File
```
python annotate.py  argument1   [argument2]
```
1. argument1: input file path
2. argument2: output file path
3. input: vcf file
4. output: annotated vcf file saved in current working directory with name 'annotated.vcf' in default.

### Check output file
```
./checkfile.py  argument
```
argument: annotated file path

## Annotated Items
1. DP: Total read depth at the locus
2. AO: Alternate allele observations, with partial observations recorded fractionally
3. AO_RO: Percentage of reads supporting the variant versus those supporting reference reads
4. Consequence: Type of variation, such as transcript_ablation, stop_lost, et,al. 
If there are multiple possibilities, most deleterious is recorded.
5. ExACFrequency: Allele frequency of variant from Broad Institute ExAC Project API

