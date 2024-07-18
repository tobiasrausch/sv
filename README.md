[![Docker CI](https://github.com/tobiasrausch/sv/workflows/Docker%20CI/badge.svg)](https://hub.docker.com/r/trausch/sv/)
[![GitHub license](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://github.com/tobiasrausch/sv/blob/master/LICENSE)
[![GitHub Releases](https://img.shields.io/github/release/tobiasrausch/sv.svg)](https://github.com/tobiasrausch/sv/releases)

# Structural variant calling tutorial using long-reads.

## Installation - Linux

Installation of required libraries depends on your linux distribution, e.g., for Ubuntu:

`apt-get update`

`apt-get install -y autoconf build-essential cmake g++ git libcurl4-gnutls-dev libbz2-dev libdeflate-dev libgl1-mesa-dev libncurses-dev liblzma-dev pkg-config zlib1g-dev`

Then clone the repository and install all dependencies using mamba:

`git clone --recursive https://github.com/tobiasrausch/sv`

`cd sv`

`make all`

Download the course data:

`make FILE=1PfCy8yESCxvI8RJsfxTbF-QsygfnKNA2 download`

## Load the software stack

Load the mamba environment with all required tools.

```bash
if [ ! -z ${CONDA_PREFIX+x} ]; then conda deactivate; fi
export PATH=`pwd`/mamba/bin:${PATH}
```

## SV Calling

### Reconstructing a derivative chromosome in cancer

In this practical we will analyze structural variants (SVs) of a [pancreatic cancer cell line](https://www.nist.gov/programs-projects/cancer-genome-bottle) from the [Genome in a Bottle project](https://www.nature.com/articles/sdata201625). The data was subsampled and subset to chr1 and chr5 to speed up all subsequent analysis. The tumor genome alignment file is named `tumor.bam` and the control genome alignment file is named `control.bam`.

### Structural variant alignment quality control

For the discovery of structural variants, you should first assess the quality of the sequencing data as, for example, local assembly methods suffer from short reads and a high sequencing error rate and read-depth methods from uneven coverage. Common quality criteria are the percentage of reads mapped, the duplicate rate, the read-length distribution and the error rate. Popular tools to compute long-read quality control metrics are [NanoPlot](https://github.com/wdecoster/NanoPlot) and [Alfred](https://github.com/tobiasrausch/alfred).

```bash
cd data/lr/
alfred qc -r genome.fa -o qc.tsv.gz -j qc.json.gz tumor.bam
zcat qc.tsv.gz | grep ^ME | datamash transpose
NanoPlot --bam tumor.bam -o qc_tumor
cat qc_tumor/NanoStats.txt
```

As you can see from the QC results, the data has been downsampled to fairly low coverage to speed up all analyses in this tutorial.
This implies that some structural variants will have only weak support. In terms of QC interpretation, there are some general things to watch out for, such as unexpected high error rates (>4%), less than 80% of reads above Q10, an N50 read length below 10Kbp or unexpected patterns in the read length histogram.
 
#### Exercises

* What is the difference between the weighted read length histogram and the non-weighted histogram (WeightedHistogramReadlength.png, Non_weightedHistogramReadlength.png)?
* How is the N50 read length calculated?
* What is the median coverage of the tumor genome?

### Germline Structural Variants

Before we dive into structural variant (SV) calling, let's get an idea of how SVs look like in long-read sequencing data. I have prepared a [BED](https://bedtools.readthedocs.io/) file with some "simple" germline structural variants like deletions and insertions and one complex example.

```bash
cat svs.bed
```

Using [IGV](http://software.broadinstitute.org/software/igv/) we can browse these SVs interactively.

```bash
igv -g genome.fa
```

Once IGV has started use 'File' and 'Load from File' to load the `tumor.bam` and `control.bam` alignment file. Then import the `svs.bed` file from your working directory using 'Regions' and 'Import Regions'.
You can then easily navigate to the structural variants with 'Regions' and 'Region Navigator'.
Select a structural variant in the region navigator and click 'View', which will center the IGV alignment view on the selected structural variant.
You can zoom in and out using the '+' and '-' signs in the toolbar at the top.
You may also want to visualize clipped reads using the 'View' menu, then click 'Preferences...' and then click the 'Alignments' tab to switch on 'Show soft-clipped reads'. For the last complex SV, you may also want to highlight one specific read named `8c0eba99-ad7d-4246-87f0-46177b9946d4` which is over 100Kbp long: Use a right click on `tumor.bam` and then 'Select by name...'. We will later use that read to visualize this complex region.

#### Exercises

* Which of the two insertions could be a mobile element insertion? What typical features of a mobile element can you observe for that insertion?
* Do you see for heterozygous SVs, potential heterozygous SNPs in the vicinity that are on the same haplotype?

### Plotting structural variants

IGV is excellent for interactive browsing but for large numbers of SVs you can use command-line tools such as [wally](https://github.com/tobiasrausch/wally) to plot multiple SVs in batch.

```bash
wally region -R svs.bed -cp -g genome.fa tumor.bam control.bam
```

### Complex structural variants

Even in germline genomes we can observe complex structural variants or potential reference errors. One example was in the `svs.bed` file.

```bash
cat svs.bed | grep "complex"
```

In IGV, we already saw that a very long read named `8c0eba99-ad7d-4246-87f0-46177b9946d4` maps to that region.
Let's do a dotplot of that read against the reference region using [wally](https://github.com/tobiasrausch/wally).

```bash
wally dotplot -g genome.fa -r 8c0eba99-ad7d-4246-87f0-46177b9946d4 -e chr1:16500000-16600000 tumor.bam
```

This dotplot suggests a repeat expansion relative to GRCh38.

#### Exercises

* The [UCSC genome browser](https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr1%3A16500000%2D16650000) suggests that there is already a patch of the reference for this area of the genome, called chr1_MU273333v1_fix.fa. You may want to try and realign the read to the patched genome using minimap2.

```bash
samtools view -F 3844 tumor.bam chr1:16500000-16650000 | grep "^8c0eba99-ad7d-4246-87f0-46177b9946d4" | cut -f 1,10 | awk '{print ">"$1"\n"$2;}' > read.fa
minimap2 -ax map-ont -L chr1_MU273333v1_fix.fa read.fa | samtools sort -o read.bam -
samtools index read.bam
wally dotplot -g chr1_MU273333v1_fix.fa -r 8c0eba99-ad7d-4246-87f0-46177b9946d4 -e chr1_MU273333v1_fix:300000-400000 read.bam
```

### Delly structural variant calling

[Delly](https://github.com/dellytools/delly) is a method for detecting structural variants using short- or long-read sequencing data.
Using the tumor and normal genome alignment, delly calculates structural variants and outputs them as a BCF file, the binary encoding of [VCF](https://samtools.github.io/hts-specs). Delly's long-read SV discovery mode uses the subcommand `lr`.

```bash
delly lr -g genome.fa -o sv.bcf tumor.bam control.bam
```

#### VCF encoding of structural variants

VCF was originally designed for small variants such as single-nucleotide variants (SNVs) and short insertions and deletions (InDels). That's why all SV callers heavily use the VCF INFO fields to encode additional information about the SV such as the structural variant end position (INFO:END) and the SV type (INFO:SVTYPE). You can look at the header of the BCF file using grep where '-A 2' includes the first two structural variant records after the header in the file:

```bash
bcftools view sv.bcf | grep "^#" -A 2
```

[Delly](https://github.com/dellytools/delly) uses the VCF:INFO fields for structural variant site information, such as how confident the structural variant prediction is and how accurate the breakpoints are. The genotyping fields contain the actual sample genotype, its genotype quality and genotype likelihoods and various count fields for the variant and reference supporting reads.
Please note that at this stage the BCF file contains germline and somatic structural variants but also false positives caused by mis-mappings or incomplete reference sequences.

#### Querying VCF files

[Bcftools](https://github.com/samtools/bcftools) offers many possibilities to query and reformat SV calls. For instance, to output a table with the chromosome, start, end, identifier and genotype of each SV we can use:

```bash
bcftools query -e 'SVTYPE=="BND"' -f "%CHROM\t%POS\t%INFO/END\t%ID[\t%GT]\n" sv.bcf | head
```

Inter-chromosomal translocations with SV type `BND` are a special case because they involve two different chromosomes.

```bash
bcftools query -i 'SVTYPE=="BND"' -f "%CHROM\t%POS\t%INFO/CHR2\t%INFO/POS2\t%ID[\t%GT]\n" sv.bcf
```

Bcftools can also filter by genomic regions. For instance, to verify delly found the Alu (mobile element) insertion we previously looked at in IGV, you can use:

```bash
grep "INS02" svs.bed
bcftools view sv.bcf chr5:56632119-56632158
```

Delly's consensus sequence (INFO:CONSENSUS) is a local assembly of all SV-supporting reads. So we can again create a dotplot using [wally](https://github.com/tobiasrausch/wally) to highlight the insertion relative to GRCh38.

```bash
bcftools query -f "%ID\t%INFO/CONSENSUS\n" sv.bcf | grep "INS00001302" | awk '{print ">"$1"\n"$2;}' > ins.fa
samtools faidx genome.fa chr5:56631000-56633000 | sed 's/^>.*$/>hg38/' >> ins.fa
wally dotplot ins.fa
```

To identify additional Alu candidate insertions, we can for instance filter by size and select insertions with a poly-A tail.

```bash
bcftools query -i 'SVTYPE=="INS" && SVLEN>250 && SVLEN<350' -f "%CHROM\t%POS\t%ALT\n" sv.bcf  | grep "AAAAAAAAAAA"
```

#### Exercises

* How many inter-chromosomal translocations were identified by delly?
* How can bcftools be used to count the number of structural variants for the different SV types (DEL, INS, DUP, INV, BND)?
* Using the bcftools query command above: How many potential Alu insertions are in forward and reverse orientation?

### Somatic structural variant filtering

Delly's initial SV calling cannot differentiate somatic and germline structural variants. We therefore now use delly's somatic filtering, which requires a sample file listing tumor and control sample names from the VCF file.

```bash
cat spl.tsv
```

There are many parameters available to tune the somatic structural variant filtering like the minimum variant allele frequency to filter out subclonal variants, for instance. Let's leave the defaults for the time being and only filter for SVs with the VCF FILTER field set to PASS.

```bash
delly filter -p -f somatic -o somatic.bcf -s spl.tsv sv.bcf
```
As expected, the somatic SVs have a homozygous reference genotype in the control sample.

```bash
bcftools query -f "%CHROM\t%POS\t%INFO/END\t%ID[\t%GT]\n" somatic.bcf
```

Using [Bcftools](https://github.com/samtools/bcftools) and [wally](https://github.com/tobiasrausch/wally) we can also easily plot the intra-chromosomal somatic SVs.

```bash
bcftools query -e 'SVTYPE=="BND"' -f "%CHROM\t%POS\t%INFO/END\t%ID\n" somatic.bcf | awk '{print $1"\t"($2-50)"\t"($3+50)"\t"$4;}' > somatic.bed
wally region -R somatic.bed -cp -g genome.fa tumor.bam control.bam
```

#### Exercises

* Do you think all somatic variants are truly somatic? Which ones are likely false positive somatic SVs?

### Reconstructing a derivative chromosome in cancer

[IGV](http://software.broadinstitute.org/software/igv/) and [wally](https://github.com/tobiasrausch/wally) are good for relatively small SVs but for large SVs like the duplication-type SV or inter-chromosomal translocations we need to integrate read-depth with structural variant predictions to get a better overview of complex somatic rearrangements. Let's first create a simple read-depth plot.

```bash
delly cnv -i 50000 -j 50000 -w 50000 -uo cnv.bcf -c cnv.cov.gz -g genome.fa -m map.fa.gz tumor.bam
Rscript cnBafSV.R cnv.cov.gz
```
Now we can overlay the somatic structural variants on top of the read-depth information.

```bash
bcftools query -f "%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%ID\t%INFO/CHR2\t%INFO/POS2\n" somatic.bcf > svs.tsv
Rscript cnBafSV.R cnv.cov.gz svs.tsv
```

Apparently, chr1 and chr5 are connected but the question is whether the left or the right end is joined from each inter-chromosomal translocation breakpoint. Delly outputs the orientation of each segment in the INFO:CT field.

```bash
bcftools query -i 'SVTYPE=="BND"' -f "%INFO/CT\n" somatic.bcf
```

In this case, `3to3` indicates that chr1p is joined with chr5p in inverted orientation. Here is a brief summary of the different [connection types](https://github.com/tobiasrausch/wally#paired-end-view).


#### Exercises

* Does the derivative chromosome containing chr1p and chr5p contain a centromere?

### Calling and phasing small variants (SNVs and InDels)

For phasing and haplo-tagging the tumor reads we first call small variants using [longshot](https://github.com/pjedge/longshot) in the matched control.

```bash
longshot --no_haps -s HG008-N-P --bam control.bam --ref genome.fa --out control.longshot.vcf
bgzip control.longshot.vcf
tabix control.longshot.vcf.gz
```

We can then perform long-read based haplotyping using [WhatsHap](https://github.com/whatshap/whatshap) and phase heterozygous variants.

```bash
whatshap phase --ignore-read-groups --reference genome.fa control.longshot.vcf.gz control.bam -o control.whatshap.vcf
bgzip control.whatshap.vcf
tabix control.whatshap.vcf.gz
```

Let's briefly inspect some phasing statistics using the [WhatsHap](https://github.com/whatshap/whatshap) `stats` subcommand.

```bash
whatshap stats --chr-lengths <(cut -f 1,2 genome.fa.fai) control.whatshap.vcf.gz
```

#### Exercises

* How many mega-basepairs is the longest phased block?
* What is the N50 phased block length?

### Haplotagging of tumor reads

We can now use the haplotypes derived from the matched control to phase the tumor reads.

```bash
whatshap haplotag -o tumor.hp.bam --reference genome.fa --ignore-read-groups --tag-supplementary control.whatshap.vcf.gz tumor.bam
samtools index tumor.hp.bam
```

The `tumor.hp.bam` file now has the `HP` tag included for all the reads that could be phased to haplotype 1 or 2. [IGV](http://software.broadinstitute.org/software/igv/) can group and color reads by `HP` tag so let's revisit the SVs from the `svs.bed` file that we used in the beginning of the tutorial.

```bash
igv -g genome.fa
```

As previously, once IGV has started use 'File' and 'Load from File' to load the `tumor.hp.bam` alignment file with the `HP` haplotype tags. Then import the `svs.bed` file from your working directory using 'Regions' and 'Import Regions'. You can then easily navigate to the structural variants with 'Regions' and 'Region Navigator'. Once you have selected an SV, right click on `tumor.hp.bam`, use 'Color alignments by' and select 'tag'. Then enter 'HP' to color by haplotype. Alternatively, you can also use 'Group alignments by' tag to sort reads according to assigned haplotype.

#### Exercises

* Did all heterozygous SVs separate correctly by haplotype?

### Variant allele frequency of small variants across somatic copy-number alterations

Let's annotate the small variants of the control genome in the matched tumor.

```bash
longshot --no_haps --potential_variants control.longshot.vcf.gz --min_cov 3 --min_alt_count 1 --min_alt_frac 0.01 -s HG008-T --bam tumor.bam --ref genome.fa --out tumor.longshot.vcf
bgzip tumor.longshot.vcf
tabix tumor.longshot.vcf.gz
```

