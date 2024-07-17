[![Docker CI](https://github.com/tobiasrausch/vcLR/workflows/Docker%20CI/badge.svg)](https://hub.docker.com/r/trausch/vcLR/)
[![GitHub license](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://github.com/tobiasrausch/vcLR/blob/master/LICENSE)
[![GitHub Releases](https://img.shields.io/github/release/tobiasrausch/vcLR.svg)](https://github.com/tobiasrausch/vcLR/releases)

# Structural variant calling tutorial using long-reads.

## Installation - Linux

Installation of required libraries depends on your linux distribution, e.g., for Ubuntu:

`apt-get update`

`apt-get install -y autoconf build-essential cmake g++ git libcurl4-gnutls-dev libbz2-dev libdeflate-dev libgl1-mesa-dev libncurses-dev liblzma-dev pkg-config zlib1g-dev`

Then clone the repository and install all dependencies using mamba:

`git clone --recursive https://github.com/tobiasrausch/vcLR`

`cd vcLR`

`make all`

Download the course data:

ToDo

## Load the software stack

Load the mamba environment with all required tools

```bash
if [ ! -z ${CONDA_PREFIX+x} ]; then conda deactivate; fi
export PATH=`pwd`/mamba/bin:${PATH}
```

## SV Calling

### Reconstructing a derivative chromosome in cancer

In this practical we will analyze germline and somatic structural variants (SVs) of a [tumor cell line](https://www.nist.gov/programs-projects/cancer-genome-bottle) from a pancreatic ductal adenocarcinoma. The data was subsampled and subset to chr1 and chr5 to speed up all subsequent analysis. The tumor genome alignment file is named `tumor.bam` and the control genome alignment file is named `control.bam`.

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
You may also want to visualize clipped reads using the 'View' menu, then click 'Preferences...' and then click the 'Alignments' tab to switch on 'Show soft-clipped reads'. For the last complex SV, you may also want to highlight one specific read named `4190c914-8c3e-413b-88e3-caa07fca2861` which is over 100Kbp long: Use a right click on `tumor.bam` and then 'Select by name...'. We will later use that read to visualize this complex region.

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

In IGV, we already saw that a very long read named '4190c914-8c3e-413b-88e3-caa07fca2861' maps to that region. Let's do a dotplot of that read againt the reference region using [wally](https://github.com/tobiasrausch/wally).

```bash
wally dotplot -g genome.fa -r 4190c914-8c3e-413b-88e3-caa07fca2861 -e chr1:16500000-16650000 tumor.bam
```

### Delly structural variant calling

[Delly](https://github.com/dellytools/delly) is a method for detecting structural variants using short- or long-read sequencing data.
Using the tumor and normal genome alignment, delly calculates structural variants and outputs them as a BCF file, the binary encoding of [VCF](https://samtools.github.io/hts-specs). `lr` is delly's long-read SV discovery mode.

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
bcftools query -i 'SVTYPE=="INS" && SVLEN>250 && SVLEN<350' -f "%CHROM\t%POS\t%ALT[\t%GT]\n" sv.bcf  | grep "AAAAAAAAAAA"
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

