[![Docker CI](https://github.com/tobiasrausch/sv/workflows/Docker%20CI/badge.svg)](https://hub.docker.com/r/trausch/sv/)
[![GitHub license](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://github.com/tobiasrausch/sv/blob/main/LICENSE)
[![GitHub Releases](https://img.shields.io/github/release/tobiasrausch/sv.svg)](https://github.com/tobiasrausch/sv/releases)

# Structural variant calling tutorial using long-reads.

In this practical we reconstruct a derivative chromosome in cancer using long reads. The data comes from the [HG008 cancer cell line](https://www.nist.gov/programs-projects/cancer-genome-bottle) of the [Cancer Genome in a Bottle](https://www.nist.gov/programs-projects/cancer-genome-bottle) project. The data was subsampled and subset to `chr1` and `chr5` to keep all analyses fast.

## Installation

[pixi](https://pixi.sh) can be used to install all required tools and download the course data.

```bash
curl -fsSL https://pixi.sh/install.sh | bash
git clone --recursive https://github.com/tobiasrausch/sv
cd sv
pixi install
pixi shell
FILE=1PfCy8yESCxvI8RJsfxTbF-QsygfnKNA2 pixi run download
```

Alternatively, you can also use the pre-built docker image with JupyterLab to run the practical in a web browser.

```bash
docker pull trausch/sv:latest
docker run -it -p 8888:8888 trausch/sv:latest
```

In JupyterLab, you then need to download the data using `FILE=1PfCy8yESCxvI8RJsfxTbF-QsygfnKNA2 pixi run download`.

In courses, I use pre-built AWS course images with the data so you can directly launch the container with the mounted data.

```bash
ssh -L 8888:localhost:8888 ubuntu@<host>
docker run -it -p 8888:8888 -v /data/lr:/opt/sv/data/lr trausch/sv:latest
```

## SV Calling

### Reconstructing a derivative chromosome in cancer

The tutorial data was subsampled and subset to chr1 and chr5.  The tumor genome alignment file is named `tumor.bam` and the control genome alignment file is named `control.bam`. The BAM files are so-called modBAM files with methylation information and reads have been tagged by parental haplotype.

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
This implies that some structural variants will have only weak support. In terms of QC interpretation, there are some general things to watch out for, such as unexpected high error rates (>4%), less than 80% of reads above Q10, an N50 read length below 10Kbp or unexpected patterns in the read length histogram. Alfred also reports phasing metrics such as `#HaploTagged`, `FractionHaploTagged`, `#PhasedBlocks` and `N50PhasedBlockLength`. You can explore the full interactive report by uploading `qc.json.gz` to the [Alfred web app](https://www.gear-genomics.com/alfred/).
 
#### Exercises

* What is the median coverage of the tumor genome?
* What fraction of tumor reads could be haplotagged and what is the N50 phased block length?
* How is the N50 read length calculated?

### Germline Structural Variants

Before we dive into structural variant (SV) calling, let's get an idea of how SVs look like in long-read sequencing data. I have prepared a [BED](https://bedtools.readthedocs.io/) file with some "simple" germline structural variants like deletions and insertions and one complex example.

```bash
cat svs.bed
```

From JupyterLab, you can open the `igv.ipynb` notebook or use the hosted [IGV web app](https://igv.org/app/). Then load `genome.fa`, `tumor.hp.bam`, `control.hp.bam` and `svs.bed` and then color and group alignments by the `HP` tag.

Alternatively, you can use [wally](https://github.com/tobiasrausch/wally) to generate plots on the command line.

```bash
wally region -R svs.bed -cp -g genome.fa tumor.hp.bam control.hp.bam
```

#### Exercises

* Which of the two insertions could be a mobile element insertion? What typical features of a mobile element can you observe for that insertion?
* For the heterozygous SVs, do nearby heterozygous SNPs "tag" the SV (same `HP`)?

### Methylation

The reads were basecalled with a 5mC model, so each CpG carries a methylation probability. [wally](https://github.com/tobiasrausch/wally) can plot this methylation information onto the alignments with its modified-base view. Each CpG is colored from blue (unmethylated) to red (methylated).

```bash
wally region -m 5mC -cp -g genome.fa -r chr1:789000-790200:methylation tumor.hp.bam
```

### Delly structural variant calling

[Delly](https://github.com/dellytools/delly) is a method for detecting structural variants using short- or long-read sequencing data.
Using the tumor and normal genome alignment, delly calculates structural variants and outputs them as a BCF file, the binary encoding of [VCF](https://samtools.github.io/hts-specs). Delly's long-read SV discovery mode uses the subcommand `lr`.

```bash
delly lr -y ont -g genome.fa -o sv.bcf tumor.hp.bam control.hp.bam
```

#### VCF encoding of structural variants

VCF was originally designed for small variants such as single-nucleotide variants (SNVs) and short insertions and deletions (InDels). That's why all SV callers heavily use the VCF INFO fields to encode additional information about the SV such as the structural variant end position (INFO:END) and the SV type (INFO:SVTYPE). You can look at the header of the BCF file using grep where '-A 2' includes the first two structural variant records after the header in the file:

```bash
bcftools view sv.bcf | grep "^#" -A 2
```

[Delly](https://github.com/dellytools/delly) uses the VCF:INFO fields for structural variant site information, such as how confident the structural variant prediction is and how accurate the breakpoints are. The genotyping fields contain the actual sample genotype, its genotype quality and genotype likelihoods and various count fields for the variant and reference supporting reads.
Please note that at this stage the BCF file contains germline and somatic structural variants but also false positives caused by mis-mappings or incomplete reference sequences.

#### Querying VCF files

[Bcftools](https://github.com/samtools/bcftools) offers many possibilities to query and reformat SV calls. For instance, to output a table with the chromosome, start, end, identifier, phaed genotype and haplotype count of each SV we can use:

```bash
bcftools query -e 'SVTYPE=="BND"' -f "%CHROM\t%POS\t%INFO/SVTYPE\t%ID[\t%GT\t%PS\t%HP]\n" sv.bcf | head
```

Inter-chromosomal translocations with SV type `BND` are a special case because they involve two different chromosomes.

```bash
bcftools query -i 'SVTYPE=="BND"' -f "%CHROM\t%POS\t%INFO/CHR2\t%INFO/POS2\t%ID[\t%GT]\n" sv.bcf
```

To confirm delly recovered the Alu (mobile-element) insertion we saw earlier:

```bash
grep "INS02" svs.bed
bcftools view sv.bcf chr5:56632119-56632158
```

Delly's consensus sequence (INFO:CONSENSUS) is a local assembly of all SV-supporting reads. So we can again create a dotplot using [wally](https://github.com/tobiasrausch/wally) to highlight the insertion relative to GRCh38.

```bash
bcftools query -f "%POS\t%ID\t%INFO/CONSENSUS\n" sv.bcf | grep "^566321" | awk '{print ">"$2"\n"$3;}' > ins.fa
samtools faidx genome.fa chr5:56631000-56633000 | sed 's/^>.*$/>hg38/' >> ins.fa
wally dotplot ins.fa
```

Delly directly annotates SV subtypes so you can for instance count the different types of mobile elements.

```bash
bcftools query -i 'SVTYPE=="INS"' -f "%SVTYPE\t%SUBTYPE\n" sv.bcf  | sort | uniq -c
```

#### Exercises

* How many inter-chromosomal translocations were identified by delly?
* How can bcftools be used to count the number of structural variants for the different SV types (DEL, INS, DUP, INV, BND)?
* How many potential Alu insertions are in forward and reverse orientation?

## Allele-specific methylation at structural variants

For each SV, delly reports methylation separately for the reference and alternative allele (`MR` vs `MA`, each with four windows around the SV start and end breakpoints). To find insertions where the alternative allele is methylated, we can use for instance:

```bash
bcftools query -f "%CHROM\t%POS\t%INFO/SVTYPE\t%ID[\t%GT\t%MR\t%MA]\n" sv.bcf \
  | awk -F'\t' '$3=="INS" && $7 ~ /,9[0-9],/'
```

A typical example is the below ~2.3 Kbp insertion where the reference is midly methylated but the inserted sequence is almost fully methylated.

```bash
wally region -m 5mC -cp -g genome.fa -r chr1:789000-790200:methylation tumor.hp.bam
```

#### Exercises

* Why might a newly inserted repeat element be methylated? What would you expect for an actively transcribed L1 element instead?
* Compare `MR` and `MA` for the same SV in tumor vs control. Is any methylation difference germline or tumor-specific?

### Somatic structural variant filtering

Delly's initial SV calling cannot differentiate somatic and germline structural variants. We therefore now use delly's somatic filtering, which requires a sample file listing tumor and control sample names from the VCF file.

```bash
cat spl.tsv
delly filter -p -f somatic -o somatic.bcf -s spl.tsv sv.bcf
```

There are many parameters available to tune the somatic structural variant filtering like the minimum variant allele frequency to filter out subclonal variants, for instance. As expected, the somatic SVs have a homozygous reference genotype in the control sample.

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
delly cnv -w 50000 -o cnv.bcf -c cnv.cov.gz -g genome.fa tumor.hp.bam
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

### B-allele frequency across somatic copy-number changes

As the control is already phased, we can annotate heterozygous SNPs with their allelic depth in the tumor to read out the so-called B-allele frequency.


```bash
bcftools view -v snps -i '(GT="0|1" || GT="1|0")' control.phased.vcf.gz -O b -o control.het.bcf
bcftools index control.het.bcf
bcftools query -f '%CHROM\t%POS\n' control.het.bcf | awk 'NR%7==1' | bgzip > sites.tsv.gz
tabix -s1 -b2 -e2 sites.tsv.gz
```

As we only need read counts and no variant calls, we can simply use `bcftools mpileup` on these variant sites.

```bash
bcftools mpileup -f genome.fa -a FORMAT/AD -T sites.tsv.gz tumor.hp.bam -Oz -o tumor.ad.vcf.gz
tabix tumor.ad.vcf.gz
bcftools merge -O b -o tumor.control.bcf tumor.ad.vcf.gz control.het.bcf
```

Now we split the tumor allelic depths by the control haplotype (`0|1` vs `1|0`) and plot he BAF signal.

```bash
bcftools query -f "%CHROM\t%POS[\t%GT\t%AD]\n" tumor.control.bcf | grep "0|1" | cut -f 1,2,4 | sed 's/,/\t/g' | awk '$3+$4>0 {print $1"\t"$2"\t"($3/($3+$4));}' > var.vaf
bcftools query -f "%CHROM\t%POS[\t%GT\t%AD]\n" tumor.control.bcf | grep "1|0" | cut -f 1,2,4 | sed 's/,/\t/g' | awk '$3+$4>0 {print $1"\t"$2"\t"($4/($3+$4));}' >> var.vaf
Rscript cnBafSV.R cnv.cov.gz svs.tsv var.vaf
```

#### Exercises

* For the somatic duplication, we have an estimated total copy-number of 3. What are the expected B-allele frequencies in that region?

