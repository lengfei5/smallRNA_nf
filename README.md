# smallRNA-meth pipeline

This pipeline analyzes smallRNA methylation data. It results in a count table.

## Get pipeline

Project: https://gitlab.com/tburk/smallRNA-meth  
Clone: `git clone https://gitlab.com/tburk/smallRNA-meth.git`

## Dependencies

The pipeline was run with tools of following versions:

* Nexflow (0.26.0): https://www.nextflow.iomo/
* Singularity (2.2): http://singularity.lbl.gov/

### Singularity container

The container can be created with following commands:

```
cd singularity
sudo singularity create -s 3000 smallRNA
sudo singularity bootstrap smallRNA smallRNA.singularity.def

#shrink it (https://github.com/singularityware/singularity/issues/623)
image=smallRNA
stripped_img=`tempfile --directory=.`
tail -n +2 $image > $stripped_img
sudo e2fsck -f $stripped_img
sudo resize2fs -M $stripped_img
shrunk_img=`tempfile --directory=.`
head -n 1 $image > $shrunk_img
cat $stripped_img >> $shrunk_img
rm $stripped_img
sudo mv $shrunk_img $image
chmod a+rx $image
```

## Parameter

Parameters are defined in `nextflow.config` and can be set as documented in Nextflow.

* outdir: Output base directory (default: ".")
* reads: Unaligned BAMs of miRNA experiment (default: "ngs_raw/*.bam")
* adapter: Adapter to trim at 3' (default: "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG")
* adapterER: Adapter error rate (default: 0.1)
* minLength: minimal length of trimmed/cut sequence (default: 18)
* maxLength: maximal length of trimmed/cut sequenc; to turn of set to -1 (default: 30)
* minAlign: minimal length of genome matchng sequence (default: 18)
* genome: FASTA of genome (default: "$baseDir/genome/WBcel235.fa")
* contamination: FASTA of contaminations or none (default: "none")
* gtf: GTF of miRNA (default: "$baseDir/annotation/WBcel235.mirBase.gtf")
* tailFraction: longest allowed tail fraction (defailt: 0.12)
* cpus: maximal number of CPUs to user (default: 8)

## Run pipeline

### C.elegans

* Genome: ftp://ftp.wormbase.org/pub/wormbase/releases/WS235/species/c_elegans/
* GTF: http://www.mirbase.org/ftp.shtml

Analysis as in publication:

```
nextflow run -with-singularity $PATH_SINGULARITY_CONTAINER/smallRNA \
   https://gitlab.com/tburk/smallRNA-meth.git \
   -r 6d254abb \
   --genome ../genome/WBcel235.fa \
   --gtf ../annotation/WBcel235.mirBase.gtf
```

### D.melanogaster

* Genome: ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r5.57_FB2014_03/fasta/
* GTF: ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r5.57_FB2014_03/gff/
  * converted to GTF and filtered for pre_miRNA features; alternative mirBase)

Analysis as in publication:

```
nextflow run -with-singularity $PATH_SINGULARITY_CONTAINER/smallRNA \
   https://gitlab.com/tburk/smallRNA-meth.git \
   -r 6d254abb \
   --genome ../genome/dmel-all-chromosome-r5.57.fasta \
   --contamination ../genome/dmS2_viruses.fa \
   --gtf ../annotation/flybase_r5.57.pre_miRNA.gtf
```

dmS2_viruses.fa contain sequences of following NCBI entries:

```
>gi|264980981|gb|GU144510.1| Mosquito nodavirus MNV-1 coat protein gene, partial cds
>gi|9629650|ref|NC_001834.1| Drosophila C virus, complete genome
>gi|22855185|ref|NC_004169.1| Drosophila x virus segment B, complete sequence
>gi|22855187|ref|NC_004177.1| Drosophila x virus segment A, complete sequence
>gi|346421290|ref|NC_007919.3| Nora virus, complete genome
>gi|253761971|ref|NC_012958.1| Drosophila A virus, complete genome
>gi|256535775|ref|NC_013135.1| Drosophila melanogaster sigma virus AP30 N, P, X, M, G and L genes, genomic RNA, isolate AP30
>gi|262225309|gb|GQ342965.1| Drosophila melanogaster American nodavirus (ANV) SW-2009a segment RNA1, complete sequence
>gi|262225312|gb|GQ342966.1| Drosophila melanogaster American nodavirus (ANV) SW-2009a segment RNA2, complete sequence
>gi|262225302|gb|GQ342962.1| Drosophila melanogaster birnavirus SW-2009a strain DBV segment A, complete sequence
>gi|262225305|gb|GQ342963.1| Drosophila melanogaster birnavirus SW-2009a strain DBV segment B, complete sequence
>gi|262225307|gb|GQ342964.1| Drosophila melanogaster tetravirus SW-2009a strain DTRV putative RNA-dependent RNA polymerase gene, complete cds
>gi|268053723|ref|NC_013499.1| Drosophila melanogaster totivirus SW-2009a, complete genome
>gi|22681055|ref|NC_004146.1| Flock house virus RNA 1, complete sequence
>gi|22711883|ref|NC_004144.1| Flock house virus, complete genome
```


## Results

* count: contains count tables for each input BAM with following columns
  * Name: miRNA name
  * GM: Genome matching
  * PM: Prefix matching
  * Total: GM + PM
* cutadapt: contains cutting log
* results:
  * countStatTable.txt: count statistics for each sample
  * countTable.txt: merged count tables (see above). Contains in addition the "miR" prefixed columns which indicating which miRNA is more abundant (TRUE = more abundant).


## Post processing

Differential miRNA analysis can be performed as shown in the "example" subfolder.


## Citation

Paper under revision.
