# MxFinder
Testing of mutual exclussion between pairs of exons of nonsense mediated decay (NMD)

MxFinder is a pipeline designed for checking of mutual exclussion between pairs of exons at the time of triggering NMD, just relying in an annotated transcriptome.

These are the main scripts of the pipeline:

- format_gtf_refseq.py: given a gtf and a gene, extract the info assocated to this gene
and transform it to bed file. There is another version available for Ensembl annotation (format_gtf.py)
- get_fasta.py: from [MoSEA](https://github.com/comprna/MoSEA) software, extract the fasta sequence associated to each exon
- create_paths.py: get all the possible paths with the info given. It returns the fasta sequence associated to each combination.
- extract_orfs.py: Get all possible ORFs associated to each sequence
- get_distance_to_ss.py: Given the ORFs per transcript (take the longest per transcript), the sequences of each transcript and the position of the ss, get the relative distance to this ss
- test_NMDs.py: apply a Fisher test for each pair of exons, in order to test independence between all the exons respect to the NMD condition

The user can run this scripts separately, or run it all in one line as follows:

```
python3 MxFinder.py <gtf-file> <gene> <python2-executable> <mosea-path> <genome-fasta> <bedtools-path> <output-path>
```

Here we show an example of execution, with mm10 annotation and a gene (Mbnl1):

```
python3 MxFinder.py ~/refseq_mm10_full.formatted.gtf Mbnl1 ~/python2.7/bin/python2.7 ~/mosea.py ~/mm10.fa ~/bedtools-2.26/bin/bedtools ~/MxFinder_mm10_output
```
