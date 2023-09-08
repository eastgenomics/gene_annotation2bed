# gene_annotation2bed

Custom script for processing a list of ids (HGNC, transcript) or coordinates with associated annotation, into a comprehensive bed file for the corresponding refseq transcripts for each ID entry.

## What are typical use cases for this script?
- Converting a list of HGNC ids + associated gene level annotation information
  into a comprehensive bed file for annotation with Ensemble's VEP.
- Other use cases include providing different inputs such a list of transcripts.
  Or using exact coordinates to flag a regions such as TERT promoter.

## What data are required for this script to run?

- List of ids and annotation information in TSV format.
- Human Genome Reference (i.e. hs37d5)
- RefSeq Transcripts file (gff3)
## IGV reports output
Example Output:
[!image] (https://raw.githubusercontent.com/eastgenomics/gene_annotation2bed/dev/data/demo/demo_igv_reports.png)

The script produces a HTML report of all the bed file entries. Displayed in IGV with the refseq track
and bed file aligned with the respecive annotation.

## Script Inputs - Defaults & Behaviour

- `Genome` (required): The genome build for the resource
- `Refseq gff` (`--gff_file`): The corresponding gff file for refseq transcripts for the genome build.
- OR the processed dataframe for the refseq gff in pickle format (--pickle).
- annotation or transcript file with the annotation information in TSV format.
- The reference fasta for using for igv reports (i.e. `-fasta hs37d5.fa.gz`), the corresponding
  index should be present in the same folder.
- The output file suffix for the outputed .bed file.
- Flanking (int): The required flanking either side of the transcripts selected.
- Assembly summary - corresponding assembly report file for the refseq.gff, this is used
  to determine the corresponding chromosome for each transcript.

## Requirements

- pysam
- pandas
- igv-reports (v)
- numpy
- re

install using `requirements.txt`. `pip install requirements.txt`

## Running Script
