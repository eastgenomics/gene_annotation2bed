# gene_annotation2bed

Purpose:
To provide bed files for custom gene-level annotation with VEP.
This custom script processes a list of ids (HGNC, transcript) or coordinates with associated annotation, into a comprehensive bed file for the corresponding refseq transcripts for each ID entry.

![Workflow diagram showing TSV containing IDs and annotation to bed file and how it is used in VEP and visualised in IGV using a VCF](https://raw.githubusercontent.com/eastgenomics/gene_annotation2bed/sprint_2/Workflow.png)

## What are typical use cases for this script?

- Converting a list of HGNC ids + associated gene level annotation information
  into a comprehensive bed file for annotation with Ensemble's VEP.
- Other use cases include providing different inputs such a list of transcripts.
  Or using exact coordinates to flag a regions such as TERT promoter.

## What data are required for this script to run?

- List of ids and annotation information in TSV format.
- Human Genome Reference (i.e. hs37d5)
- RefSeq Transcripts file (gff3) from 001_reference

## IGV reports output

IGV report:
![image](<https://raw.githubusercontent.com/eastgenomics/gene_annotation2bed/dev/data/demo/demo_igv_reports.png>)

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


## Strucute of code

├── data
│   ├── demo
│   │   ├── after_table.tsv.gz
│   │   ├── before_table.tsv.gz
│   │   ├── demo_igv_reports.png
│   │   └── initial_table.tsv.gz
│   ├── GCF_000001405.25_GRCh37.p13_assembly_report.txt
│   ├── GCF_000001405.25_GRCh37.p13_genomic.gff
│   ├── hg19
│   │   ├── ncbiRefSeq.txt.gz
│   │   ├── ncbiRefSeq.txt.gz.tbi
│   │   ├── refGene.txt.gz
│   │   └── refGene.txt.gz.tbi
│   └── hg38
│       ├── ncbiRefSeq.txt.gz
│       ├── ncbiRefSeq.txt.gz.tbi
│       ├── refGene.txt.gz
│       └── refGene.txt.gz.tbi
├── gene_annotation2bed.py (MAIN SCRIPT)
├── LICENSE
├── output_new_test.vcf
├── README.md
├── requirements.txt
├── scripts
│   ├── construct_vcf.py
│   └── igv_report.py
├── testing_090124
├── tests
│   ├── __init__.py
│   ├── test_construct_vcf.py
│   ├── test_data
│   │   ├── coordinates_anno_test.tsv
│   │   ├── example_bed_hg38.bed
│   │   ├── expected_output.vcf
│   │   ├── hgcn_ids_anno_test.tsv
│   │   ├── hs37d5.fa
│   │   ├── hs37d5.fa.fai
│   │   ├── refseq_gff_preprocessed.pkl
│   │   ├── test_empty_attributes.gff
│   │   ├── test_empty.gff
│   │   ├── test_missing_attributes.gff
│   │   └── transcripts_anno_test.tsv
│   ├── test_gene_annotation2bed.py
│   ├── test_gff_parsing.py
│   └── test_igv_report.py
├── tracks_config.json (not used)
├── utils
│   ├── configure_gff.py
│   ├── gff2pandas.py
│   ├── __init__.py
└── Workflow.png
