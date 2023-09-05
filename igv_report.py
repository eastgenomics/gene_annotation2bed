import json
import subprocess


def create_igv_report(bed_file, maf_file, genome, reference_file, info_columns, title, output_file):
    """
    Create an IGV report from a bed file.
    Parameters
    ----------
    bed_file : str
        file path to bed file.
    genome : str
        genome name. i.e hg19.
    info_columns : list
        list of column names to be displayed in the report.
    title : str
        title of the report.
    output_file : str
        file path to save the report.

    Returns
    -------
    Prints the standard output and error of the subprocess.
    Creates an IGV report.
    """
    tracks_config = [
        {
            "name": 'Genes',
            "type": '',
            "url": f'data/{genome}/refGene.txt.gz',
            "indexURL": f'data/{genome}/refGene.txt.gz.tbi'
        },
        {
            "name": 'BED',
            "type": '',
            "url": f'{bed_file}.gz', # /home/rswilson1/Documents/Programming_project/gene_annotation2bed/output_hg19_overlap_test7.bed.gz
            "indexURL": f'{bed_file}.gz.tbi'
        },
    ]
    # Serializing json
    tracks_json = json.dumps(tracks_config, indent=4)

    # Writing to sample.json
    with open("tracks_config.json", "w") as outfile:
        outfile.write(tracks_json)

    bgzip_result = subprocess.run(["bgzip", bed_file])
    print(bgzip_result.returncode)
    print("Standard Output:", bgzip_result.stdout)
    print("Standard Error:", bgzip_result.stderr)
    index_result = subprocess.run(["tabix", f"{bed_file}.gz"])
    print(index_result.returncode)
    print("Standard Output:", index_result.stdout)
    print("Standard Error:", index_result.stderr)
    if reference_file:
        maf_based_cmd = [
            "create_report",
            maf_file,
            "--fasta", reference_file,  # /home/rswilson1/Documents/MSC_dissertation/MELT_example_data/hs37d5.fa
            "--sequence", "1",
            "--begin", "2",
            "--end", "3",
            "--tracks", maf_file,
            "--track-config", "tracks_config.json",
            "--title", title,
            "--output", output_file
        ]

    elif genome:
        maf_based_cmd = [
            "create_report",
            maf_file,
            "--genome", genome,
            "--tracks", maf_file,
            "--track-config", "tracks_config.json",
            "--title", title,
            "--output", output_file
        ]

    # tsv_based_cmd = [
    #     "create_report",
    #     "/home/rswilson1/Documents/Programming_project/gene_annotation2bed/output_hg19_overlap_test5_edit.maf",
    #     "--sequence", "1",
    #     "--begin", "2",
    #     "--end", "3",
    #     "--fasta", "/home/rswilson1/Documents/MSC_dissertation/MELT_example_data/hs37d5.fa",
    #     "--zero_based", "true",
    #     #"--info-columns", "chr", "start", "end", "annotation", "gene", # info_columns
    #     "--tracks", bed_file,
    #     "--track-config", "TracksConfig.json",
    #     #"--genome", genome,
    #     "--title", title,
    #     "--output", output_file
    # ]

    result = subprocess.run(maf_based_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print(result.returncode)
    print("Standard Output:", result.stdout)
    print("Standard Error:", result.stderr)
    print(f"IGV Reports visualization saved as '{output_file}'")


if __name__ == "__main__":
    bed_file = "data/test.maflite.maf"
    genome = "hg19"
    title = f"TEST"
    output = "test.html"
    create_igv_report(bed_file, genome, title, output)
