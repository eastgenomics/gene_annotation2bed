import subprocess

def create_igv_report(bed_file, genome, info_columns, title, output_file):
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
    bed_based_cmd = [
        "create_report",
        bed_file,
        "--fasta", "/home/rswilson1/Documents/MSC_dissertation/MELT_example_data/hs37d5.fa",
        "--tracks", bed_file,
        "--track-config", "TracksConfig.json",
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

    result = subprocess.run(bed_based_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
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
