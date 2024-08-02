"""
Generate an IGV report from a bed file.
"""

import json
import subprocess


def create_igv_report(bed_file: str, maf_file: str,
                      genome: str, reference_file: str,
                      title: str,
                      output_file: str) -> None:
    """
    Create an IGV report from a bed file.
    Parameters
    ----------
    bed_file : str
        file path to bed file.
    maf_file : str
        file path to maf file.
    genome : str
        genome name. i.e hg19.
    title : str
        title of the report.
    output_file : str
        file path to save the report.

    Other required files:
    ---------------------
    ncibRefSeq.txt.gz and ncibRefSeq.txt.gz.tbi
    Download from:
    https://s3.amazonaws.com/igv.org.genomes/hg19/ncbiRefSeq.txt.gz.tbi
    https://s3.amazonaws.com/igv.org.genomes/hg19/ncbiRefSeq.txt.gz
    hg38 also available by changing hg19 to hg38 in the above links.

    Returns
    -------
    None.
    Prints the standard output and error of the subprocess.
    Creates an IGV report file.
    """
    tracks_config = [
        {
            "name": 'Genes',
            "type": '',
            "url": f'data/{genome}/ncbiRefSeq.txt.gz',
            "indexURL": f'data/{genome}/ncbiRefSeq.txt.gz.tbi'
        },
        {
            "name": 'BED',
            "type": '',
            "url": f'{bed_file}.sorted.gz',
            "indexURL": f'{bed_file}.gz.tbi'
        },
    ]
    # Serializing json
    tracks_json = json.dumps(tracks_config, indent=4)

    # Writing to sample.json
    with open("tracks_config.json", "w") as outfile:
        outfile.write(tracks_json)
    bed_file_sorted = f"{bed_file}.sorted"
    sort_result = subprocess.run(["sort", "-k1,1", "-k2,2n", "-k3,3n", bed_file, "-o", bed_file_sorted])
    print("Standard Output:", sort_result.stdout)
    print("Standard Error:", sort_result.stderr)
    bgzip_result = subprocess.run(["bgzip", bed_file_sorted])
    print(bgzip_result.returncode)
    print("Standard Output:", bgzip_result.stdout)
    print("Standard Error:", bgzip_result.stderr)
    index_result = subprocess.run(["tabix", "-p", "bed", f"{bed_file_sorted}.gz"])
    print(index_result.returncode)
    print("Standard Output:", index_result.stdout)
    print("Standard Error:", index_result.stderr)

    print(f"Using provided reference {reference_file}")
    maf_based_cmd = [
        "create_report",
        maf_file,
        "--fasta", reference_file,
        "--sequence", "1",
        "--begin", "2",
        "--end", "3",
        "--tracks", maf_file,
        "--track-config", "tracks_config.json",
        "--title", title,
        "--output", output_file,
        "--zero_based", "true" # add this flag if the bed file is zero based
    ]

    result = subprocess.run(maf_based_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print(result.returncode)
    print("Standard Output:", result.stdout)
    print("Standard Error:", result.stderr)
    print(f"IGV Reports visualization saved as '{output_file}'")


if __name__ == "__main__":
    bed_file_str = "data/test.maflite.maf"
    genome_str = "hg19"
    title_str = f"TEST"
    output_str = "test.html"
    create_igv_report(bed_file, genome, title, output)
