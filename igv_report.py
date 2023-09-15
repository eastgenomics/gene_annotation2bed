import json
import subprocess


def create_igv_report(bed_file: str, maf_file: str,
                      genome: str, reference_file: str,
                      info_columns: list, title: str,
                      output_file: str) -> None:
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
    None.
    Prints the standard output and error of the subprocess.
    Creates an IGV report file.
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
            "url": f'{bed_file}.gz',
            "indexURL": f'{bed_file}.gz.tbi'
        },
    ]
    # Serializing json
    tracks_json = json.dumps(tracks_config, indent=4)

    # Writing to sample.json
    with open("tracks_config.json", "w") as outfile:
        outfile.write(tracks_json)
    sort_result = subprocess.run(["sort", "-k1,1", "-k2,2n", "-k3,3n", bed_file])
    print("Standard Output:", sort_result.stdout)
    print("Standard Error:", sort_result.stderr)
    bgzip_result = subprocess.run(["bgzip", bed_file])
    print(bgzip_result.returncode)
    print("Standard Output:", bgzip_result.stdout)
    print("Standard Error:", bgzip_result.stderr)
    index_result = subprocess.run(["tabix", f"{bed_file}.gz"])
    print(index_result.returncode)
    print("Standard Output:", index_result.stdout)
    print("Standard Error:", index_result.stderr)
    if reference_file:
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
            "--output", output_file
        ]

    elif genome:
        print(f"Using genome reference {genome}")
        maf_based_cmd = [
            "create_report",
            maf_file,
            "--genome", genome,
            "--tracks", maf_file,
            "--track-config", "tracks_config.json",
            "--title", title,
            "--output", output_file
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
