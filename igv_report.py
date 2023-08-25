import subprocess

def create_igv_report(bed_file, genome, info_columns, title, output):
    cmd = [
        "create_report",
        bed_file,
        "--fasta", "/home/rswilson1/Documents/MSC_dissertation/MELT_example_data/hs37d5.fa",
        "--tracks", bed_file,
        "--track-config", "TracksConfig.json",
        #"--genome", genome,
        "--title", title,
        "--output", output
    ]

    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print(result.returncode)
    print("Standard Output:", result.stdout)
    print("Standard Error:", result.stderr)
    print(f"IGV Reports visualization saved as '{output}'")

if __name__ == "__main__":
    bed_file = "test/data/junctions/Introns.38.bed"
    genome = "hg38"
    info_columns = []
    title = f"TEST"
    output = "examples/example_junctions.html"

    create_igv_report(bed_file, genome, info_columns, title, output)
