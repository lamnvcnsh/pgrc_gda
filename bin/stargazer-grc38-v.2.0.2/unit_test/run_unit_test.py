import subprocess
import os

def run_stargazer(test_dir, gene):
    subprocess.call([
        "stargazer",
        "-o", f"getrm137-{gene}-vdr",
        "-d", "wgs",
        "-t", gene,
        "-c", "vdr",
        "-i", f"{test_dir}/getrm137_{gene}.vcf",
        "-g", f"{test_dir}/getrm137_{gene}.gdf",
    ])

def run_test(gene, gt_file1, gt_file2):
    subprocess.call([
        "pypgx", "compare-stargazer-calls",
        "-r", gt_file1,
        "-t", gt_file2,
        "-o", f"{gene}.txt"
    ])

def main():
    test_dir = os.path.dirname(os.path.realpath(__file__))
    genes = ["cyp2d6", "cyp2c19", "ugt1a1"]
    for gene in genes:
        run_stargazer(test_dir, gene)
        run_test(gene,
                 f"{test_dir}/getrm137_{gene}.txt",
                 f"getrm137-{gene}-vdr/genotype-calls.tsv")

if __name__ == "__main__":
    main()
