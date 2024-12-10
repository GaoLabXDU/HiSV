import argparse
import cooler
import pandas as pd
import os

def process_cooler_file(cool_file, out_dir, chrom1, chrom2):
    # Load cooler file and get chroms
    clr = cooler.Cooler(cool_file)
    chroms = clr.chromnames
    filtered_chroms = [chrom for chrom in chroms if chrom not in ['chrY', 'chrM']]
    # Ensure the chrom exists
    if chrom1 in filtered_chroms and chrom2 in filtered_chroms:
        # Ensure the output directory exists
        os.makedirs(out_dir, exist_ok=True)
        print(f"Processing {chrom1} and {chrom2} ...")
        hic_chr = clr.matrix(balance=False, as_pixels=True, join=True).fetch(chrom1, chrom2)
        hic_chr = hic_chr.iloc[:, [0, 1, 3, 4, 6]]
        names = ['chrom_1', 'pos_1', 'chrom_2', 'pos_2', 'count']
        hic_chr.columns = names
        # Generate output file name based on the chromosome
        outfile = os.path.join(out_dir, f"{chrom1}_{chrom2}.txt")
        # Save
        hic_chr.to_csv(outfile, sep='\t', index=None, header=True)
        print(f"{chrom1}_{chrom2} processed and saved to {outfile}.")
    else:
        print("Error chrom pairs.")
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert cooler file to bed files for each chromosome.")
    parser.add_argument("cool_file", type=str, help="Path to the .cool file.")
    parser.add_argument("out_dir", type=str, help="Directory to save the output .txt files.")
    parser.add_argument("chrom1", type=str, help="Chromosomal pairs.")
    parser.add_argument("chrom2", type=str, help="Chromosomal pairs.")
    args = parser.parse_args()
    # Process the cooler file
    process_cooler_file(args.cool_file, args.out_dir, args.chrom1, args.chrom2)