import os, glob
import argparse
import csv
import re
import sys


def main():

    args = get_args()
    filenames = get_input(args.i)

    min = args.min
    plddt = args.plddt
    if args.o:
        if os.path.exists(args.o):
            sys.exit("Given output file already exists!")
        elif args.o.endswith(".fasta") or args.o.endswith(".tsv"):
            output_file = open(args.o, "w", newline="")
            if args.o.endswith(".tsv"):
                header = 0
        else:
            sys.exit("Invalid file format!")
    else:
        output_file = None

    for filename in filenames:
        chains = {}
        name = filename.split("\\")[-1].split(".cif")[0]
        with open(filename, "r") as file:
            for line in file:
                if line.startswith("ATOM"):
                    match = re.match(
                        r"^ATOM.{8}CA.{4}(\w{3}).(\w).{3}(\d{1,4}).{30,34}(\d\d\.\d\d)",
                        line,
                    )
                    if match:
                        chain_id = match.group(2)
                        if chain_id not in chains:
                            chains[chain_id] = {"res": [], "pos": [], "plddt": []}
                        chains[chain_id]["res"].append(
                            convert_three_to_one_letter_aa(match.group(1))
                        )
                        chains[chain_id]["pos"].append(match.group(3))
                        chains[chain_id]["plddt"].append(match.group(4))
            seq = {}
            residues = {}
            positions = {}
            ranges = {}
            length = {}
            for chain_id, value in chains.items():
                seq[chain_id] = "".join(
                    res.lower() if float(score) <= plddt else res
                    for res, score in zip(value["res"], value["plddt"])
                )
                temp = "".join(
                    res.lower() if float(score) <= plddt else " "
                    for res, score in zip(value["res"], value["plddt"])
                )
                temp = " ".join(l if len(l) >= min else " " for l in temp.split())
                residues[chain_id] = ",".join(temp.split())
                temp = " ".join(
                    pos if float(plddt) <= 50 else ""
                    for pos, plddt in zip(value["pos"], value["plddt"])
                )
                temp = ",".join(temp.split())
                positions[chain_id] = temp
                ranges[chain_id] = convert_to_ranges(temp, min)[0]
                length[chain_id] = convert_to_ranges(temp, min)[1]

        if output_file:
            if args.o.endswith(".fasta"):
                for chain in chains:
                    output_file.write(">" + name + "_" + chain + "\n")
                    output_file.write(seq[chain] + "\n")
            if args.o.endswith(".tsv"):
                writer = csv.writer(output_file, delimiter="\t")
                if header == 0:
                    filednames = ["ID", "Chain", "residue", "range", "length"]
                    writer.writerow(filednames)
                    header = 1
                writer.writerows(
                    zip(
                        [name] * len(chains.keys()),
                        chains.keys(),
                        residues.values(),
                        ranges.values(),
                        length.values(),
                    )
                )
        else:
            print(f"Name:{name}")
            print(f"Sequence:{seq}")
            print(f"Residues:{residues}")
            print(f"Positions:{positions}")
            print(f"Ranges:{ranges}")
            print(f"Length:{length}")

    if output_file:
        output_file.close()


def convert_three_to_one_letter_aa(seq):
    """Returns the one letter codes of amino acid residues. Takes a string of three letter amino acid residues"""
    amino_acid_dict = {
        "ALA": "A",
        "ARG": "R",
        "ASN": "N",
        "ASP": "D",
        "CYS": "C",
        "GLU": "E",
        "GLN": "Q",
        "GLY": "G",
        "HIS": "H",
        "ILE": "I",
        "LEU": "L",
        "LYS": "K",
        "MET": "M",
        "PHE": "F",
        "PRO": "P",
        "SER": "S",
        "THR": "T",
        "TRP": "W",
        "TYR": "Y",
        "VAL": "V",
        "ASX": "B",
        "GLX": "Z",
        "XLE": "J",
        "XAA": "X",
        "SEC": "U",
        "PYL": "O",
    }
    if len(seq) % 3 != 0:
        raise ValueError("Sequence length should be a multiple of three")
    else:
        return "".join(amino_acid_dict[seq[i : i + 3]] for i in range(0, len(seq), 3))


def b_less_than(plddt, lst, value):
    return [lst[i] for i in range(len(plddt)) if plddt[i] <= str(value)]


def convert_to_ranges(values, min=1):
    """Returns a range of numbers and their associated length.
    Takes a string of numbers,internally seperated by commas, and an optional argument about the minimum length of the range
    """
    values = values.split(",")
    numbers = list(map(int, values))
    numbers.sort()
    ranges = []
    start = end = numbers[0]
    for num in numbers[1:]:
        if num != end + 1:
            ranges.append((start, end))
            start = end = num
        else:
            end = num
    ranges.append((start, end))

    range_str = " ".join(
        f"{s}-{e}" if e != s and e - s + 1 >= min else f"{s}" for s, e in ranges
    )
    range_str = ",".join(range_str.split())

    len_str = " ".join(str(e - s + 1) if e - s + 1 >= min else "" for s, e in ranges)
    len_str = ",".join(len_str.split())
    return range_str, len_str


def get_args():
    """Parses command line arguments that user can give"""
    parser = argparse.ArgumentParser(
        description="Arguments to determine amino acid residues having pLDDT values less than a specified one"
    )
    parser.add_argument(
        "-i",
        required=False,
        type=str,
        help="Enter directory to look for AlphaFold model cif files",
    )
    parser.add_argument(
        "-plddt",
        default=50,
        type=float,
        help="Look for residues having pLDDT <= given_value",
    )
    parser.add_argument(
        "-o",
        required=False,
        type=str,
        help="Specify a '.tsv' file or '.fasta' file to save the output, otherwise prints the result",
    )
    parser.add_argument(
        "-min",
        default=1,
        type=int,
        help="Minimum continous length over which the plddt value should apply",
    )
    args = parser.parse_args()
    return args


def get_input(input_folder):
    """Takes the path to a directory to look for input files,otherwise takes input from user of the filepath"""
    if input_folder:
        path = input_folder
        for filename in glob.glob(os.path.join(path, "*.cif")):
            yield filename
    else:
        filename = input("Enter filename:")
        yield filename


if __name__ == "__main__":
    main()
