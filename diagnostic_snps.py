from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from collections import defaultdict
from copy import copy
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description="Outputs all diagnostic SNPs of an aligmnent in terms of positions in a reference sequence. Does strict SNPs (that only occur in one sequence) or all SNPs. Ignores uncalled bases.")

    parser.add_argument(
        "-i", "--input", type=str, help="Alignment file path."
    )

    parser.add_argument(
        "-ref", "--reference", type=str, required=True, help="Reference sequence ID as given in alignment file."
    )

    parser.add_argument("-s", "--stop", type=int, required=False, help="Stop at position X of the reference sequence.")

    parser.add_argument(
        "--strict",
        action="store_true",
        required=False,
        help="If invoked, this will restrict diagnostic SNPs to SNPs that are present in exactly one sequence.",
    )

    parser.add_argument(
        "--ignore", nargs="+", help="ID(s) of the sequence(s) you want to ignore."
    )

    args = parser.parse_args()

    # The next 6 lines are just for input handling
    # If neither -i nor data from stdin is given throw an error
    if args.input is None and sys.stdin.isatty():
        parser.error("Either -i/--input or stdin must be provided")

    if args.input:
        alignment = AlignIO.read(args.input, "fasta")
    else:
        alignment = AlignIO.read(sys.stdin, "fasta")

    og_alignment = copy(alignment) #needed for find_reference alignment if ignored individual is also ref
    # Filtering
    if args.ignore:
        filtered_records = [
            record for record in alignment if record.id not in args.ignore
        ]
        og_len = len(alignment)
        alignment = MultipleSeqAlignment(filtered_records)
        print(f"Removed {og_len - len(alignment)} sequence(s).", file=sys.stderr)


    # Get the number of sequences and length of the alignment
    alignment_length = alignment.get_alignment_length()

    # Initialize a dictionary to store the counts of each nucleotide at each position
    position_counts = defaultdict(lambda: defaultdict(int))

    # Count the nucleotides at each position
    for record in alignment:
        for i in range(alignment_length):
            nucleotide = record.seq[i]
            if nucleotide in ["n", "N"]:
                continue
            position_counts[i][nucleotide] += 1

    # Identify variable positions
    variable_positions = [
        pos for pos, counts in position_counts.items() if len(counts) > 1
    ]

    # Function to determine if a position is diagnostic
    def is_diagnostic(pos):
        nucleotide_groups = defaultdict(list)
        for record in alignment:
            nucleotide_groups[record.seq[pos]].append(record.id)
        # Check if any nucleotide is unique to a single sequence or group of sequences
        for nucleotide, seq_ids in nucleotide_groups.items():
            if len(seq_ids) == 1:
                return True
        return False

    # Find diagnostic positions
    if args.strict:
        diagnostic_positions = [pos for pos in variable_positions if is_diagnostic(pos)]
    else: # all SNPS are considered diagnostic SNPs
        diagnostic_positions = variable_positions

    diagnostic_positions.sort()

    def find_reference_alignment(ref_name: str, alignment) -> str:
        """Extracts the named reference from the MSA."""

        # for record in alignment:
        for record in alignment:#SeqIO.parse(alignment, "fasta"):
            if record.id == ref_name:

                if len(record.seq) == 0:
                    raise ValueError(
                        f"No sequence found for Reference ID {ref_name} in alignment."
                    )
                return record.seq  # upper()

        raise ValueError(
            f"Reference ID {ref_name} not found in alignment."
        )

    aligned_ref = find_reference_alignment(args.reference, og_alignment)

    def aln_to_ref(aligned_ref):
        aln_to_ref = {}
        index = 0
        for i, nuc in enumerate(aligned_ref):
            if nuc != "-":
                index += 1
            aln_to_ref[i] = index
        return aln_to_ref

    aln_to_reference = aln_to_ref(aligned_ref)

    def print_nucleotides(diagnostic_positions, alignment):
        print("pos", end="\t")

        for record in alignment:
            print(record.id, end="\t")
        print()
        
        for diagnostic_position in diagnostic_positions:
            aln_pos = aln_to_reference[diagnostic_position]
            if args.stop:
                if aln_pos > args.stop:
                    exit()
            print(aln_pos, end="\t")
            for record in alignment:
                print(record.seq[diagnostic_position].upper(), end="\t")
            print()

    print_nucleotides(diagnostic_positions, alignment)


if __name__ == "__main__":
    main()
