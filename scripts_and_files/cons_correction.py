#!/usr/bin/env python3
# Usage
"""
cons_correction.py

Refines a consensus sequence generated from an aligned set of OTUs by
re-evaluating only positions where the consensus differs from a reference
sequence.

The correction follows the hierarchy:
    MISMATCHES
    OTU plurality > Reference-guided correction > OTU size
    GAPS
    OTU plurality > OTU size (default - can be changed with --ref-gap)

where:
    - OTU plurality is based on the most represented nucleotide (or gap)
      among valid OTUs at a given alignment position.
    - The reference sequence is used only when enabled and only if at least
      one OTU supports the reference state.
    - OTU size (number of reads supporting an OTU) is used as a confidence
      measure when neither plurality nor reference resolves the position.

Leading and trailing alignment padding are ignored by default. The
--extend-ends option allows correction outside the original consensus
boundaries, enabling recovery of terminal nucleotides supported by the OTUs.

Input:
    An aligned multi-FASTA containing:
        - one consensus sequence;
        - one reference sequence;
        - two or more OTUs with headers formatted as:

              >OTUx;size=y

Output:
    A FASTA file containing only the corrected consensus sequence.

Usage:
    python cons_correction.py \
        -i alignment.fasta \
        -r reference \
        -o corrected_consensus.fasta

Optional arguments:
    --extend-ends (default True)
        Allow correction outside the original consensus boundaries,
        potentially extending the consensus using OTU-supported terminal
        nucleotides.
    
    --ref-mismatch (default True)
        Use the reference sequence to resolve nucleotide mismatches.

    --ref-gap (default False)
        Use the reference sequence to resolve gap positions.  
"""

import argparse
from collections import defaultdict, Counter

# IUPAC ambiguity codes
IUPAC = {
    frozenset(["A", "G"]): "R",
    frozenset(["C", "T"]): "Y",
    frozenset(["G", "C"]): "S",
    frozenset(["A", "T"]): "W",
    frozenset(["G", "T"]): "K",
    frozenset(["A", "C"]): "M",
    frozenset(["A", "C", "G"]): "V",
    frozenset(["A", "C", "T"]): "H",
    frozenset(["A", "G", "T"]): "D",
    frozenset(["C", "G", "T"]): "B",
    frozenset(["A", "C", "G", "T"]): "N",
}

DNA_BASES = set(["A", "C", "G", "T"])

# Option parser
def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ("true", "t", "yes", "y", "1"):
        return True
    if v.lower() in ("false", "f", "no", "n", "0"):
        return False
    raise argparse.ArgumentTypeError("Boolean value expected.")

# FASTA reader
def read_fasta(path):
    seqs = {}
    header = None
    seq = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    seqs[header] = "".join(seq)
                header = line[1:]
                seq = []
            else:
                seq.append(line.upper())
        if header:
            seqs[header] = "".join(seq)

    return seqs

# Valid range detection (ignore padding gaps - namely the sites at the termini of the alignment that are not actually part of the sequence and that could introduce errors in the consensus generation)
def valid_range(seq):
    start = None
    end = None

    for i, c in enumerate(seq):
        if c != "-":
            start = i
            break

    for i in range(len(seq) - 1, -1, -1):
        if seq[i] != "-":
            end = i
            break

    if start is None:
        return None  # fully gaps (should not happen)
    return start, end

# OTU size parser
def parse_otu_size(header):
    try:
        return int(header.split("size=")[1])
    except Exception:
        raise ValueError(f"Cannot parse OTU size from header: {header}")

# Main function
def main():
    parser = argparse.ArgumentParser(
        description="Correct a consensus sequence using evidence from aligned OTUs and the reference sequence.",
        epilog=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-r", "--ref", required=True)
    parser.add_argument(
        "--ref-mismatch",
        type=str2bool,
        default=True,
        help="Use the reference to resolve nucleotide mismatches"
    )
    parser.add_argument(
        "--ref-gap",
        type=str2bool,
        default=False,
        help="Use the reference to resolve gap positions"
    )
    parser.add_argument(
        "--extend-ends",
        type=str2bool,
        default=True,
        help="Allow correction outside the original consensus boundaries"
    )
    args = parser.parse_args()

    seqs = read_fasta(args.input)

    ref_seq = None
    consensus_header = None
    consensus_seq = None
    otus = []

    for h, s in seqs.items():
        if h == args.ref:
            ref_seq = s
        elif h.startswith("OTU"):
            otus.append((h, s, parse_otu_size(h)))
        else:
            consensus_header = h
            consensus_seq = s

    if ref_seq is None:
        raise RuntimeError("Reference not found")
    if consensus_seq is None:
        raise RuntimeError("Consensus not found")

    aln_len = len(ref_seq)
    corrected = list(consensus_seq)

    # Precompute valid ranges
    consensus_range = valid_range(consensus_seq)
    otu_ranges = {}
    for h, s, size in otus:
        otu_ranges[h] = valid_range(s)

    # ----------------------------
    # Per-site correction - read the alignment column by column looking for suspicious sites
    # ----------------------------
    for i in range(aln_len):

        # Skip the column if it falls outside the consensus valid range,
        # unless end extension is enabled.
        if not args.extend_ends:
            if consensus_range is None or not (
                consensus_range[0] <= i <= consensus_range[1]
            ):
                continue

        R = ref_seq[i]
        C = corrected[i]

        # Not a suspicious site
        if C == R:
            continue

        # Collect valid OTUs at this site
        valid = []
        for h, s, size in otus:
            vr = otu_ranges[h]
            if vr and vr[0] <= i <= vr[1]:
                valid.append((s[i], size))

        # No valid OTUs -> keep original consensus
        if not valid:
            continue

        bases = [b for b, _ in valid]

        # Determine case type
        gap_case = (R == "-" or C == "-")
        mismatch_case = (R != "-" and C != "-")

        # ----------------------------
        # Step 1 - OTU plurality
        # ----------------------------
        counts = Counter(bases)

        max_count = max(counts.values())
        top_bases = [b for b, c in counts.items() if c == max_count]

        if len(top_bases) == 1:
            corrected[i] = top_bases[0]
            continue

        # ----------------------------
        # Step 2 - Reference anchor
        # ----------------------------
        if mismatch_case and args.ref_mismatch:
            if R in bases:
                corrected[i] = R
                continue

        if gap_case and args.ref_gap:
            if R in bases:
                corrected[i] = R
                continue

        # ----------------------------
        # Step 3 - Largest OTU fallback
        # ----------------------------
        size_by_base = defaultdict(int)

        for b, size in valid:
            size_by_base[b] = max(size_by_base[b], size)

        max_size = max(size_by_base.values())
        top_bases = [b for b, sz in size_by_base.items() if sz == max_size]

        if len(top_bases) == 1:
            corrected[i] = top_bases[0]
        else:
            # gap + base -> gap wins
            if "-" in top_bases:
                corrected[i] = "-"
            else:
                corrected[i] = IUPAC.get(
                    frozenset(b for b in top_bases if b in DNA_BASES),
                    "N",
                )

    # ----------------------------
    # Write output
    
    with open(args.output, "w") as out:
        out.write(f">{consensus_header}_corrected\n")
        out.write("".join(corrected) + "\n")


if __name__ == "__main__":
    main()
