import xml.etree.ElementTree as ET
import re
import argparse

# Namespace used by BLAST XML
NS = {"ns": "http://www.ncbi.nlm.nih.gov"}

# --------------------------------------------------
# GAP-BASED TRIMMING (coordinate-based)
# --------------------------------------------------

def gap_trim_coordinates(seq, window=200, min_gap=5):
    """
    Identify left and right trimming coordinates based on long gap stretches.
    Returns (left_cut, right_cut).
    """
    left_cut = 0
    right_cut = len(seq)

    # Left side: rightmost long gap in first window
    left_region = seq[:window]
    matches = list(re.finditer(rf"-{{{min_gap + 1},}}", left_region))
    if matches:
        left_cut = matches[-1].end()

    # Right side: leftmost long gap in last window
    right_region = seq[-window:]
    matches = list(re.finditer(rf"-{{{min_gap + 1},}}", right_region))
    if matches:
        right_cut = len(seq) - window + matches[0].start()

    return left_cut, right_cut


# --------------------------------------------------
# MIDLINE-BASED TRIMMING (unchanged logic)
# --------------------------------------------------

def trim_left_by_midline(qseq, midline, window=50):
    search_region = midline[:window]
    bad_positions = [i for i, c in enumerate(search_region) if c == " "]

    if bad_positions:
        cut = bad_positions[-1] + 1
        return qseq[cut:], midline[cut:]

    return qseq, midline


def trim_right_by_midline(qseq, midline, window=50):
    search_region = midline[-window:]
    bad_positions = [i for i, c in enumerate(search_region) if c == " "]

    if bad_positions:
        cut = len(midline) - window + bad_positions[0]
        return qseq[:cut], midline[:cut]

    return qseq, midline


# --------------------------------------------------
# MAIN
# --------------------------------------------------

def main(input_xml, output_fasta):
    tree = ET.parse(input_xml)
    root = tree.getroot()

    # First Search only
    search = root.find(".//ns:Search", NS)
    header = search.findtext("ns:query-title", namespaces=NS)

    # First HSP only
    hsp = search.find(".//ns:Hsp", NS)

    qseq = hsp.findtext("ns:qseq", namespaces=NS)
    midline = hsp.findtext("ns:midline", namespaces=NS)

    # --------------------------------------------------
    # Stage 1: GAP-based trimming (apply to both)
    # --------------------------------------------------
    left_cut, right_cut = gap_trim_coordinates(qseq)

    qseq = qseq[left_cut:right_cut]
    midline = midline[left_cut:right_cut]

    # --------------------------------------------------
    # Stage 2: MIDLINE-based trimming
    # --------------------------------------------------
    qseq, midline = trim_left_by_midline(qseq, midline)
    qseq, midline = trim_right_by_midline(qseq, midline)

    # --------------------------------------------------
    # Final cleanup
    # --------------------------------------------------
    qseq = qseq.replace("-", "")

    # --------------------------------------------------
    # Output
    # --------------------------------------------------
    with open(output_fasta, "w") as out:
        out.write(f">{header}\n")
        out.write(qseq + "\n")


# --------------------------------------------------
# CLI
# --------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Trim BLAST XML query sequence using gap-based trimming followed by midline mismatch trimming (first HSP only)"
    )
    parser.add_argument(
        "-i", "--input", required=True, help="Input BLAST XML file (outfmt 16)"
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Output trimmed FASTA file"
    )

    args = parser.parse_args()
    main(args.input, args.output)
