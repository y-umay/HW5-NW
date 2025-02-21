# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
     # Read human BRD2 sequence
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    
    # Read other species' BRD2 sequences
    species = {
        "Gallus gallus": read_fasta("./data/Gallus_gallus_BRD2.fa"),
        "Mus musculus": read_fasta("./data/Mus_musculus_BRD2.fa"),
        "Balaeniceps rex": read_fasta("./data/Balaeniceps_rex_BRD2.fa"),
        "Tursiops truncatus": read_fasta("./data/tursiops_truncatus_BRD2.fa")
    }
   
    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)

    alignment_scores = {}
    for name, (seq, _) in species.items():
        score, _, _ = nw.align(hs_seq, seq)
        alignment_scores[name] = score

    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    
    sorted_species = sorted(alignment_scores.items(), key=lambda x: x[1], reverse=True)
    print("Species ordered by similarity to Homo sapiens BRD2:")
    for name, score in sorted_species:
        print(f"{name}: {score}")


if __name__ == "__main__":
    main()
