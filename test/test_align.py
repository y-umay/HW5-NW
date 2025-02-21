# Importing Dependencies
import pytest
import os
from align import NeedlemanWunsch, read_fasta
import numpy as np

'''
To run python -m pytest -v test/* 
'''

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    """
    Steps:
    1. Verify the existence of required files (substitution matrix and FASTA sequences).
    2. Perform sequence alignment using the custom Needleman-Wunsch implementation.
    3. Validate the alignment score and sequences against known expected values.
    """
    # Validate substitution matrix existence
    sub_matrix_path = "./substitution_matrices/BLOSUM62.mat"
    if not os.path.isfile(sub_matrix_path):
        pytest.fail(f"Substitution matrix not found at {sub_matrix_path}. Please ensure the file is present.")

    # Validate existence of FASTA files
    for fasta_file in ["./data/test_seq1.fa", "./data/test_seq2.fa"]:
        if not os.path.isfile(fasta_file):
            pytest.fail(f"Required FASTA file not found: {fasta_file}")

    # Read sequences from FASTA files
    seq1, _ = read_fasta("./data/test_seq1.fa")  
    seq2, _ = read_fasta("./data/test_seq2.fa") 

    # Perform custom alignment
    nw = NeedlemanWunsch(sub_matrix_path, gap_open=-10, gap_extend=-1)
    custom_score, aligned_seq1, aligned_seq2 = nw.align(seq1, seq2)

    # Expected results (manually verified or provided by assignment specs)
    expected_score = 14.0  # Replace with correct expected alignment score
    expected_aligned_seq1 = "MYQR"  # Replace with correct expected alignment
    expected_aligned_seq2 = "M-QR"  # Replace with correct expected alignment

    # Assertions to validate alignment correctness
    assert round(custom_score, 2) == expected_score, (
        f"Alignment score mismatch: expected {expected_score}, got {custom_score}"
    )
    assert aligned_seq1 == expected_aligned_seq1, (
        f"Aligned sequence 1 mismatch: expected {expected_aligned_seq1}, got {aligned_seq1}"
    )
    assert aligned_seq2 == expected_aligned_seq2, (
        f"Aligned sequence 2 mismatch: expected {expected_aligned_seq2}, got {aligned_seq2}"
    )

def test_nw_backtrace():
    """
    Unit test for NW backtrace using test_seq3.fa (MAVHQLIRRP) and test_seq4.fa (MQLIRHP).

    Steps:
    1. Verify the presence of substitution matrix and FASTA files.
    2. Perform alignment with the custom Needleman-Wunsch implementation.
    3. If Biopython is available, validate the alignment score using Bio.Align.PairwiseAligner.
    """
    # Validate substitution matrix existence
    sub_matrix_path = "./substitution_matrices/BLOSUM62.mat"
    if not os.path.isfile(sub_matrix_path):
        pytest.fail(f"Substitution matrix not found at {sub_matrix_path}. Please ensure the file is present.")
    
    # Validate existence of FASTA files
    for fasta_file in ["./data/test_seq3.fa", "./data/test_seq4.fa"]:
        if not os.path.isfile(fasta_file):
            pytest.fail(f"Required FASTA file not found: {fasta_file}")

    # Read sequences from FASTA files
    seq3, _ = read_fasta("./data/test_seq3.fa")  
    seq4, _ = read_fasta("./data/test_seq4.fa")  

    # Perform custom alignment
    nw = NeedlemanWunsch(sub_matrix_path, gap_open=-10, gap_extend=-1)
    custom_score, aligned_seq3, aligned_seq4 = nw.align(seq3, seq4)

    # Expected results (provided by assignment specs or manually verified)
    expected_score = 27.0
    expected_aligned_seq3 = "MAVHQLIRRP"
    expected_aligned_seq4 = "M---QLIRHP"

    # Assertions to validate alignment correctness
    assert round(custom_score, 2) == expected_score, (
        f"Alignment score mismatch: expected {expected_score}, got {custom_score}"
    )
    assert aligned_seq3 == expected_aligned_seq3, (
        f"Aligned sequence 3 mismatch: expected {expected_aligned_seq3}, got {aligned_seq3}"
    )
    assert aligned_seq4 == expected_aligned_seq4, (
        f"Aligned sequence 4 mismatch: expected {expected_aligned_seq4}, got {aligned_seq4}"
    )
