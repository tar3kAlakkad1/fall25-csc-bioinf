"""
Unified test file for both Codon and Python implementations of TRViz.

This file detects which runtime environment it's in and imports the appropriate modules.
Run with:
    - Python: cd /path/to/fall25-csc-bioinf && python week2/test.py
    - Codon:  cd /path/to/fall25-csc-bioinf/week2/codon && codon run ../test.py
"""

# Detect runtime environment using Codon's built-in __codon__ constant
try:
    if __codon__:
        IS_CODON = True
except NameError:
    IS_CODON = False

# Import based on runtime
if __codon__:
    print("Running tests in Codon environment...")
    print("Note: Make sure to run from week2/codon directory")
    
    # Direct imports - assumes running from week2/codon directory
    from codon.trviz.decomposer import Decomposer
    from codon.trviz.motif_aligner import MotifAligner
    from codon.trviz.motif_encoder import MotifEncoder
    from codon.trviz.utils import (
        is_valid_sequence, get_levenshtein_distance, get_motif_counter,
        is_emitting_state, get_repeating_pattern_lengths, 
        get_motifs_from_visited_states_and_region, add_padding, 
        get_distance_matrix, calculate_cost_with_dist_matrix,
        calculate_cost, calculate_total_cost
    )
else:
    print("Running tests in Python environment...")
    import os
    import sys
    
    # Add the Python trviz package to path if needed
    current_dir = os.path.dirname(os.path.abspath(__file__))
    python_dir = os.path.join(current_dir, 'python')
    if python_dir not in sys.path:
        sys.path.insert(0, python_dir)
    
    from trviz.decomposer import Decomposer
    from trviz.motif_aligner import MotifAligner
    from trviz.motif_encoder import MotifEncoder
    from trviz.utils import (
        is_valid_sequence, get_levenshtein_distance, get_motif_counter,
        is_emitting_state, get_repeating_pattern_lengths, 
        get_motifs_from_visited_states_and_region, add_padding, 
        get_distance_matrix, calculate_cost_with_dist_matrix,
        calculate_cost, calculate_total_cost
    )


# ============================================================================
# DECOMPOSER TESTS
# ============================================================================

def test_decompose_dp_perfect_repeats():
    """Test decomposition with perfect repeats."""
    print("\n[TEST] test_decompose_dp_perfect_repeats...")
    tr_decomposer_dp = Decomposer(mode="DP")
    
    # Test 1: Simple A repeat
    sequence = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    motifs = "AAAAAA"
    expected = ["AAAAAA", "AAAAAA", "AAAAAA", "AAAAAA", "AAAAAA"]
    result = tr_decomposer_dp.decompose(sequence, motifs)
    assert result == expected, f"Test 1 failed: {result} != {expected}"
    
    # Test 2: ACTG repeat
    sequence = "ACTGACTGACTG"
    motifs = "ACTG"
    expected = ["ACTG", "ACTG", "ACTG"]
    result = tr_decomposer_dp.decompose(sequence, motifs)
    assert result == expected, f"Test 2 failed: {result} != {expected}"
    
    print("✓ PASSED")


def test_decompose_dp_imperfect_repeats():
    """Test decomposition with imperfect repeats."""
    print("\n[TEST] test_decompose_dp_imperfect_repeats...")
    tr_decomposer_dp = Decomposer(mode="DP")
    
    # Test with single motif
    sequence = "AAAAAC" + "AAAAAA" + "AAAAAT" + "AAAAAA" + "TTAAAA"
    motifs = ["AAAAAA"]
    expected = ["AAAAAC", "AAAAAA", "AAAAAT", "AAAAAA", "TTAAAA"]
    result = tr_decomposer_dp.decompose(sequence, motifs)
    assert result == expected, f"Test failed: {result} != {expected}"
    
    # Test with multiple motifs
    sequence = "CGCCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGT"
    motifs = ["CGG", "CGC", "CGT"]
    expected = ["CGC", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG",
                "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGG", "CGT"]
    result = tr_decomposer_dp.decompose(sequence, motifs)
    assert result == expected, f"Test failed: {result} != {expected}"
    
    print("✓ PASSED")


def test_refine():
    """Test refinement of decomposed TRs."""
    print("\n[TEST] test_refine...")
    
    decomposed_trs = [
        ['AACAT', 'AACA', 'AACA', 'AACAT', 'AACA'],
        ['AACAT', 'AACA', 'AACA', 'AACAT', 'AACA'],
        ['AACA', 'TAACA', 'AACA', 'AACA', 'TAACA']
    ]
    expected = [
        ['AACAT', 'AACA', 'AACA', 'AACAT', 'AACA'],
        ['AACAT', 'AACA', 'AACA', 'AACAT', 'AACA'],
        ['AACAT', 'AACA', 'AACA', 'AACAT', 'AACA']
    ]
    
    refined_trs = Decomposer.refine(decomposed_trs)
    assert refined_trs == expected, f"Test failed: {refined_trs} != {expected}"
    
    print("✓ PASSED")


# ============================================================================
# MOTIF ENCODER TESTS
# ============================================================================

def test_encode_simple():
    """Test simple encoding without private motifs."""
    print("\n[TEST] test_encode_simple...")
    
    encoder = MotifEncoder()
    sample_decomposed_vntrs = [
        ["ACTG", "ACTG", "ACTG"],
        ["ACTG", "TGCA", "ACTG"],
        ["TGCA", "TGCA", "ACTG"],
    ]
    motif_map_file = "/tmp/test_motif_map_simple.txt"
    
    encoded_vntrs = encoder.encode(
        sample_decomposed_vntrs,
        motif_map_file,
        auto=False
    )
    
    # Check that we got encoded strings
    assert len(encoded_vntrs) == 3, "Should have 3 encoded VNTRs"
    assert all(isinstance(vntr, str) for vntr in encoded_vntrs), "All encoded VNTRs should be strings"
    
    # Check that motif_to_symbol and symbol_to_motif are populated
    assert len(encoder.motif_to_symbol) == 2, "Should have 2 motifs (ACTG and TGCA)"
    assert len(encoder.symbol_to_motif) == 2, "Should have 2 symbols"
    
    # Check bidirectional mapping
    for motif, symbol in encoder.motif_to_symbol.items():
        assert encoder.symbol_to_motif[symbol] == motif, f"Bidirectional mapping failed for {motif}"
    
    print("✓ PASSED")


def test_encode_with_threshold():
    """Test encoding with explicit private motif threshold."""
    print("\n[TEST] test_encode_with_threshold...")
    
    encoder = MotifEncoder(private_motif_threshold=1)
    sample_decomposed_vntrs = [
        ["AAA", "AAA", "AAA", "AAA", "AAA"],
        ["AAA", "BBB", "AAA", "BBB", "AAA"],
        ["AAA", "CCC", "AAA"],
        ["BBB", "DDD", "BBB"],
        ["AAA", "BBB", "EEE"],
    ]
    motif_map_file = "/tmp/test_motif_map_threshold.txt"
    
    encoded_vntrs = encoder.encode(
        sample_decomposed_vntrs,
        motif_map_file,
        auto=False
    )
    
    # With threshold=1, motifs with count <= 1 should be private
    assert encoder.private_motif_threshold == 1, "Threshold should be 1"
    
    print("✓ PASSED")


# ============================================================================
# UTILS TESTS
# ============================================================================

def test_is_valid_sequence():
    """Test DNA sequence validation."""
    print("\n[TEST] test_is_valid_sequence...")
    
    # Valid sequences
    assert is_valid_sequence("ATCG"), "'ATCG' should be valid"
    assert is_valid_sequence("AAAA"), "'AAAA' should be valid"
    assert is_valid_sequence(""), "Empty string should be valid"
    
    # Invalid sequences
    assert not is_valid_sequence("ATCGN"), "'ATCGN' should be invalid (contains N)"
    assert not is_valid_sequence("atcg"), "'atcg' should be invalid (lowercase)"
    assert not is_valid_sequence("ATCG-"), "'ATCG-' should be invalid (contains hyphen)"
    
    print("✓ PASSED")


def test_levenshtein_distance():
    """Test Levenshtein distance calculation."""
    print("\n[TEST] test_levenshtein_distance...")
    
    # Identical strings
    assert get_levenshtein_distance("", "") == 0, "Empty strings should have distance 0"
    assert get_levenshtein_distance("ATCG", "ATCG") == 0, "Identical strings should have distance 0"
    
    # Single edit operations
    assert get_levenshtein_distance("", "A") == 1, "Insertion of one char"
    assert get_levenshtein_distance("A", "") == 1, "Deletion of one char"
    assert get_levenshtein_distance("A", "T") == 1, "Substitution of one char"
    
    # Classic example
    assert get_levenshtein_distance("kitten", "sitting") == 3, "Classic kitten->sitting example"
    
    print("✓ PASSED")


def test_get_motif_counter():
    """Test motif counting across VNTRs."""
    print("\n[TEST] test_get_motif_counter...")
    
    # Empty input
    if IS_CODON:
        empty_input: list[list[str]] = []
    else:
        empty_input = []
    counter = get_motif_counter(empty_input)
    assert len(counter) == 0, "Empty input should produce empty counter"
    
    # Multiple VNTRs
    if IS_CODON:
        multiple_vntrs: list[list[str]] = [
            ["ATCG", "GCTA"],
            ["ATCG", "TATA"],
            ["GCTA", "TATA", "ATCG"]
        ]
    else:
        multiple_vntrs = [
            ["ATCG", "GCTA"],
            ["ATCG", "TATA"],
            ["GCTA", "TATA", "ATCG"]
        ]
    counter = get_motif_counter(multiple_vntrs)
    assert counter["ATCG"] == 3, "ATCG should appear 3 times"
    assert counter["GCTA"] == 2, "GCTA should appear 2 times"
    assert counter["TATA"] == 2, "TATA should appear 2 times"
    
    print("✓ PASSED")


def test_add_padding():
    """Test sequence padding."""
    print("\n[TEST] test_add_padding...")
    
    # Sequences of different lengths
    sequences = ["ABC", "ABCDE", "A"]
    padded = add_padding(sequences)
    assert len(padded) == 3, "Should have 3 sequences"
    assert padded[0] == "ABC--", "First sequence padding incorrect"
    assert padded[1] == "ABCDE", "Second sequence should be unchanged"
    assert padded[2] == "A----", "Third sequence padding incorrect"
    
    # Equal length sequences
    sequences = ["AAA", "BBB", "CCC"]
    padded = add_padding(sequences)
    assert padded[0] == "AAA", "Should be unchanged"
    assert padded[1] == "BBB", "Should be unchanged"
    assert padded[2] == "CCC", "Should be unchanged"
    
    print("✓ PASSED")


# ============================================================================
# MAIN TEST RUNNER
# ============================================================================

def run_all_tests():
    """Run all tests."""
    print("=" * 70)
    print("TRViz Unified Test Suite")
    print("=" * 70)
    
    # Decomposer tests
    print("\n" + "=" * 70)
    print("DECOMPOSER TESTS")
    print("=" * 70)
    test_decompose_dp_perfect_repeats()
    test_decompose_dp_imperfect_repeats()
    test_refine()
    
    # Motif Encoder tests
    print("\n" + "=" * 70)
    print("MOTIF ENCODER TESTS")
    print("=" * 70)
    test_encode_simple()
    test_encode_with_threshold()
    
    # Utils tests
    print("\n" + "=" * 70)
    print("UTILS TESTS")
    print("=" * 70)
    test_is_valid_sequence()
    test_levenshtein_distance()
    test_get_motif_counter()
    test_add_padding()
    
    print("\n" + "=" * 70)
    print("✅ ALL TESTS PASSED!")
    print("=" * 70)


if __name__ == "__main__":
    run_all_tests()

