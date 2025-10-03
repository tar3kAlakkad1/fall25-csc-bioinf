import pytest
from pathlib import Path

from trviz.motif_encoder import MotifEncoder
from trviz.utils import INDEX_TO_CHR, PRIVATE_MOTIF_LABEL


@pytest.fixture
def encoder():
    """Create a MotifEncoder instance."""
    return MotifEncoder()


@pytest.fixture
def sample_decomposed_vntrs_simple():
    """Simple decomposed VNTRs with a few motifs."""
    return [
        ["ACTG", "ACTG", "ACTG"],
        ["ACTG", "TGCA", "ACTG"],
        ["TGCA", "TGCA", "ACTG"],
    ]


@pytest.fixture
def sample_decomposed_vntrs_with_rare():
    """Decomposed VNTRs with frequent and rare motifs."""
    return [
        ["AAA", "AAA", "AAA", "AAA", "AAA"],  # AAA: 5
        ["AAA", "BBB", "AAA", "BBB", "AAA"],  # BBB: 2
        ["AAA", "CCC", "AAA"],                 # CCC: 1 (rare)
        ["BBB", "DDD", "BBB"],                 # DDD: 1 (rare)
        ["AAA", "BBB", "EEE"],                 # EEE: 1 (rare)
    ]


@pytest.fixture
def sample_decomposed_vntrs_many_motifs():
    """Decomposed VNTRs with many unique motifs (for testing ASCII limit)."""
    # Create enough unique motifs to exceed INDEX_TO_CHR length
    motifs = [f"M{i:03d}" for i in range(len(INDEX_TO_CHR) + 10)]
    return [[motif] for motif in motifs]


class TestMotifEncoderBasic:
    """Test basic encoding functionality."""
    
    def test_encode_simple_no_private(self, encoder, sample_decomposed_vntrs_simple, tmp_path):
        """Test encoding simple decomposed VNTRs without private motifs."""
        motif_map_file = tmp_path / "motif_map.txt"
        
        encoded_vntrs = encoder.encode(
            sample_decomposed_vntrs_simple,
            str(motif_map_file),
            auto=False
        )
        
        # Check that we got encoded strings
        assert len(encoded_vntrs) == 3
        assert all(isinstance(vntr, str) for vntr in encoded_vntrs)
        
        # Check that motif_to_symbol and symbol_to_motif are populated
        assert encoder.motif_to_symbol is not None
        assert encoder.symbol_to_motif is not None
        assert len(encoder.motif_to_symbol) == 2  # ACTG and TGCA
        assert len(encoder.symbol_to_motif) == 2
        
        # Check bidirectional mapping
        for motif, symbol in encoder.motif_to_symbol.items():
            assert encoder.symbol_to_motif[symbol] == motif
        
        # Check that the motif map file was created
        assert motif_map_file.exists()
        
        # Check score matrix is populated
        assert encoder.score_matrix is not None
        assert 'gap_open' in encoder.score_matrix
        assert 'gap_extension' in encoder.score_matrix
    
    def test_encode_with_auto_private_threshold(self, encoder, sample_decomposed_vntrs_with_rare, tmp_path):
        """Test encoding with automatic private motif threshold."""
        motif_map_file = tmp_path / "motif_map.txt"
        
        encoded_vntrs = encoder.encode(
            sample_decomposed_vntrs_with_rare,
            str(motif_map_file),
            auto=True
        )
        
        # Check that we got encoded strings
        assert len(encoded_vntrs) == 5
        
        # With auto and only 5 unique motifs (less than INDEX_TO_CHR length),
        # all motifs should fit without needing private labels
        # So threshold should be 0
        assert encoder.private_motif_threshold == 0
        
        # All motifs should get normal labels
        assert len(encoder.motif_to_symbol) == 5
    
    def test_encode_with_explicit_threshold(self, sample_decomposed_vntrs_with_rare, tmp_path):
        """Test encoding with explicit private motif threshold."""
        encoder = MotifEncoder(private_motif_threshold=1)
        motif_map_file = tmp_path / "motif_map.txt"
        
        encoded_vntrs = encoder.encode(
            sample_decomposed_vntrs_with_rare,
            str(motif_map_file),
            auto=False
        )
        
        # With threshold=1, motifs with count <= 1 should be private
        assert encoder.private_motif_threshold == 1
        
        # Check that rare motifs (count=1) are private
        rare_motifs = ["CCC", "DDD", "EEE"]
        for motif in rare_motifs:
            assert encoder.motif_to_symbol[motif] == PRIVATE_MOTIF_LABEL
    
    def test_encoded_strings_match_input_structure(self, encoder, sample_decomposed_vntrs_simple, tmp_path):
        """Test that encoded strings preserve the structure of input."""
        motif_map_file = tmp_path / "motif_map.txt"
        
        encoded_vntrs = encoder.encode(
            sample_decomposed_vntrs_simple,
            str(motif_map_file),
            auto=False
        )
        
        # Check that the length of encoded strings matches the length of decomposed VNTRs
        for i, (original, encoded) in enumerate(zip(sample_decomposed_vntrs_simple, encoded_vntrs)):
            assert len(encoded) == len(original), f"VNTR {i}: expected {len(original)} motifs, got {len(encoded)}"
    
    def test_encoded_strings_use_correct_symbols(self, encoder, sample_decomposed_vntrs_simple, tmp_path):
        """Test that encoded strings use the correct symbols for each motif."""
        motif_map_file = tmp_path / "motif_map.txt"
        
        encoded_vntrs = encoder.encode(
            sample_decomposed_vntrs_simple,
            str(motif_map_file),
            auto=False
        )
        
        # Verify that we can decode back to original motifs
        for i, (original, encoded) in enumerate(zip(sample_decomposed_vntrs_simple, encoded_vntrs)):
            decoded = [encoder.symbol_to_motif[symbol] for symbol in encoded]
            assert decoded == original, f"VNTR {i}: decoding mismatch"


class TestMotifMapFile:
    """Test motif map file generation."""
    
    def test_motif_map_file_contents(self, encoder, sample_decomposed_vntrs_simple, tmp_path):
        """Test that the motif map file has correct contents."""
        motif_map_file = tmp_path / "motif_map.txt"
        
        encoder.encode(
            sample_decomposed_vntrs_simple,
            str(motif_map_file),
            auto=False
        )
        
        # Read the motif map file
        with open(motif_map_file, 'r') as f:
            lines = f.readlines()
        
        # Should have one line per motif
        assert len(lines) == 2  # ACTG and TGCA
        
        # Each line should have: motif\tsymbol\tcount
        for line in lines:
            parts = line.strip().split('\t')
            assert len(parts) == 3
            motif, symbol, count = parts
            
            # Verify that the symbol matches the encoding
            assert encoder.motif_to_symbol[motif] == symbol
            
            # Verify that the count is correct
            assert int(count) == encoder.motif_counter[motif]
    
    def test_motif_map_file_ordered_by_frequency(self, encoder, sample_decomposed_vntrs_with_rare, tmp_path):
        """Test that motifs in the map file are ordered by frequency."""
        motif_map_file = tmp_path / "motif_map.txt"
        
        encoder.encode(
            sample_decomposed_vntrs_with_rare,
            str(motif_map_file),
            auto=False
        )
        
        # Read the motif map file
        with open(motif_map_file, 'r') as f:
            lines = f.readlines()
        
        # Extract counts
        counts = []
        for line in lines:
            parts = line.strip().split('\t')
            counts.append(int(parts[2]))
        
        # Verify that counts are in descending order
        assert counts == sorted(counts, reverse=True), "Motifs should be ordered by frequency"


class TestFindPrivateMotifThreshold:
    """Test the find_private_motif_threshold method."""
    
    def test_find_threshold_without_label_count(self, sample_decomposed_vntrs_with_rare):
        """Test finding private motif threshold without label_count constraint."""
        # With no label_count, all motifs should fit, so threshold should be 0
        threshold = MotifEncoder.find_private_motif_threshold(sample_decomposed_vntrs_with_rare)
        
        # We have 5 unique motifs, which is less than INDEX_TO_CHR length
        # So threshold should be 0
        assert threshold == 0
    
    def test_find_threshold_with_label_count(self, sample_decomposed_vntrs_with_rare):
        """Test finding private motif threshold with label_count constraint."""
        # With label_count=3, we can only encode 2 normal motifs (1 reserved for private)
        # So the threshold should be set to exclude the 3rd most common motif and below
        threshold = MotifEncoder.find_private_motif_threshold(
            sample_decomposed_vntrs_with_rare,
            label_count=3
        )
        
        # We have: AAA (10), BBB (5), CCC (1), DDD (1), EEE (1)
        # With label_count=3, maximum_label_count=2, so we keep first 2 motifs (AAA, BBB)
        # The threshold is the count of the 3rd motif (CCC), which is 1
        assert threshold == 1
    
    def test_find_threshold_with_small_label_count(self, sample_decomposed_vntrs_with_rare):
        """Test finding private motif threshold with very small label_count."""
        # With label_count=2, we can only encode 1 normal motif
        threshold = MotifEncoder.find_private_motif_threshold(
            sample_decomposed_vntrs_with_rare,
            label_count=2
        )
        
        # We have: AAA (10), BBB (5), CCC (1), DDD (1), EEE (1)
        # With label_count=2, maximum_label_count=1, so we keep first 1 motif (AAA)
        # The threshold is the count of the 2nd motif (BBB), which is 5
        # At index 1: 1+1=2 > 1, so threshold = 5
        assert threshold == 5


class TestErrorHandling:
    """Test error handling in MotifEncoder."""
    
    def test_too_many_motifs_raises_error(self, encoder, sample_decomposed_vntrs_many_motifs, tmp_path):
        """Test that encoding too many motifs raises ValueError."""
        motif_map_file = tmp_path / "motif_map.txt"
        
        with pytest.raises(ValueError, match="Too many unique motifs"):
            encoder.encode(
                sample_decomposed_vntrs_many_motifs,
                str(motif_map_file),
                auto=False
            )
    
    def test_too_many_normal_motifs_with_private_raises_error(self, tmp_path):
        """Test that too many normal motifs (even with private threshold) raises ValueError."""
        # Create enough motifs that even after removing private ones, we still exceed the limit
        num_motifs = len(INDEX_TO_CHR) + 5
        motifs = [f"M{i:03d}" for i in range(num_motifs)]
        # Each motif appears twice, so they're all "normal" with threshold=1
        decomposed_vntrs = [[motif, motif] for motif in motifs]
        
        encoder = MotifEncoder(private_motif_threshold=1)
        motif_map_file = tmp_path / "motif_map.txt"
        
        with pytest.raises(ValueError, match="Too many unique motifs"):
            encoder.encode(
                decomposed_vntrs,
                str(motif_map_file),
                auto=False
            )


class TestSymbolMapping:
    """Test symbol to motif mapping."""
    
    def test_symbol_to_motif_bidirectional(self, encoder, sample_decomposed_vntrs_simple, tmp_path):
        """Test that symbol_to_motif and motif_to_symbol are consistent."""
        motif_map_file = tmp_path / "motif_map.txt"
        
        encoder.encode(
            sample_decomposed_vntrs_simple,
            str(motif_map_file),
            auto=False
        )
        
        # Check bidirectional consistency
        for motif, symbol in encoder.motif_to_symbol.items():
            assert encoder.symbol_to_motif[symbol] == motif
        
        for symbol, motif in encoder.symbol_to_motif.items():
            assert encoder.motif_to_symbol[motif] == symbol
    
    def test_symbols_are_valid_ascii(self, encoder, sample_decomposed_vntrs_simple, tmp_path):
        """Test that all symbols are valid ASCII characters from INDEX_TO_CHR."""
        motif_map_file = tmp_path / "motif_map.txt"
        
        encoder.encode(
            sample_decomposed_vntrs_simple,
            str(motif_map_file),
            auto=False
        )
        
        # Check that all symbols are in INDEX_TO_CHR or PRIVATE_MOTIF_LABEL
        for symbol in encoder.symbol_to_motif.keys():
            assert symbol in INDEX_TO_CHR or symbol == PRIVATE_MOTIF_LABEL
    
    def test_private_motifs_all_map_to_same_symbol(self, sample_decomposed_vntrs_with_rare, tmp_path):
        """Test that all private motifs map to PRIVATE_MOTIF_LABEL."""
        encoder = MotifEncoder(private_motif_threshold=1)
        motif_map_file = tmp_path / "motif_map.txt"
        
        encoder.encode(
            sample_decomposed_vntrs_with_rare,
            str(motif_map_file),
            auto=False
        )
        
        # Find all private motifs (count <= 1)
        private_motifs = [motif for motif, count in encoder.motif_counter.items() if count <= 1]
        
        # Check that they all map to PRIVATE_MOTIF_LABEL
        for motif in private_motifs:
            assert encoder.motif_to_symbol[motif] == PRIVATE_MOTIF_LABEL


class TestScoreMatrix:
    """Test score matrix generation."""
    
    def test_score_matrix_created(self, encoder, sample_decomposed_vntrs_simple, tmp_path):
        """Test that score matrix is created after encoding."""
        motif_map_file = tmp_path / "motif_map.txt"
        
        encoder.encode(
            sample_decomposed_vntrs_simple,
            str(motif_map_file),
            auto=False
        )
        
        assert encoder.score_matrix is not None
        assert isinstance(encoder.score_matrix, dict)
    
    def test_score_matrix_has_gap_penalties(self, encoder, sample_decomposed_vntrs_simple, tmp_path):
        """Test that score matrix includes gap penalties."""
        motif_map_file = tmp_path / "motif_map.txt"
        
        encoder.encode(
            sample_decomposed_vntrs_simple,
            str(motif_map_file),
            auto=False
        )
        
        assert 'gap_open' in encoder.score_matrix
        assert 'gap_extension' in encoder.score_matrix
        assert isinstance(encoder.score_matrix['gap_open'], (int, float))
        assert isinstance(encoder.score_matrix['gap_extension'], (int, float))
    
    def test_score_matrix_has_entries_for_all_symbols(self, encoder, sample_decomposed_vntrs_simple, tmp_path):
        """Test that score matrix has entries for all symbols."""
        motif_map_file = tmp_path / "motif_map.txt"
        
        encoder.encode(
            sample_decomposed_vntrs_simple,
            str(motif_map_file),
            auto=False
        )
        
        # Check that all symbols have entries in the score matrix
        for symbol in encoder.symbol_to_motif.keys():
            assert symbol in encoder.score_matrix
            assert isinstance(encoder.score_matrix[symbol], dict)
            
            # Check that each symbol has scores against all other symbols
            for other_symbol in encoder.symbol_to_motif.keys():
                assert other_symbol in encoder.score_matrix[symbol]


class TestEncodeWithLabelCount:
    """Test encoding with label_count parameter."""
    
    def test_encode_with_label_count(self, encoder, sample_decomposed_vntrs_with_rare, tmp_path):
        """Test encoding with explicit label_count."""
        motif_map_file = tmp_path / "motif_map.txt"
        
        # Use label_count=3 to force private motif threshold
        encoded_vntrs = encoder.encode(
            sample_decomposed_vntrs_with_rare,
            str(motif_map_file),
            label_count=3,
            auto=False
        )
        
        # Check that threshold was set correctly
        # With label_count=3, we can keep 2 normal motifs (AAA, BBB)
        # The threshold is the count of the 3rd motif (CCC), which is 1
        assert encoder.private_motif_threshold == 1
        
        # Check that we only have 2 normal motifs + 1 private label = 3 symbols
        assert len(encoder.symbol_to_motif) == 3
        
        # Check that PRIVATE_MOTIF_LABEL is one of the symbols
        assert PRIVATE_MOTIF_LABEL in encoder.symbol_to_motif


class TestEdgeCases:
    """Test edge cases."""
    
    def test_encode_empty_vntrs(self, encoder, tmp_path):
        """Test encoding with empty VNTRs list."""
        motif_map_file = tmp_path / "motif_map.txt"
        
        encoded_vntrs = encoder.encode(
            [],
            str(motif_map_file),
            auto=False
        )
        
        assert encoded_vntrs == []
        assert encoder.motif_to_symbol == {}
        assert encoder.symbol_to_motif == {}
    
    def test_encode_single_motif(self, encoder, tmp_path):
        """Test encoding with single motif repeated."""
        decomposed_vntrs = [
            ["AAA", "AAA", "AAA"],
            ["AAA", "AAA"],
        ]
        motif_map_file = tmp_path / "motif_map.txt"
        
        encoded_vntrs = encoder.encode(
            decomposed_vntrs,
            str(motif_map_file),
            auto=False
        )
        
        # Should have exactly one symbol
        assert len(encoder.symbol_to_motif) == 1
        assert len(encoder.motif_to_symbol) == 1
        
        # All encoded VNTRs should use the same symbol
        symbol = list(encoder.symbol_to_motif.keys())[0]
        assert encoded_vntrs[0] == symbol * 3
        assert encoded_vntrs[1] == symbol * 2
    
    def test_encode_vntrs_with_single_occurrence_each(self, encoder, tmp_path):
        """Test encoding VNTRs where each motif occurs exactly once."""
        decomposed_vntrs = [
            ["AAA"],
            ["BBB"],
            ["CCC"],
        ]
        motif_map_file = tmp_path / "motif_map.txt"
        
        encoded_vntrs = encoder.encode(
            decomposed_vntrs,
            str(motif_map_file),
            auto=False
        )
        
        # Each motif should get its own symbol
        assert len(encoder.symbol_to_motif) == 3
        assert len(encoder.motif_to_symbol) == 3
        
        # Each encoded VNTR should be a single character
        assert all(len(vntr) == 1 for vntr in encoded_vntrs)
        
        # All encoded VNTRs should use different symbols
        assert len(set(encoded_vntrs)) == 3

