import pytest
import os
import tempfile
import shutil
from unittest.mock import patch, MagicMock, call
from io import StringIO

import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from trviz.motif_aligner import MotifAligner


class TestMotifAlignerDispatcher:
    """Test that the align method correctly dispatches to the requested tool."""

    def test_dispatcher_routes_to_mafft(self):
        """Test that align() routes to _align_motifs_with_mafft when tool='mafft'."""
        aligner = MotifAligner()
        sample_ids = ["sample1", "sample2"]
        encoded_vntrs = ["ACGT", "ACGG"]
        
        with patch.object(aligner, '_align_motifs_with_mafft', return_value=(sample_ids, encoded_vntrs)) as mock_mafft:
            result = aligner.align(sample_ids, encoded_vntrs, tool='mafft')
            mock_mafft.assert_called_once()
            assert result == (sample_ids, encoded_vntrs)

    def test_dispatcher_routes_to_muscle(self):
        """Test that align() routes to _align_motifs_with_muscle when tool='muscle'."""
        aligner = MotifAligner()
        sample_ids = ["sample1", "sample2"]
        encoded_vntrs = ["ACGT", "ACGG"]
        
        with patch.object(aligner, '_align_motifs_with_muscle', return_value=(sample_ids, encoded_vntrs)) as mock_muscle:
            result = aligner.align(sample_ids, encoded_vntrs, tool='muscle')
            mock_muscle.assert_called_once()
            assert result == (sample_ids, encoded_vntrs)

    def test_dispatcher_routes_to_clustalo(self):
        """Test that align() routes to _align_motifs_with_clustalo when tool='clustalo'."""
        aligner = MotifAligner()
        sample_ids = ["sample1", "sample2"]
        encoded_vntrs = ["ACGT", "ACGG"]
        
        with patch.object(aligner, '_align_motifs_with_clustalo', return_value=(sample_ids, encoded_vntrs)) as mock_clustalo:
            result = aligner.align(sample_ids, encoded_vntrs, tool='clustalo')
            mock_clustalo.assert_called_once()
            assert result == (sample_ids, encoded_vntrs)

    def test_dispatcher_routes_to_star(self):
        """Test that align() routes to _align_motifs_with_star when tool='star'."""
        aligner = MotifAligner()
        sample_ids = ["sample1", "sample2"]
        encoded_vntrs = ["ACGT", "ACGG"]
        
        with patch.object(aligner, '_align_motifs_with_star', return_value=(sample_ids, encoded_vntrs)) as mock_star:
            result = aligner.align(sample_ids, encoded_vntrs, tool='star')
            mock_star.assert_called_once()
            assert result == (sample_ids, encoded_vntrs)


class TestMafftAlignment:
    """Test the MAFFT alignment workflow."""

    def test_mafft_not_installed(self):
        """Test that ValueError is raised when MAFFT is not installed."""
        with patch('shutil.which', return_value=None):
            with pytest.raises(ValueError, match="MAFFT is not installed"):
                MotifAligner._align_motifs_with_mafft(
                    ["sample1"], ["ACGT"], None, None, "/tmp"
                )

    def test_mafft_success_without_score_matrix(self):
        """Test successful MAFFT alignment without score matrix."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Setup mock alignment output
            aln_output = f"{tmpdir}/alignment_output.fa"
            with open(aln_output, "w") as f:
                f.write(">sample1\nACGT\n>sample2\nACGG\n")
            
            sample_ids = ["sample1", "sample2"]
            encoded_vntrs = ["ACGT", "ACGG"]
            
            with patch('shutil.which', return_value='/usr/bin/mafft'):
                with patch('os.system', return_value=0) as mock_system:
                    result_ids, result_seqs = MotifAligner._align_motifs_with_mafft(
                        sample_ids, encoded_vntrs, None, None, tmpdir
                    )
                    
                    # Verify os.system was called
                    mock_system.assert_called_once()
                    call_args = mock_system.call_args[0][0]
                    assert 'mafft' in call_args
                    assert 'alignment_input.fa' in call_args
                    
                    # Verify results
                    assert result_ids == ["sample1", "sample2"]
                    assert result_seqs == ["ACGT", "ACGG"]

    def test_mafft_success_with_score_matrix(self):
        """Test successful MAFFT alignment with score matrix."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Setup mock alignment output
            vid = "test_vntr"
            aln_output = f"{tmpdir}/{vid}_alignment_output.fa"
            with open(aln_output, "w") as f:
                f.write(">sample1\nA-CGT\n>sample2\nACGG-\n")
            
            sample_ids = ["sample1", "sample2"]
            encoded_vntrs = ["ACGT", "ACGG"]
            score_matrix = {
                'A': {'A': 5, 'C': -1, 'G': -1, 'T': -1},
                'C': {'A': -1, 'C': 5, 'G': -1, 'T': -1},
                'G': {'A': -1, 'C': -1, 'G': 5, 'T': -1},
                'T': {'A': -1, 'C': -1, 'G': -1, 'T': 5},
                'gap_open': 10,
                'gap_extension': 1
            }
            
            with patch('shutil.which', return_value='/usr/bin/mafft'):
                with patch('subprocess.Popen') as mock_popen:
                    mock_process = MagicMock()
                    mock_process.returncode = 0
                    mock_process.wait = MagicMock()
                    mock_popen.return_value = mock_process
                    
                    result_ids, result_seqs = MotifAligner._align_motifs_with_mafft(
                        sample_ids, encoded_vntrs, vid, score_matrix, tmpdir
                    )
                    
                    # Verify subprocess.Popen was called
                    mock_popen.assert_called_once()
                    call_args = mock_popen.call_args[0][0]
                    assert 'mafft' in call_args
                    assert '--textmatrix' in call_args
                    assert '--op 10' in call_args
                    assert '--ep 1' in call_args
                    
                    # Verify results
                    assert result_ids == ["sample1", "sample2"]
                    assert result_seqs == ["A-CGT", "ACGG-"]

    def test_mafft_missing_output_file(self):
        """Test that FileNotFoundError is raised when alignment output is missing."""
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_ids = ["sample1", "sample2"]
            encoded_vntrs = ["ACGT", "ACGG"]
            
            with patch('shutil.which', return_value='/usr/bin/mafft'):
                with patch('os.system', return_value=0):
                    # Don't create the output file - should raise error
                    with pytest.raises(FileNotFoundError, match="motif alignment was not performed correctly"):
                        MotifAligner._align_motifs_with_mafft(
                            sample_ids, encoded_vntrs, None, None, tmpdir
                        )

    def test_mafft_preserve_order(self):
        """Test that preserve_order parameter affects MAFFT command."""
        with tempfile.TemporaryDirectory() as tmpdir:
            aln_output = f"{tmpdir}/alignment_output.fa"
            with open(aln_output, "w") as f:
                f.write(">sample1\nACGT\n")
            
            sample_ids = ["sample1"]
            encoded_vntrs = ["ACGT"]
            
            with patch('shutil.which', return_value='/usr/bin/mafft'):
                with patch('os.system', return_value=0) as mock_system:
                    # Test with preserve_order=False
                    MotifAligner._align_motifs_with_mafft(
                        sample_ids, encoded_vntrs, None, None, tmpdir, preserve_order=False
                    )
                    
                    call_args = mock_system.call_args[0][0]
                    assert '--reorder' in call_args


class TestLoadAlignedTrs:
    """Test the load_aligned_trs static method."""

    def test_load_aligned_trs_success(self):
        """Test successful loading of aligned VNTRs from FASTA file."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            f.write(">sample1\n")
            f.write("ACGT-TGA\n")
            f.write(">sample2\n")
            f.write("ACGG-TGA\n")
            f.write(">sample3\n")
            f.write("AC-GTTGA\n")
            temp_file = f.name
        
        try:
            sample_ids, aligned_seqs = MotifAligner.load_aligned_trs(temp_file)
            
            assert sample_ids == ["sample1", "sample2", "sample3"]
            assert aligned_seqs == ["ACGT-TGA", "ACGG-TGA", "AC-GTTGA"]
        finally:
            os.unlink(temp_file)

    def test_load_aligned_trs_empty_file(self):
        """Test that ValueError is raised when loading empty alignment file."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            # Write empty file
            temp_file = f.name
        
        try:
            with pytest.raises(ValueError, match="No aligned alleles found"):
                MotifAligner.load_aligned_trs(temp_file)
        finally:
            os.unlink(temp_file)

    def test_load_aligned_trs_single_sequence(self):
        """Test loading a single aligned sequence."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
            f.write(">only_sample\n")
            f.write("AAACCCTTTGGG\n")
            temp_file = f.name
        
        try:
            sample_ids, aligned_seqs = MotifAligner.load_aligned_trs(temp_file)
            
            assert sample_ids == ["only_sample"]
            assert aligned_seqs == ["AAACCCTTTGGG"]
        finally:
            os.unlink(temp_file)


class TestStarAlignment:
    """Test the center-star alignment method."""

    def test_star_alignment_basic(self):
        """Test basic center-star alignment with simple sequences."""
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_ids = ["center", "seq1", "seq2"]
            labeled_vntrs = ["ACGT", "ACGG", "ACTT"]
            
            # Mock the internal MAFFT calls
            with patch.object(MotifAligner, '_align_motifs_with_mafft') as mock_mafft:
                # First call: align center with seq1
                mock_mafft.side_effect = [
                    (["center", "seq1"], ["ACG-T", "ACGG-"]),
                    (["center", "seq2"], ["ACGT", "ACTT"])
                ]
                
                result_ids, result_seqs = MotifAligner._align_motifs_with_star(
                    sample_ids, labeled_vntrs, "test_vntr", None, tmpdir
                )
                
                # Verify MAFFT was called for pairwise alignments
                assert mock_mafft.call_count == 2
                
                # Verify we got all sample IDs back
                assert len(result_ids) == 3
                assert "center" in result_ids
                assert "seq1" in result_ids
                assert "seq2" in result_ids
                
                # Verify we got alignment results
                assert len(result_seqs) == 3

    def test_star_alignment_with_center_only(self):
        """Test star alignment with only center sequence (edge case)."""
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_ids = ["center", "seq1"]
            labeled_vntrs = ["ACGT", "ACGG"]
            
            with patch.object(MotifAligner, '_align_motifs_with_mafft') as mock_mafft:
                mock_mafft.return_value = (["center", "seq1"], ["ACGT", "ACGG"])
                
                result_ids, result_seqs = MotifAligner._align_motifs_with_star(
                    sample_ids, labeled_vntrs, "test_vntr", None, tmpdir
                )
                
                # Should have made one pairwise alignment call
                assert mock_mafft.call_count == 1
                assert len(result_ids) == 2
                assert len(result_seqs) == 2

    def test_star_alignment_with_gaps(self):
        """Test star alignment where gap insertion is needed."""
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_ids = ["center", "seq1", "seq2"]
            labeled_vntrs = ["ACGT", "A-GT", "ACGG"]
            
            with patch.object(MotifAligner, '_align_motifs_with_mafft') as mock_mafft:
                # Simulate alignments where gaps need to be merged
                mock_mafft.side_effect = [
                    (["center", "seq1"], ["AC-GT", "A--GT"]),
                    (["center", "seq2"], ["ACGT", "ACGG"])
                ]
                
                result_ids, result_seqs = MotifAligner._align_motifs_with_star(
                    sample_ids, labeled_vntrs, "test_vntr", None, tmpdir
                )
                
                # All sequences should be aligned
                assert len(result_ids) == 3
                assert len(result_seqs) == 3
                # All sequences should have the same length after alignment
                assert len(set(len(seq) for seq in result_seqs)) == 1


class TestScoreMatrixWriting:
    """Test score matrix file generation for MAFFT."""

    def test_score_matrix_file_creation(self):
        """Test that score matrix is correctly written for MAFFT."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create alignment output to avoid FileNotFoundError
            aln_output = f"{tmpdir}/alignment_output.fa"
            with open(aln_output, "w") as f:
                f.write(">sample1\nA\n")
            
            sample_ids = ["sample1"]
            encoded_vntrs = ["A"]
            score_matrix = {
                'A': {'A': 5, 'B': -1},
                'B': {'A': -1, 'B': 5},
                'gap_open': 10,
                'gap_extension': 1
            }
            
            with patch('shutil.which', return_value='/usr/bin/mafft'):
                with patch('subprocess.Popen') as mock_popen:
                    mock_process = MagicMock()
                    mock_process.returncode = 0
                    mock_process.wait = MagicMock()
                    mock_popen.return_value = mock_process
                    
                    MotifAligner._align_motifs_with_mafft(
                        sample_ids, encoded_vntrs, None, score_matrix, tmpdir
                    )
                    
                    # Verify the matrix file would have been created and used
                    call_args = mock_popen.call_args[0][0]
                    assert 'matrixfile.txt' in call_args
                    assert '--textmatrix' in call_args

    def test_score_matrix_file_removed_on_success(self):
        """Test that score matrix file is removed after successful alignment."""
        with tempfile.TemporaryDirectory() as tmpdir:
            aln_output = f"{tmpdir}/alignment_output.fa"
            with open(aln_output, "w") as f:
                f.write(">sample1\nA\n")
            
            sample_ids = ["sample1"]
            encoded_vntrs = ["A"]
            score_matrix = {
                'A': {'A': 5},
                'gap_open': 10,
                'gap_extension': 1
            }
            
            with patch('shutil.which', return_value='/usr/bin/mafft'):
                with patch('subprocess.Popen') as mock_popen:
                    mock_process = MagicMock()
                    mock_process.returncode = 0
                    mock_process.wait = MagicMock()
                    mock_popen.return_value = mock_process
                    
                    with patch('os.remove') as mock_remove:
                        MotifAligner._align_motifs_with_mafft(
                            sample_ids, encoded_vntrs, None, score_matrix, tmpdir
                        )
                        
                        # Verify os.remove was called on the matrix file
                        mock_remove.assert_called_once()
                        removed_file = mock_remove.call_args[0][0]
                        assert 'matrixfile.txt' in removed_file


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

