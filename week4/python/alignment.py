from typing import Final

import numpy as np


# Scoring scheme
MATCH_SCORE: Final[int] = 3
MISMATCH_PENALTY: Final[int] = -3
GAP_PENALTY: Final[int] = -2

# Affine gap penalties
GAP_OPEN_PENALTY: Final[int] = -5
GAP_EXTEND_PENALTY: Final[int] = -1


def global_align(seq1: str, seq2: str) -> int:
    """
    Compute the global alignment score (Needleman–Wunsch) for two sequences.

    - Uses a dynamic programming matrix implemented with numpy.ndarray.
    - Scoring: +3 match, -3 mismatch, -2 gap (linear gap penalty).
    - Returns only the final optimal alignment score (no traceback).

    Args:
        seq1: First sequence (string of characters)
        seq2: Second sequence (string of characters)

    Returns:
        The optimal global alignment score as an integer.
    """
    num_rows: int = len(seq1) + 1
    num_cols: int = len(seq2) + 1

    # DP matrix where dp[i, j] is the best score aligning seq1[:i] with seq2[:j]
    dp: np.ndarray = np.empty((num_rows, num_cols), dtype=int)

    # Initialize first row and column with cumulative gap penalties
    dp[0, 0] = 0
    if num_rows > 1:
        dp[1:, 0] = np.arange(1, num_rows) * GAP_PENALTY
    if num_cols > 1:
        dp[0, 1:] = np.arange(1, num_cols) * GAP_PENALTY

    # Fill the DP matrix
    for i in range(1, num_rows):
        for j in range(1, num_cols):
            is_match = seq1[i - 1] == seq2[j - 1]
            subst_score = MATCH_SCORE if is_match else MISMATCH_PENALTY

            score_diag = dp[i - 1, j - 1] + subst_score
            score_up = dp[i - 1, j] + GAP_PENALTY
            score_left = dp[i, j - 1] + GAP_PENALTY

            dp[i, j] = max(score_diag, score_up, score_left)

    return int(dp[num_rows - 1, num_cols - 1])


def local_align(seq1: str, seq2: str) -> int:
    """
    Compute the local alignment score (Smith–Waterman) for two sequences.

    - Uses a dynamic programming matrix implemented with numpy.ndarray.
    - Scoring: +3 match, -3 mismatch, -2 gap (linear gap penalty).
    - Each cell value is floored at 0 (no negative scores).
    - Returns the maximum score found anywhere in the matrix (no traceback).

    Args:
        seq1: First sequence (string of characters)
        seq2: Second sequence (string of characters)

    Returns:
        The optimal local alignment score as an integer.
    """
    num_rows: int = len(seq1) + 1
    num_cols: int = len(seq2) + 1

    # DP matrix where dp[i, j] is the best local score ending at (i, j)
    dp: np.ndarray = np.zeros((num_rows, num_cols), dtype=int)

    max_score: int = 0

    # Fill the DP matrix; first row/col remain 0 for Smith–Waterman
    for i in range(1, num_rows):
        for j in range(1, num_cols):
            is_match = seq1[i - 1] == seq2[j - 1]
            subst_score = MATCH_SCORE if is_match else MISMATCH_PENALTY

            score_diag = dp[i - 1, j - 1] + subst_score
            score_up = dp[i - 1, j] + GAP_PENALTY
            score_left = dp[i, j - 1] + GAP_PENALTY

            best_here = max(0, score_diag, score_up, score_left)
            dp[i, j] = best_here
            if best_here > max_score:
                max_score = best_here

    return int(max_score)


def semi_global_align(seq1: str, seq2: str) -> int:
    """
    Compute the semi-global (fitting) alignment score for two sequences.

    - No penalty for starting gaps: first row and first column initialized to 0.
    - Scoring: +3 match, -3 mismatch, -2 gap (linear gap penalty).
    - Returns the maximum score in the last row or last column.

    Args:
        seq1: First sequence (string of characters)
        seq2: Second sequence (string of characters)

    Returns:
        The optimal semi-global alignment score as an integer.
    """
    num_rows: int = len(seq1) + 1
    num_cols: int = len(seq2) + 1

    # DP matrix where dp[i, j] is the best score aligning seq1[:i] with seq2[:j]
    dp: np.ndarray = np.empty((num_rows, num_cols), dtype=int)

    # No penalties for starting gaps
    dp[0, 0] = 0
    if num_rows > 1:
        dp[1:, 0] = 0
    if num_cols > 1:
        dp[0, 1:] = 0

    # Fill the DP matrix (no floor at 0 unlike local alignment)
    for i in range(1, num_rows):
        for j in range(1, num_cols):
            is_match = seq1[i - 1] == seq2[j - 1]
            subst_score = MATCH_SCORE if is_match else MISMATCH_PENALTY

            score_diag = dp[i - 1, j - 1] + subst_score
            score_up = dp[i - 1, j] + GAP_PENALTY
            score_left = dp[i, j - 1] + GAP_PENALTY

            dp[i, j] = max(score_diag, score_up, score_left)

    # Best score ending anywhere in last row or last column
    last_row_max: int = int(dp[num_rows - 1, :].max())
    last_col_max: int = int(dp[:, num_cols - 1].max())
    return max(last_row_max, last_col_max)


def affine_global_align(seq1: str, seq2: str) -> int:
    """
    Compute the global alignment score with an affine gap penalty.

    - Uses three DP matrices: M (match/mismatch), Ix (gap in seq2), Iy (gap in seq1).
    - Scoring: +3 match, -3 mismatch.
    - Gap penalties: -5 open, -1 extend, with cost(open + (k-1)*extend) for k-length gap.
    - Returns only the final optimal alignment score (no traceback).

    Args:
        seq1: First sequence
        seq2: Second sequence

    Returns:
        The optimal global alignment score as an integer.
    """
    m: int = len(seq1)
    n: int = len(seq2)

    # Use a large negative sentinel for -inf in integer matrices
    NEG_INF: int = -(10**12)

    # DP matrices (int64 for safety)
    M: np.ndarray = np.empty((m + 1, n + 1), dtype=np.int64)
    Ix: np.ndarray = np.empty((m + 1, n + 1), dtype=np.int64)
    Iy: np.ndarray = np.empty((m + 1, n + 1), dtype=np.int64)

    M.fill(NEG_INF)
    Ix.fill(NEG_INF)
    Iy.fill(NEG_INF)

    # Initialization
    M[0, 0] = 0
    # Leading gaps in seq2 (vertical moves)
    for i in range(1, m + 1):
        Ix[i, 0] = GAP_OPEN_PENALTY + (i - 1) * GAP_EXTEND_PENALTY
        # M[i,0] and Iy[i,0] remain NEG_INF
    # Leading gaps in seq1 (horizontal moves)
    for j in range(1, n + 1):
        Iy[0, j] = GAP_OPEN_PENALTY + (j - 1) * GAP_EXTEND_PENALTY
        # M[0,j] and Ix[0,j] remain NEG_INF

    # Fill
    for i in range(1, m + 1):
        a = seq1[i - 1]
        for j in range(1, n + 1):
            b = seq2[j - 1]
            subst = MATCH_SCORE if a == b else MISMATCH_PENALTY

            # Match/mismatch state
            best_prev = M[i - 1, j - 1]
            if Ix[i - 1, j - 1] > best_prev:
                best_prev = Ix[i - 1, j - 1]
            if Iy[i - 1, j - 1] > best_prev:
                best_prev = Iy[i - 1, j - 1]
            M[i, j] = best_prev + subst

            # Gap in seq2 (vertical move): extend/delete in seq2
            open_gap_vert = M[i - 1, j] + GAP_OPEN_PENALTY
            extend_gap_vert = Ix[i - 1, j] + GAP_EXTEND_PENALTY
            Ix[i, j] = (
                open_gap_vert if open_gap_vert >= extend_gap_vert else extend_gap_vert
            )

            # Gap in seq1 (horizontal move): extend/insert in seq1
            open_gap_horiz = M[i, j - 1] + GAP_OPEN_PENALTY
            extend_gap_horiz = Iy[i, j - 1] + GAP_EXTEND_PENALTY
            Iy[i, j] = (
                open_gap_horiz
                if open_gap_horiz >= extend_gap_horiz
                else extend_gap_horiz
            )

    best_final = M[m, n]
    if Ix[m, n] > best_final:
        best_final = Ix[m, n]
    if Iy[m, n] > best_final:
        best_final = Iy[m, n]
    return int(best_final)


__all__ = ["global_align", "local_align", "semi_global_align", "affine_global_align"]
