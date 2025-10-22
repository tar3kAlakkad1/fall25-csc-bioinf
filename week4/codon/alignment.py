import numpy as np
from numpy import ndarray

# Scoring scheme (linear gap)
MATCH_SCORE: int = 3
MISMATCH_PENALTY: int = -3
GAP_PENALTY: int = -2

# Affine gap penalties
GAP_OPEN_PENALTY: int = -5
GAP_EXTEND_PENALTY: int = -1


def global_align(seq1: str, seq2: str) -> int:
    """
    Compute the global alignment score (Needleman–Wunsch) for two sequences.

    Uses a 2D DP matrix of ints with linear gap penalty.
    Returns only the final optimal alignment score (no traceback).
    """
    num_rows: int = len(seq1) + 1
    num_cols: int = len(seq2) + 1

    # DP matrix where dp[i, j] is the best score aligning seq1[:i] with seq2[:j]
    dp: ndarray[int, 2] = np.empty((num_rows, num_cols), dtype=int)

    # Initialize first row and column with cumulative gap penalties
    dp[0, 0] = 0
    # First column
    i: int = 1
    while i < num_rows:
        dp[i, 0] = dp[i - 1, 0] + GAP_PENALTY
        i += 1
    # First row
    j: int = 1
    while j < num_cols:
        dp[0, j] = dp[0, j - 1] + GAP_PENALTY
        j += 1

    # Fill the DP matrix
    i = 1
    while i < num_rows:
        a = seq1[i - 1]
        j = 1
        while j < num_cols:
            b = seq2[j - 1]
            subst_score: int = MATCH_SCORE if a == b else MISMATCH_PENALTY

            score_diag: int = dp[i - 1, j - 1] + subst_score
            score_up: int = dp[i - 1, j] + GAP_PENALTY
            score_left: int = dp[i, j - 1] + GAP_PENALTY

            # max of three ints
            best: int = score_diag if score_diag >= score_up else score_up
            best = best if best >= score_left else score_left
            dp[i, j] = best
            j += 1
        i += 1

    return int(dp[num_rows - 1, num_cols - 1])


def local_align(seq1: str, seq2: str) -> int:
    """
    Compute the local alignment score (Smith–Waterman) for two sequences.

    Uses a 2D DP matrix of ints with linear gap penalty.
    Each cell value is floored at 0 (no negative scores).
    Returns the maximum score found anywhere in the matrix (no traceback).
    """
    num_rows: int = len(seq1) + 1
    num_cols: int = len(seq2) + 1

    # DP matrix where dp[i, j] is the best local score ending at (i, j)
    dp: ndarray[int, 2] = np.zeros((num_rows, num_cols), dtype=int)

    max_score: int = 0

    i: int = 1
    while i < num_rows:
        a = seq1[i - 1]
        j: int = 1
        while j < num_cols:
            b = seq2[j - 1]
            subst_score: int = MATCH_SCORE if a == b else MISMATCH_PENALTY

            score_diag: int = dp[i - 1, j - 1] + subst_score
            score_up: int = dp[i - 1, j] + GAP_PENALTY
            score_left: int = dp[i, j - 1] + GAP_PENALTY

            # best_here = max(0, score_diag, score_up, score_left)
            best_here: int = 0
            if score_diag > best_here:
                best_here = score_diag
            if score_up > best_here:
                best_here = score_up
            if score_left > best_here:
                best_here = score_left

            dp[i, j] = best_here
            if best_here > max_score:
                max_score = best_here
            j += 1
        i += 1

    return int(max_score)


def semi_global_align(seq1: str, seq2: str) -> int:
    """
    Compute the semi-global (fitting) alignment score for two sequences.

    No penalty for starting gaps: first row and first column initialized to 0.
    Returns the maximum score in the last row or last column.
    """
    num_rows: int = len(seq1) + 1
    num_cols: int = len(seq2) + 1

    dp: ndarray[int, 2] = np.empty((num_rows, num_cols), dtype=int)

    # Initialize first row/column with 0s (no starting gap penalties)
    dp[0, 0] = 0
    i: int = 1
    while i < num_rows:
        dp[i, 0] = 0
        i += 1
    j: int = 1
    while j < num_cols:
        dp[0, j] = 0
        j += 1

    # Fill the DP matrix
    i = 1
    while i < num_rows:
        a = seq1[i - 1]
        j = 1
        while j < num_cols:
            b = seq2[j - 1]
            subst_score: int = MATCH_SCORE if a == b else MISMATCH_PENALTY

            score_diag: int = dp[i - 1, j - 1] + subst_score
            score_up: int = dp[i - 1, j] + GAP_PENALTY
            score_left: int = dp[i, j - 1] + GAP_PENALTY

            # max(score_diag, score_up, score_left)
            best: int = score_diag if score_diag >= score_up else score_up
            best = best if best >= score_left else score_left
            dp[i, j] = best
            j += 1
        i += 1

    # Best score ending anywhere in last row or last column
    last_row_best: int = dp[num_rows - 1, 0]
    j = 1
    while j < num_cols:
        if dp[num_rows - 1, j] > last_row_best:
            last_row_best = dp[num_rows - 1, j]
        j += 1

    last_col_best: int = dp[0, num_cols - 1]
    i = 1
    while i < num_rows:
        if dp[i, num_cols - 1] > last_col_best:
            last_col_best = dp[i, num_cols - 1]
        i += 1

    return last_row_best if last_row_best >= last_col_best else last_col_best


def affine_global_align(seq1: str, seq2: str) -> int:
    """
    Compute the global alignment score with an affine gap penalty.

    Uses three DP matrices (float) to support -inf initialization:
      M (match/mismatch), Ix (gap in seq2), Iy (gap in seq1).
    Gap penalties: open + (k-1)*extend for a k-length gap.
    Returns only the final optimal alignment score (no traceback).
    """
    m: int = len(seq1)
    n: int = len(seq2)

    neg_inf: float = -float("inf")

    # DP matrices (float to allow use of -inf)
    M: ndarray[float, 2] = np.empty((m + 1, n + 1), dtype=float)
    Ix: ndarray[float, 2] = np.empty((m + 1, n + 1), dtype=float)
    Iy: ndarray[float, 2] = np.empty((m + 1, n + 1), dtype=float)

    # Fill with -inf
    M.fill(neg_inf)
    Ix.fill(neg_inf)
    Iy.fill(neg_inf)

    # Initialization
    M[0, 0] = 0.0

    # Leading gaps in seq2 (vertical moves)
    i: int = 1
    while i <= m:
        # GAP_OPEN_PENALTY + (i - 1) * GAP_EXTEND_PENALTY
        Ix[i, 0] = float(GAP_OPEN_PENALTY) + float(i - 1) * float(GAP_EXTEND_PENALTY)
        # M[i,0] and Iy[i,0] remain -inf
        i += 1

    # Leading gaps in seq1 (horizontal moves)
    j: int = 1
    while j <= n:
        Iy[0, j] = float(GAP_OPEN_PENALTY) + float(j - 1) * float(GAP_EXTEND_PENALTY)
        # M[0,j] and Ix[0,j] remain -inf
        j += 1

    # Fill
    i = 1
    while i <= m:
        a = seq1[i - 1]
        j = 1
        while j <= n:
            b = seq2[j - 1]
            subst: float = float(MATCH_SCORE if a == b else MISMATCH_PENALTY)

            # Match/mismatch state
            best_prev: float = M[i - 1, j - 1]
            if Ix[i - 1, j - 1] > best_prev:
                best_prev = Ix[i - 1, j - 1]
            if Iy[i - 1, j - 1] > best_prev:
                best_prev = Iy[i - 1, j - 1]
            M[i, j] = best_prev + subst

            # Gap in seq2 (vertical move)
            open_gap_vert: float = M[i - 1, j] + float(GAP_OPEN_PENALTY)
            extend_gap_vert: float = Ix[i - 1, j] + float(GAP_EXTEND_PENALTY)
            Ix[i, j] = (
                open_gap_vert if open_gap_vert >= extend_gap_vert else extend_gap_vert
            )

            # Gap in seq1 (horizontal move)
            open_gap_horiz: float = M[i, j - 1] + float(GAP_OPEN_PENALTY)
            extend_gap_horiz: float = Iy[i, j - 1] + float(GAP_EXTEND_PENALTY)
            Iy[i, j] = (
                open_gap_horiz
                if open_gap_horiz >= extend_gap_horiz
                else extend_gap_horiz
            )

            j += 1
        i += 1

    best_final: float = M[m, n]
    if Ix[m, n] > best_final:
        best_final = Ix[m, n]
    if Iy[m, n] > best_final:
        best_final = Iy[m, n]

    return int(best_final)


__all__ = [
    "global_align",
    "local_align",
    "semi_global_align",
    "affine_global_align",
]
