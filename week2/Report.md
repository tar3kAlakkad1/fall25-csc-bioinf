## Documenting steps (for bash script later)
Python:
1) copied over trviz necessary files
2) added `pytest` to `requirments.txt`
3) ran `pytest` and passed all tests

Codon:
1) Ported `motif_encoder.py` to Codon (`motif_encoder.codon`)
   - Replaced Python's `Dict` and `List` with Codon's `dict` and `list` syntax
   - Fixed string formatting to use f-strings instead of `.format()`
   - Handled `Counter` from `collections` which Codon supports
   - Removed `Optional` for class attributes where possible, initialized with empty containers
   - Handled type inference issues with mixed dictionaries in `get_score_matrix()`
   
2) Ported `test_motif_encoder.py` to Codon (`test_motif_encoder.codon`)
   - Converted pytest-based tests to plain function-based tests with assertions
   - Removed pytest fixtures and replaced with inline data structures
   - Removed `tmp_path` fixture and used `/tmp/` directory directly
   - All 17 tests pass successfully

3) Updated `utils.codon` to fix type inference in `get_score_matrix()`
   - Changed return type to `dict[str, dict[str, float]]` to satisfy Codon's type system
   - Stored gap penalties as nested dictionaries instead of direct float values
   - This maintains compatibility while satisfying static typing requirements

Remember: for codon tests, may need: export CODON_PYTHON=/Users/tarekalakkadp/.pyenv/versions/3.11.9/lib/libpython3.11.dylib
