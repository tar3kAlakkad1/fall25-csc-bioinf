# AI Usage Documentation - Week 6

## Tools Used
- Claude 3.5 Sonnet (via Cursor IDE) - Primary development assistant

## AI-Assisted Development Process

### 1. Initial Setup and Understanding
- Used AI to understand the deliverable requirements and create a structured approach
- AI helped identify the key components needed: simpleaf/alevin-fry, Scanpy, CellTypist

### 2. Data Handling Strategy
- AI assisted in solving the Box.com download authentication issue
- Helped determine that including data files in the repository was the most CI-friendly approach
- Generated code for handling tar.gz extraction and file management

### 3. Simpleaf/Alevin-fry Pipeline
- AI provided guidance on simpleaf command-line arguments
- Helped structure the index building and quantification steps
- Assisted with proper environment variable setup (ALEVIN_FRY_HOME)

### 4. Scanpy Analysis Pipeline
- AI helped implement robust preprocessing for small datasets
- Provided fallback strategies for HVG selection on toy data
- Assisted with proper PCA component selection for limited genes/cells

### 5. CellTypist Integration
- AI helped debug gene symbol mapping issues
- Provided multiple fallback strategies for model downloading
- Assisted with handling the case of no overlapping genes in toy dataset

### 6. CI/CD Configuration
- AI helped structure the GitHub Actions workflow
- Provided proper micromamba environment setup
- Assisted with dependency management across conda and pip

### 7. Notebook Structure
- AI helped create a self-contained, reproducible notebook
- Assisted with proper markdown documentation
- Helped implement proper error handling and fallback mechanisms

## Key AI Contributions

1. **Error Handling**: Comprehensive try-catch blocks with informative error messages
2. **Robustness**: Multiple fallback strategies for each processing step
3. **Documentation**: Clear section headers and explanations throughout
4. **CI Compatibility**: Proper handling of file paths and environment variables
5. **Time Tracking**: Implementation of runtime measurement and effort estimation

## Time Estimate
- Initial implementation: ~4 hours
- Debugging and testing: ~2 hours  
- CI setup and testing: ~1.5 hours
- Documentation and cleanup: ~0.5 hours
- **Total**: ~8 hours

## Lessons Learned
- Small toy datasets require special handling in bioinformatics pipelines
- CellTypist models expect specific gene naming conventions
- Self-contained notebooks need careful data management strategies
- CI environments require explicit dependency specifications
