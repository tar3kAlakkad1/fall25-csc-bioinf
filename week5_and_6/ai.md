# AI Assistance Documentation - Week 5 & 6

## Tools Used

- Claude Code (claude.ai/code) - Primary development assistant
- Claude.ai (Sonnet 3.5 with extended thinking) - Research and troubleshooting
- GitHub Copilot - Code completion and suggestions

## Planning Phase

### Initial Setup

- **Task:** Created comprehensive plan (Plan.md)
- **AI Contribution:** 100% - Claude Code generated the complete implementation plan
- **Details:** Structured 15-step plan covering all deliverable requirements

## AI Assistance by Section

### Section 1: Reference Genome Download

- **Tasks:**

1. Navigate to the [Genome Browser](https://genome.ucsc.edu/) and search for CYP2C8, CYP2C9, CYP2C19.
2. Saw `chr10` in the top left of the screen (see screenshot below) and assumed this is the desired chromosome, as instructor mentions in the deliverable description that we should "focus on the chromosome that contains these genes!"
3. Asked Claude.ai (model: Sonnet 3.5 with extended thinking enabled) to determine the best way to download the target chromosome.
4. Downloaded the necessary file by running:

```
curl https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr10.fa.gz --output chr10.fa.gz
gunzip chr10.fa.gz
```

- **AI Help:** For determining how to download the necessary files and creating index files
- **Code Generated:** The `curl` command, samtools faidx, and minimap2 indexing commands
- **Percentage AI:** 25% of this task was completed with AI assistance

### Section 2: Sequencing Data Download

- **Tasks:** Download Illumina and PacBio FASTQ files, handle compression formats
- **AI Help:** Suggested using environment variables for URLs, handling both .gz and .bz2 formats
- **Code Generated:** Complete bash script for flexible data download with fallbacks
- **Percentage AI:** 40% - AI provided the robust download logic with error handling

### Section 3: Alignment with minimap2

- **Tasks:** Align short and long reads with appropriate parameters
- **AI Help:** Recommended correct presets (-ax sr for Illumina, -ax map-pb for PacBio CLR)
- **Code Generated:** Complete alignment pipeline with sorting and indexing
- **Percentage AI:** 70% - AI provided the full pipeline structure

### Section 4: Variant Calling

- **Tasks:** Call variants using bcftools restricted to gene regions
- **AI Help:** Suggested using genes.bed file for efficiency, proper bcftools mpileup + call pipeline
- **Code Generated:** Full variant calling commands with compression and indexing
- **Percentage AI:** 40% - AI generated most of the variant calling pipeline but had to be modified / fixed due to errors

### Section 5: Variant Phasing

- **Tasks:** Phase variants with HapCUT2 and convert to VCF format
- **AI Help:** Provided complete HapCUT2 workflow including extractHAIRS parameters for different technologies
- **Code Generated:** Complete phasing pipeline with technology-specific parameters
- **Percentage AI:** 90% - Complex HapCUT2 pipeline largely AI-generated

### Section 6: VCF Comparison

- **Tasks:** Compare VCFs between technologies, calculate concordance
- **AI Help:** Suggested bcftools isec for comparison, provided statistics calculations
- **Code Generated:** VCF comparison and per-gene analysis scripts
- **Percentage AI:** 50% - AI provided comparison methodology and scripts

### Section 7: IGV Visualization

- **Tasks:** Generate IGV batch scripts and automate screenshot capture
- **AI Help:** Created complete IGV batch automation system with loci selection
- **Code Generated:** Python code for batch script generation, bash commands for IGV execution
- **Percentage AI:** 40% - IGV automation is complex and was primarily AI-generated, but code was not working. So, I had to run some terminal commands and install IGV with brew via `brew install igv`, then generated manual screenshots with `igv -b week5_and_6/igv/manual_batch.txt`, and finally updated the code to just read and load those files.

### Section 8: Star-Allele Determination

- **Tasks:** Determine CYP2C star-alleles based on detected variants
- **AI Help:** Provided PharmVar database knowledge, star-allele definitions, and clinical interpretations
- **Code Generated:** Complete Python analysis with star-allele database and matching logic
- **Percentage AI:** 40% - Comprehensive star-allele analysis system was AI-designed but had to be modified again as well.

## Debugging and Troubleshooting

### Tool Installation Issues

- **Problem:** HapCUT2 not available via conda on macOS ARM64
- **AI Assistance:** Provided local compilation instructions and dynamic library path fixes
- **Solution:** Built HapCUT2 locally and set DYLD_FALLBACK_LIBRARY_PATH for htslib

### Pipeline Errors

- **Problem:** IGV batch mode timing out when generating screenshots
- **AI Assistance:** Suggested alternative approaches including manual batch files
- **Solution:** Created manual_batch.txt and ran IGV locally with `igv -b manual_batch.txt`

### Parameter Selection

- **Task:** Choosing correct minimap2 presets for different sequencing technologies
- **AI Help:** Explained difference between map-pb (CLR) and map-hifi (CCS) for PacBio
- **Final Choice:** Used map-pb for CLR reads with configurable PACBIO_PRESET variable

## Parameter and Method Choices

### Alignment Parameters

- minimap2 for Illumina: `-ax sr` (short-read preset)
- minimap2 for PacBio: `-ax map-pb` (PacBio CLR preset)
- **AI Assistance:** Provided correct preset parameters and explained their purposes

### Variant Calling Settings

- Tool chosen: bcftools mpileup + call
- Quality threshold: Default (QUAL > 20)
- Depth threshold: Default (adequate for this dataset)
- **AI Assistance:** Suggested restricting to genes.bed for efficiency

### Phasing Tool

- Tool chosen: HapCUT2
- **AI Assistance:** Provided complete workflow including technology-specific extractHAIRS parameters

## Overall Assessment

**Total AI Contribution:** Approximately 55-60%

**Areas of significant AI help:**

1. HapCUT2 phasing pipeline - complex tool with limited documentation
2. IGV batch automation - sophisticated scripting for GUI tool automation
3. Star-allele determination - required PharmVar database knowledge and clinical interpretation
4. Clinical recommendations - drug-specific dosing guidelines based on genotypes

**Areas of independent work:**

1. Initial gene location research using UCSC Genome Browser
2. Manual IGV screenshot generation when automation failed
3. Troubleshooting tool installation on macOS ARM64
4. Data file organization and notebook structure

## Time Estimate

**Total time:** ~ 24 hours

- Planning and research: 2 hours
- Tool installation and setup: ~ 45 minutes
- Pipeline implementation: 13 hours
- Debugging and documentation: 8 hours

## Notes

- AI was essential for understanding complex bioinformatics tools (HapCUT2, bcftools)
- Most bash scripts were AI-generated but manually tested and refined
- Star-allele clinical interpretations heavily relied on AI's PharmVar knowledge
- IGV automation was the most AI-dependent section due to complexity
