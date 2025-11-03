# Week 5 & 6 Deliverable - Setup Complete

## Files Created

1. **Plan.md** - Comprehensive 15-step implementation plan
   - Detailed instructions for each section
   - Tool installation guides
   - Expected outputs
   - Timeline estimates (~25-35 hours total)

2. **week5.ipynb** - Jupyter notebook skeleton
   - Pre-structured with all required sections
   - Ready for implementation with Claude/Cursor
   - Contains TODO markers for each step

3. **ai.md** - AI assistance tracking template
   - Ready to document AI contributions
   - Structured by section
   - Includes percentage tracking

4. **.gitignore** (root) - Updated to exclude large data files
   - Prevents committing BAM/VCF/FASTQ files
   - Includes all bioinformatics file formats

## Quick Start

### 1. Review the Plan
```bash
cat Plan.md
```
Read through the complete implementation plan to understand the workflow.

### 2. Open Jupyter Notebook
```bash
cd week5_and_6
jupyter notebook week5.ipynb
```

### 3. Use Claude/Cursor for Implementation
- Work through each section sequentially
- Use Plan.md as a reference guide
- Update ai.md as you use AI assistance

## What You Need to Do

### Immediate Actions

1. **Get Data URLs** (CRITICAL)
   - Check Piazza/assignment for Illumina and PacBio download links
   - Update notebook cells with actual URLs

2. **Set Up CI** (Before starting work)
   - Update `.github/workflows/actions.yml`
   - Add tool installations (minimap2, samtools, bcftools, HapCUT2)
   - Set timeout to 60 minutes
   - Add jupyter execution step

### Implementation Order

Follow the plan sections 1-15:

**Week 1 (Days 1-3):**
- Setup environment
- Download reference genome (chr10)
- Download sequencing data
- Configure CI

**Week 1-2 (Days 4-6):**
- Implement alignment with minimap2
- Call variants with bcftools
- Generate statistics

**Week 2 (Days 7-9):**
- Phase variants with HapCUT2
- Compare VCFs
- Identify discordant variants
- Create IGV screenshots

**Week 2 (Days 9-10):**
- Research PharmVar
- Determine star-alleles
- Write conclusions
- Test and debug

## Key Requirements

### Self-Contained Notebook
- All downloads must happen in notebook cells
- No manual steps allowed
- Running `jupyter execute week5.ipynb` should produce all results

### Outputs Required
1. Reference: `chr10.fa`
2. Alignments: `illumina.bam`, `pacbio.bam` (+ `.bai` indices)
3. Variants: `illumina.vcf.gz`, `pacbio.vcf.gz`
4. Phased: `illumina_phased.vcf.gz`, `pacbio_phased.vcf.gz`
5. Analysis: VCF comparison statistics
6. Figures: IGV screenshots (2-3 per gene)
7. Results: Star-allele calls with justification

### Grading (6 points total)
- 1 pt: Reference genome download
- 1 pt: Alignment
- 1 pt: Variant calling
- 1 pt: Phasing and conversion
- 1 pt: Variant analysis
- 1 pt: Star-allele calls
- +0.5 pt bonus: Automated IGV screenshots

## Important Notes

### Gene Locations (all on chr10)
- **CYP2C19:** chr10:94,760,653-94,853,205
- **CYP2C9:** chr10:96,696,685-96,748,843
- **CYP2C8:** chr10:96,796,649-96,829,254

### Tools Needed
- minimap2 (alignment)
- samtools (BAM manipulation)
- bcftools (variant calling)
- HapCUT2 (phasing)
- IGV or igv-reports (visualization)
- Python packages: pysam, pandas, biopython, matplotlib

### CI Configuration
Add to `.github/workflows/actions.yml`:
```yaml
- name: Execute notebook
  run: jupyter nbconvert --to notebook --execute week5_and_6/week5.ipynb
```

Set appropriate timeout (pipeline takes ~30-60 minutes).

## Resources

### Documentation
- **Plan.md** - Complete implementation guide
- **ai.md** - Track AI assistance here

### External Resources
- UCSC Genome Browser: https://genome.ucsc.edu/
- PharmVar: https://www.pharmvar.org/
- minimap2: https://github.com/lh3/minimap2
- HapCUT2: https://github.com/vibansal/HapCUT2
- IGV: https://igv.org/

## Troubleshooting

### If tools aren't installed
```bash
# Using conda
conda install -c bioconda minimap2 samtools bcftools

# Or via apt (Linux)
sudo apt-get install minimap2 samtools bcftools
```

### If HapCUT2 isn't available
```bash
git clone https://github.com/vibansal/HapCUT2
cd HapCUT2
make
# Add to PATH
```

### If notebook fails to execute
- Test each section individually
- Check tool availability
- Verify data file paths
- Look for missing dependencies

## Next Steps

1. Review Plan.md thoroughly
2. Obtain data download URLs from Piazza
3. Configure CI workflow
4. Start implementation with Section 1
5. Update ai.md as you progress
6. Test incrementally
7. Submit before deadline (Sun, Nov 2, 23:59)

## Questions?

- Post on Piazza
- Attend office hours
- Check Plan.md for detailed guidance

---

**Ready to start? Open week5.ipynb and begin with Section 1!**

Good luck! 🧬
