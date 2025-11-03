# Week 5 & 6 Implementation Plan (Notebook-first, CI-friendly)

## Scope and key decisions

- Final notebook path: `week5_and_6/week5.ipynb`.
- IGV screenshots: manual (no headless automation). Include clear instructions + narrative prompts.
- Focus: hg38/GRCh38 `chr10` only. All downloads happen inside the notebook (self-contained).
- Tools: minimap2, samtools, bcftools, HapCUT2, whatshap (for HapCUT2→VCF conversion), IGV (manual use on desktop).

## Expected inputs and outputs

- Inputs (fetched in-notebook): `chr10.fa(.gz)`, Illumina interleaved FASTQ, PacBio FASTQ.
- Outputs:
  - Alignments: `illumina.sorted.bam(.bai)`, `pacbio.sorted.bam(.bai)`
  - Variants: `illumina.vcf.gz(.tbi)`, `pacbio.vcf.gz(.tbi)`
  - Phased: `illumina_phased.vcf.gz(.tbi)`, `pacbio_phased.vcf.gz(.tbi)`
  - Comparison: `vcf_compare/` with shared/unique sets and summary table
  - Figures: manual IGV PNGs (2–3 per gene)

## Notebook structure (section headers)

1) Environment + deps check (install via mamba/apt/pip if missing)

2) Parameters (URLs, threads, regions) + reproducibility seed

3) Reference genome (download chr10, index with `samtools faidx` and `minimap2 -d`)

4) Data download (Illumina interleaved, PacBio; decompress as needed)

5) (If interleaved) deinterleave Illumina → `illumina_R1.fq`/`illumina_R2.fq` (bbmap `reformat.sh` or fastp fallback)

6) Alignment (minimap2): Illumina with `-x sr`, PacBio with `-x map-pb` (or `-x map-hifi` if HiFi)

7) BAM postprocess: sort, index, flagstat; region subset (optional speed-up)

8) Variant calling (bcftools mpileup+call per sample; bgzip+tabix)

9) Phasing (HapCUT2: extractHAIRS + HAPCUT2; keep block output)

10) Convert HapCUT2 blocks → phased VCF (whatshap `hapcut2vcf`; bgzip+tabix)

11) Compare phased VCFs (bcftools isec; shared/unique counts; per-gene stats)

12) Manual IGV review instructions + placeholders for image uploads and written discussion

13) Star-allele interpretation (PharmVar lookup + phased haplotypes; write justification)

14) Runtime summary, pitfalls, and time estimate

15) Appendix: reproducibility notes, citations, and commands used

## Key commands (to embed in notebook cells)

- Reference (download + index):
  ```bash
  curl -L https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr10.fa.gz -o chr10.fa.gz
  gunzip -f chr10.fa.gz
  samtools faidx chr10.fa
  minimap2 -d chr10.mmi chr10.fa
  ```

- Illumina deinterleave (interleaved → R1/R2):
  ```bash
  mamba install -y -c bioconda bbmap || conda install -y -c bioconda bbmap
  reformat.sh in=illumina.fq out1=illumina_R1.fq out2=illumina_R2.fq overwrite=t
  ```

(Fallback) fastp:

  ```bash
  mamba install -y -c bioconda fastp || conda install -y -c bioconda fastp
  fastp --in1 illumina.fq --interleaved_in -o illumina_R1.fq -I illumina_R2.fq -w ${THREADS}
  ```

- Alignment:
  ```bash
  minimap2 -t ${THREADS} -ax sr chr10.mmi illumina_R1.fq illumina_R2.fq | samtools sort -@ ${THREADS} -o illumina.sorted.bam
  samtools index illumina.sorted.bam
  # If PacBio HiFi: -x map-hifi; otherwise -x map-pb
  minimap2 -t ${THREADS} -ax map-pb chr10.mmi pacbio.fq | samtools sort -@ ${THREADS} -o pacbio.sorted.bam
  samtools index pacbio.sorted.bam
  ```

- Variant calling (per sample):
  ```bash
  bcftools mpileup -Ou -f chr10.fa sample.sorted.bam \
    | bcftools call -mv -Oz -o sample.vcf.gz
  tabix -p vcf sample.vcf.gz
  ```

- Phasing (HapCUT2):
  ```bash
  extractHAIRS --bam sample.sorted.bam --VCF sample.vcf.gz --out sample.fragments
  HAPCUT2 --fragments sample.fragments --VCF sample.vcf.gz --output sample.hapcut
  ```

- Convert to phased VCF (whatshap):
  ```bash
  whatshap hapcut2vcf -o sample_phased.vcf sample.vcf.gz sample.hapcut
  bgzip -f sample_phased.vcf && tabix -p vcf sample_phased.vcf.gz
  ```

- Compare VCFs:
  ```bash
  mkdir -p vcf_compare
  bcftools isec -p vcf_compare illumina_phased.vcf.gz pacbio_phased.vcf.gz
  bcftools stats illumina_phased.vcf.gz pacbio_phased.vcf.gz > vcf_compare/stats.txt
  ```

## Regions for the three genes (hg38; verify once in notebook)

Create `genes.bed` in-notebook for speed/scoping. Example coordinates (verify via UCSC in code):

- CYP2C19: chr10:94,760,653-94,853,205
- CYP2C9:  chr10:96,696,685-96,748,843
- CYP2C8:  chr10:96,796,649-96,829,254

Then subset when desired, e.g. `bcftools view -R genes.bed` for per-gene tables.

## IGV (manual) – what to do and what to write

- Install/run IGV locally; load `chr10.fa`, `illumina.sorted.bam`, `pacbio.sorted.bam`, and both phased VCFs.
- For each gene: pick 2–3 discordant variants from `vcf_compare` outputs.
- Navigate to each locus; adjust coverage/zoom; take PNG snapshots.
- Save images under `igv/GENE_pos.png` (not tracked; paste into notebook cell outputs).
- In notebook, for each screenshot write:
  - Locus and gene
  - Which calls are discordant and why
  - Read-level evidence per technology (mapping quality, strand, homopolymers, indels, read length)
  - Your conclusion: artifact vs true variant (and rationale)

## CI notes (GitHub Actions)

- Use the provided template; add steps to:
  - Setup mamba with `conda-forge` + `bioconda`
  - Install: `minimap2 samtools bcftools hapcut2 whatshap bbmap fastp` + `python -m pip install -r week5_and_6/requirements.txt`
  - Execute: `jupyter nbconvert --to notebook --execute week5_and_6/week5.ipynb`
  - Cache conda pkgs; set a 60–90 min timeout

## Time estimates

- Setup + downloads: 1–2 h
- Alignment: 1–2 h
- Variant calling: 0.5–1 h
- Phasing + conversion: 0.5–1 h
- Comparison + per-gene summaries: 0.5–1 h
- Manual IGV + writeup: 1–2 h
- Star-alleles + writeup: 0.5–1 h

## Risks and mitigations

- Interleaved FASTQ: deinterleave with `bbmap reformat.sh` or `fastp`.
- PacBio type unknown: try `-x map-hifi` first for HiFi; fallback to `-x map-pb` if needed.
- CI timeouts: restrict to `genes.bed` early; reduce threads; enable caching.
- HapCUT2→VCF: use `whatshap hapcut2vcf` to avoid custom scripts.

### To-dos

- [ ] Create week5_and_6/week5.ipynb with all sections and parameters
- [ ] Add conda installs for minimap2/samtools/bcftools/hapcut2/whatshap/bbmap/fastp
- [ ] Download chr10.fa(.gz), index with samtools and minimap2
- [ ] Fetch Illumina (interleaved) and PacBio FASTQ; decompress
- [ ] Split illumina.fq into illumina_R1.fq and illumina_R2.fq
- [ ] Align both technologies with minimap2; sort/index BAMs
- [ ] Call variants per sample (bcftools); bgzip+tabix
- [ ] Run HapCUT2 (extractHAIRS + HAPCUT2) per sample
- [ ] Convert HapCUT2 blocks to phased VCF with whatshap
- [ ] Compute shared/unique variants; per-gene summaries
- [ ] Capture manual IGV screenshots; paste and discuss in notebook
- [ ] Assign star-alleles per gene using PharmVar and phased VCFs
- [ ] Update CI to install tools and execute notebook
- [ ] Polish notebook, add timing, ensure reproducibility
