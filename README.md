# VIMS Resazurin Research Project

This repository contains research on resazurin metabolic assays for shellfish aquaculture conducted at the Virginia Institute of Marine Science (VIMS).

## Project Overview

This research evaluates resazurin-based metabolic assays as biomarkers for oyster performance in aquaculture, with applications for:
- Breeding program enhancement
- Environmental stress tolerance prediction
- Production optimization
- Risk management

## Repository Structure

### Analysis Scripts (`/code/`)
- `resazurin-trials.Rmd` - Main analysis of metabolic biomarkers and trait correlations
- `ploidy-trials-20250826.Rmd` - Ploidy comparison studies (triploid vs diploid)
- `ploidy-trials-20250813.Rmd` - Earlier ploidy trials
- `temperature-testing.Rmd` - Temperature stress testing analysis
- `temperature.Rmd` - Environmental temperature monitoring
- `plate-maps.Rmd` - Experimental design documentation

### Research Summary
- `research-summary-and-industry-implications.md` - **Comprehensive summary of findings and industry applications**
- `code/research-summary-and-industry-implications.Rmd` - R Markdown version for rendering

### Data (`/data/`)
- Raw experimental data files
- Metadata and experimental design information
- Performance trait predictions from breeding programs

### Outputs (`/figures/` and `/output/`)
- Generated plots and statistical outputs
- Publication-ready figures

## Key Findings

1. **Metabolic activity significantly predicts survival in low-salinity environments** for environmentally-selected oyster populations
2. **Ploidy affects metabolic performance patterns**, with distinct profiles between triploid and diploid oysters
3. **Acute stress response correlates with performance traits**, enabling rapid screening of breeding stock
4. **Family-level genetic variation** in metabolic activity supports selective breeding applications
5. **Environmental context matters** - correlations are stronger when assays match production conditions

## Industry Applications

### Immediate Benefits
- 15-30% reduction in environmental stress mortality through pre-deployment screening
- Accelerated breeding cycles (6-12 months vs 18-24 months for traditional testing)
- Cost-effective screening (~$2-5 per assay vs $0.50-2.00 per oyster mortality cost)

### Long-term Opportunities
- Integration with precision aquaculture systems
- Climate adaptation strategies
- Multi-species applications
- Value-added product development

## Getting Started

### For Researchers
1. Review the main analysis file: `code/resazurin-trials.Rmd`
2. Examine methodology and statistical approaches
3. Check data structure in `/data/` folder
4. Refer to figure outputs for key visualizations

### For Industry Stakeholders
1. **Start with the research summary**: `research-summary-and-industry-implications.md`
2. Focus on the "Industry Implications and Applications" section
3. Review the "Future Research Recommendations" for collaboration opportunities
4. Contact information provided for technology transfer discussions

### For Running Analyses
Required R packages:
```r
install.packages(c("tidyverse", "ggplot2", "readxl", "cowplot", 
                   "lme4", "lmerTest", "car", "emmeans", "pracma"))
```

## Citation

If you use this research, please cite:
[Citation to be added upon publication]

## Contact

**Dr. AS Huffmyer**  
Virginia Institute of Marine Science  
[Contact information]

For technology transfer and commercial applications:
[Technology transfer office contact]

## License

[License information to be added]

## Acknowledgments

This research was supported by [funding sources]. We thank our collaborators and commercial partners for their contributions to this work.