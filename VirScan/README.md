# VirScan: Longitudinal Analysis of EBV Proteins in Young Adults

## Overview
This dataset and analysis focus on young adults who experienced **infectious mononucleosis (mono)** caused by primary **Epstein-Barr Virus (EBV)** infection. The study tracks these individuals longitudinally to determine:
- **Which EBV proteins are expressed** at different stages of infection.
- **Potential vaccine candidates** based on protein expression profiles.

The data were generated using **VirScan**, a high-throughput method for analyzing antibody responses to viral infections.

---

## Objectives
1. Identify **EBV proteins** expressed during various stages of infection.
2. Assess the **immune response** to these proteins over time.
3. Evaluate the feasibility of targeting these proteins as **vaccine candidates**.

---

## Dataset Description
The dataset contains the results of antibody profiling against EBV proteins. It includes:
- **Participant Information**: Young adults with primary EBV infection.
- **Time Points**: Longitudinal samples collected at various stages post-infection.
- **Protein Data**: Antibody responses to a comprehensive library of EBV proteins.

### Directory Structure
- `data/`: Raw and processed VirScan datasets.
- `scripts/`: Analysis scripts to process and visualize the data.
- `results/`: Figures, processed data, and key findings.

---

## Key Findings
1. Certain EBV proteins are consistently expressed during specific stages of infection.
2. Immune responses to some proteins show promise for vaccine development.

---

## How to Use the Data
### Accessing the Data
1. Navigate to the `data/` directory.
2. Download the datasets for analysis.

### Running the Analysis
1. Clone this repository:
   ```bash
   git clone https://github.com/ganta86/computational-virology.git
   cd computational-virology/VirScan
   ```
2. Run the analysis scripts in the `scripts/` directory:
   ```bash
   python analyze_virscan.py
   ```

### Results
- Processed results can be found in the `results/` directory, including visualizations and statistical analyses.

---

## Citation
If you use this dataset or analysis in your work, please cite:
> Ganta, et al. (2025). Longitudinal Analysis of EBV Proteins in Young Adults Using VirScan. Computational Virology Repository: [https://github.com/ganta86/computational-virology](https://github.com/ganta86/computational-virology)

---

## Contact
For questions or issues regarding this dataset or analysis, please contact:
- **Name**: [Your Full Name]
- **Email**: [Your Email Address]
- **GitHub**: [https://github.com/ganta86](https://github.com/ganta86)
