# Epidemiology Simulation of Dengue

A group project simulating the dengue epidemic using a mathematical SIRSI model integrated with mosquito population genetics. MATLAB is used for all numerical simulations and visualizations.

---

## 📖 Overview

This project simulates the spread of dengue, combining epidemiological models and population genetics to analyze how insecticide resistance in mosquitoes affects disease dynamics. The SIRSI compartmental model is used for both human and mosquito populations.

---

## 🚀 Features

- **Three evolutionary scenarios** of insecticide resistance genes in vectors:
  - Equal genotype mortality (Hardy-Weinberg equilibrium)
  - Homozygous recessive resistance (aa) dominates
  - Resistance suppressed (aa less fit than AA/Aa)
- **Coupled ODE system** for compartmental populations, solved by Euler’s method
- **Visual plots** for:
  - Human populations (susceptible, infected, recovered)
  - Vector populations (genotype-wise trends)
  - Genotype and allele frequency evolution

---

## ⚙️ Getting Started

### Requirements
- MATLAB (2019a or later recommended)

### How to Run
1. Clone/download the repository and extract files
2. Open MATLAB
3. Run the provided script

Edit parameters in the script as desired.

---

## 🧪 Methodology

- SIRSI mathematical model: S, I, R compartments for humans; S, I compartments for three vector genotypes (AA, Aa, aa)
- Simulation assumptions: random mating, only biological selection, horizontal disease transmission, no migration
- Ordinary differential equations solved with Euler’s method

---

## 💡 Key Insights

- Presence of insecticide resistance genes increases the number of infected humans and vectors, impacting epidemic severity
- Genotype and allele frequencies shift due to natural selection
---

## 📝 Improvements

- Model vector life stages and various dengue serotypes
- Add migration between populations
- Include additional vector traits

---

## 👥 Contributors

- Abiruth S 
- Rithish A B 
- S.N. Anand Mahadev 
- Sriranjana C 
- Thrishala S N
  
**Faculty Guide:** Dr. Ambika P S, Amrita Vishwa Vidyapeetham

---

## 📚 References

- [PMC5501426](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5501426/)
- [PMC6136472](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6136472/)



