# Metabolomics Analysis of Three Distinct Challenges: Fasting, Physical Activity and Oral Lipid Tolerance Test (OLTT)

With this code a subset of the HuMet dataset was analyzed, focusing on **plasma** samples across three platforms:

- Metabolon HD4 (nt-ms)
- Biocrates p150 (t-ms)
- In-house biochemistry (chem.)

### The pipeline contains the following steps:
1. Since the original dataset lacked challenge information, it was assigned based on the time column:

- **Fasting**: Time points 1–10
- **Physical Activity**: Time points 33–39
- **Oral Lipid Tolerance Test (OLTT)**: Time points 40–50
  
2. Metabolites with >30% missing values were removed
3. NA values were imputed using *randomForest* R package
4. Z-Score Normalisation
5. Hypothesis Testing with Anova-Like Test
6. Analysis of metabolites unique for each challenge (o-acetylhomoserine, fumarate, taurocholate)
7. Analysis of inter-individual variation
8. Analysis of log2 fold changes of different metabolites
  
