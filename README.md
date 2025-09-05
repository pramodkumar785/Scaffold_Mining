# Scaffold Mining using K-means Clustering

This repository contains a small Python-based pipeline for **scaffold mining and clustering of chemical compounds**.  
It was developed as part of my Master's research work in **computational drug discovery**.

The scripts use **RDKit**, **scikit-learn**, and **matplotlib** to cluster compounds based on molecular descriptors, apply **Lipinski’s rule of five** as a drug-likeness filter, and retrieve/visualize compounds of interest.

---

## 🚀 Features
- **K-means clustering** of compounds using molecular descriptors (MolWt, LogP, HBD, HBA).  
- **Lipinski’s filtering** to select drug-like molecules.  
- **Cluster visualization**: scatter plots with compound numbers.  
- **Compound retrieval**: extract molecules by index with SMILES, MolBlock, and 2D structure images.  
- Simple and modular — each script does one step in the workflow.

---

## 📂 Workflow

1. **1_Kmeans_Clustering.py**  
   - Input: `.sdf` file of compounds.  
   - Extracts descriptors (MolWt, LogP).  
   - Runs K-means clustering.  
   - Saves a new `.sdf` with **Cluster labels** in the properties.

2. **2_Kmeans_Plot_Lipinski.py**  
   - Input: clustered `.sdf` file.  
   - Filters molecules using **Lipinski’s Rule of Five**.  
   - Performs clustering on filtered compounds.  
   - Generates a **scatter plot (MolWt vs LogP)** with compound numbers.

3. **3_Extraction&Retrieval.py**  
   - Input: clustered `.sdf` file.  
   - Retrieve a compound by its number.  
   - Outputs:  
     - SMILES string  
     - Molecular weight  
     - MolBlock representation  
     - 2D structure image (`.png`)

---

## ⚙️ Dependencies
Install the following Python libraries before running the scripts:

- [RDKit](https://www.rdkit.org/)  
- NumPy  
- Pandas  
- scikit-learn  
- matplotlib  
- Pillow (PIL)  

---

## 🖥️ Example Usage

```bash
# Step 1: Assign clusters to molecules
python 1_Kmeans_Clustering.py molecules.sdf

# Step 2: Filter with Lipinski’s rule and plot clusters
python 2_Kmeans_Plot_Lipinski.py clustered_molecules.sdf

# Step 3: Extract and save a specific compound
python 3_Extraction&Retrieval.py clustered_molecules.sdf


📌 Notes

The pipeline is designed for educational and research purposes.

Input files must be in .sdf format.

Default number of clusters is 3, but this can be changed in the scripts.

Code is lightweight and can be adapted for larger datasets or integrated into other drug discovery workflows.

✨ Author

Developed by Hridhya Nair
