import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from sklearn.cluster import KMeans

#1. Load the Data file in SDF file format
file_name = input("Enter the SDF file name(with extension): ")

#Read the SDF file
supplier = Chem.SDMolSupplier(file_name)

#Extract molecular structures
data = [mol for mol in supplier if mol is not None]

#Extract descriptors by creating a list 
descriptors = []
for mol in data:
    molecular_weight = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    descriptors.append([molecular_weight, logp])

#Converting descriptors to arrays using NumPy
X = np.array(descriptors)

#Choose the Number of clusters to be formed in the data
num_clusters = 3

#To perform K-means Clustering by setting n_init to suppress the warning
kmeans = KMeans(n_clusters=num_clusters, n_init=10)
kmeans.fit(X)

#Cluster labels
cluster_labels = kmeans.labels_

#Add cluster labels to the SDF data
for mol, cluster_label in zip(data, cluster_labels):
    mol.SetProp("Cluster", str(int(cluster_label)))

#to save cluster assigned data to an SDF file
output_name = input("Enter the name of the output SDF file: ")
writer = Chem.SDWriter(output_name)
for mol in data:
    writer.write(mol)
writer.close()
