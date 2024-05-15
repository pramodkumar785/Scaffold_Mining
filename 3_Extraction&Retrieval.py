from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
from PIL import Image

# Function to extract a molecule by compound number from an SDF file
def extract_molecule(sdf_file_path, compound_number):
    # Load the SDF file
    suppl = Chem.SDMolSupplier(sdf_file_path)

    # Extract the molecule with the specified compound number
    compound = next((mol for idx, mol in enumerate(suppl, start=1) if idx == compound_number), None)

    return compound

# Main function for testing
def main():
    # Get the path to the SDF file with cluster assignments from user input (including the extension)
    sdf_file_path = input("Enter the name of the SDF file with cluster assignments (including the extension): ")

    # Get the compound number to extract
    compound_number = int(input("Enter the Compound Number to extract: "))

    # Extract the molecule
    extracted_molecule = extract_molecule(sdf_file_path, compound_number)

    # Check if the molecule was found
    if extracted_molecule is not None:
        # Get additional details
        smiles_string = Chem.MolToSmiles(extracted_molecule)
        molecular_weight = Descriptors.ExactMolWt(extracted_molecule)

        # Print details
        print(f"Molecule with Compound Number {compound_number}:")
        print(f"SMILES String: {smiles_string}")
        print(f"Molecular Weight: {molecular_weight} Da")
        print("Structure (MDL MolBlock):")
        print(Chem.MolToMolBlock(extracted_molecule))

        # Visualize the structure using the AllChem module
        AllChem.Compute2DCoords(extracted_molecule)
        img = Chem.Draw.MolToImage(extracted_molecule, size=(300, 300))

        # Save the image using PIL (Python Imaging Library)
        img.save(f"compound_{compound_number}.png")
        print(f"Structure image saved as compound_{compound_number}.png")
    else:
        print(f"Compound with Compound Number {compound_number} not found.")

if __name__ == "__main__":
    main()
