import streamlit as st
from rdkit import Chem
import pandas as pd
from rdkit import DataStructs

def _smiles_to_num_dict():
    cox_molecules_data = pd.read_csv("data/Cox-molecules/Cox-molecules-overview.csv")
    novel_molecules_data = pd.read_csv("data/novel-molecules/novel-molecules-data.csv")
    smiles_to_num_dict = {}
    for smiles, molecule_number in zip(cox_molecules_data["smiles"].values, cox_molecules_data["molecule_number"].values):
        smiles_to_num_dict[smiles] = molecule_number

    for smiles, molecule_number in zip(novel_molecules_data["smiles"].values, novel_molecules_data["molecule_number"].values):
        smiles_to_num_dict[smiles] = molecule_number
        
    return smiles_to_num_dict

def _contains_boronic_acid(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        boronic_acid_pattern = Chem.MolFromSmarts('B(O)O')
        matches = mol.GetSubstructMatches(boronic_acid_pattern)
        return len(matches) > 0
    else:
        return False
    
def _find_most_similar(smiles, smiles_to_num_dict):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        query_fp = Chem.RDKFingerprint(mol)
        most_similar_smiles = None
        highest_similarity = 0.0
        for ref_smiles in smiles_to_num_dict.keys():
            ref_mol = Chem.MolFromSmiles(ref_smiles)
            if ref_mol is not None:
                ref_fp = Chem.RDKFingerprint(ref_mol)
                similarity = DataStructs.TanimotoSimilarity(query_fp, ref_fp)
                if similarity > highest_similarity:
                    highest_similarity = similarity
                    most_similar_smiles = ref_smiles
        return most_similar_smiles
    else:
        return None
    
def main():
    st.title("Protodeboronation Prediction")
    st.write("Protodeboronation can be a big problem in cross-coupling reactions that use boronic acids (see example below). This website helps you by predicting the rate of protodeboronation for your molecule of interest! If you found this work helpful, please consider citing our [our paper](https://doi.org/10.1021/acs.jpca.2c08250)")
    # Displaying an image
    image = 'figures/pdb_example.png'
    st.image(image, caption='Suzuki reaction with protodeboronation.', use_column_width=True)
    
    # Instructions expander
    with st.expander("How do I find the SMILES of my molecule?"):
        st.write('- Step 1: Google the name of your molecule (e.g. "(3,5-Dimethyl-1,2-oxazol-4-yl)boronic Acid")')
        st.write("- Step 2: Click a suitable link (e.g. PubChem, Wikipedia, Sigma Aldrich)")
        st.write('- Step 3: Click "Ctrl + F" ("Cmd + F" on mac) and search for "SMILES"')
        st.write("To learn more about molecular representation and SMILES, check out [this paper]( https://doi.org/10.1002/wcms.1603).")
    
    # Text input field
    smiles = st.text_input("Enter boronic acid SMILES:", value="Cc1noc(C)c1B(O)O")
    
    #canon smiles
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        smiles = Chem.MolToSmiles(mol)
        # Check if the smiles is a boronic acid
        contains_boronic_acid = _contains_boronic_acid(smiles)
    else:
        smiles = None

        
    # Button
    if st.button("Go"):
        if smiles is None:
            st.write("Please enter a valid SMILES!")
            return
        if not contains_boronic_acid:
            st.write("Please enter a boronic acid!")
            return
        else:
            smiles_to_num_dict = _smiles_to_num_dict()
            
            # if smiles in smiles_to_num_dict:
            if smiles in smiles_to_num_dict.keys():
                st.write("Boronic acid found in our dataset!")
                num = smiles_to_num_dict[smiles]
                if num <=100:
                    st.image(f"figures/Cox/{num}.png", caption='Protodeboronation prediction', use_column_width=True)
                else:
                    st.image(f"figures/novel/{num}.png", caption='Protodeboronation prediction', use_column_width=True)
            
            # If we don't have data for that BA, we need to find the most similar molecule
            else:
                st.write("Boronic acid not found in our dataset. Finding the most similar molecule...")
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    most_similar = _find_most_similar(smiles, smiles_to_num_dict)
                    most_similar_num = smiles_to_num_dict[most_similar]
                    
                    if most_similar is not None:
                        col1, col2 = st.columns(2)
                        
                        with col1:
                            st.write("Queried Structure")
                            mol = Chem.MolFromSmiles(smiles)
                            if mol:
                                image = Chem.Draw.MolToImage(mol)
                                st.image(image, caption='Queried Structure', use_column_width=True)
                            else:
                                st.write("Invalid SMILES input!")

                        with col2:
                            st.write("Most Similar Structure")
                            mol = Chem.MolFromSmiles(most_similar)
                            if mol:
                                image = Chem.Draw.MolToImage(mol)
                                st.image(image, caption='Most Similar Structure', use_column_width=True)
                            else:
                                st.write("No similar structure found!")
                                
                    st.write("Protodeboronation prediction for the most similar molecule in our dataset:")
                    if most_similar_num <=100:
                        st.image(f"figures/Cox/{most_similar_num}.png", caption='Protodeboronation prediction', use_column_width=True)
                    else:
                        st.image(f"figures/novel/{most_similar_num}.png", caption='Protodeboronation prediction', use_column_width=True)
                        
        st.write("NB: Protodeboronation rate is on a logarithmic scale. See the corresponding half-life below:")
        st.markdown("- log(k) = -10 &rarr; 200 years")
        st.markdown("- log(k) = -8 &rarr; 2 years")
        st.markdown("- log(k) = -6 &rarr; 10 days")
        st.markdown("- log(k) = -4 &rarr; 2 hours")
        st.markdown("- log(k) = -2 &rarr; 1 minute")
        st.markdown("- log(k) = 0 &rarr; 1 second")
        st.markdown("- log(k) = 2 &rarr; 0.01 seconds")
                    
                    
                    

if __name__ == "__main__":
    main()