import pandas as pd
from os import walk
from rdkit import Chem
from rdkit.Chem import AllChem


class NotConvergedError(Exception):
    pass

def num_to_atom(atom_nums_list):
    # Maps an atomic number to the corresponding atom string
    atom_symb = []
    
    for i in range(len(atom_nums_list)):
        if atom_nums_list[i] == 1:
            atom_symb += ['H']
        if atom_nums_list[i] == 5:
            atom_symb += ['B']
        if atom_nums_list[i] == 6:
            atom_symb += ['C']
        if atom_nums_list[i] == 7:  
            atom_symb += ['N']
        if atom_nums_list[i] == 8:
            atom_symb += ['O']
        if atom_nums_list[i] == 9:
            atom_symb += ['F']
        if atom_nums_list[i] == 15:
            atom_symb += ['P']
        if atom_nums_list[i] == 16:
            atom_symb += ['S']
        if atom_nums_list[i] == 17:
            atom_symb += ['Cl']
        if atom_nums_list[i] == 35:
            atom_symb += ['Br']
        if atom_nums_list[i] == 53:
            atom_symb += ['I']
    
    if len(atom_symb) != len(atom_nums_list):
        print('Unknown atom used')
    len(atom_symb) == len(atom_nums_list)
            
    return atom_symb

def extract_coordinates(log_file_path):
    """
    input: filepath to a .log file
    returns: list of strings containing the final coordinates
    """
    f = open(f"{log_file_path}", "r")
    output_file_lines = f.readlines()


    # Print out which molecules didn't converge so I can fix them
    # There's an error if it doesn't end in 'Normal termination of Gaussian 16'
    if output_file_lines[-1][:34] != ' Normal termination of Gaussian 16':  
        print('Not optimised correctly: ', log_file_path)      

    
    try:
        # Find latest coordinates
        for i in range(len(output_file_lines)):

            #find charge:
            if output_file_lines[i][:10] == ' Charge = ':
                charge = int(output_file_lines[i][9:12])

            if output_file_lines[i][25:46] == 'Standard orientation:':
                #read off coordinates in the lines below here
                coordinates = output_file_lines[i:i+100]

        # cut away the header
        coordinates = coordinates[5:]

        # cut away footer
        for i in range(len(coordinates)):
            if coordinates[i][3] == '-':
                coordinates = coordinates[:i]
                break

        l_l_l = [x.split('    ') for x in coordinates]

        df = pd.DataFrame(l_l_l)

        atom_nums = list(df[3])
        atom_list = [int(x) for x in atom_nums]

        atom_symb_list = num_to_atom(atom_list)

        # extract the actual coordinates:
        coordinates = [x[36:-2] for x in coordinates]

        #pad the atom_symb_list to coordinates:

        for i in range(len(coordinates)):
            coordinates[i] = atom_symb_list[i]+ coordinates[i]+'\n'

        return coordinates, charge
    except IndexError:
        print('Could not find coordinates in: ', log_file_path)




def write_com_file(log_file_path, file_name, output_folder_path, header):
    try:
        coordinates, charge = extract_coordinates(log_file_path)
        
        ##n M06L/6-311++G** SCRF=(Solvent=Water) Opt=(TS,CalcFC,noeigentest)
    
        #header
        _header = []
        _header += ['%NProcShared=2\n']
        _header += header



        _header += f"""
 Title

{charge} 1
"""

        # create a new .com file with the modified input_file_lines
        if file_name[-6:-4] == '_2':
            f_name = file_name[:-6]
        else:
            f_name = file_name[:-4]
        with open(f"{output_folder_path}/{f_name}.com", 'w') as new_input_file:
            new_input_file.writelines(_header)
            new_input_file.writelines(coordinates)
            new_input_file.writelines('\n')
            new_input_file.close()
    except TypeError:
        error = 1


def prepare_coordinates_for_SI(mol_num, mechanism, role, log_file_path):
    """
    mol_num: e.g. 61 (molecule number)
    mechanism: e.g. 'k2'
    role: e.g. 'TS' or 'int'
    """
    list_of_lines = []
    list_of_lines += [r'\begin{lstlisting}']
    list_of_lines += [f'Molecule {mol_num}: {mechanism} - {role}']
    
    #The basis set used is simplified to B3LYP in the case of k2cat.
    #The solvent used is simplified to in vacuo (as opposed to water) in the case of k2cat
    
    if mechanism != 'k2cat':
        if role == 'TS' or role == 'Transition State':
            list_of_lines += ['#n M06L/6-311++G** SCRF=(Solvent=Water) Opt=(TS,CalcFC,noeigentest)']
        elif role != 'TS':
            list_of_lines += ['#n M06L/6-311++G** SCRF=(Solvent=Water) Opt']
    elif mechanism == 'k2cat':
        if role == 'TS' or role == 'Transition State':
            list_of_lines += ['#n B3LYP/6-31G(d) Opt=(TS,CalcFC,noeigentest)']
        elif role != 'TS':
            list_of_lines += ['#n B3LYP/6-31G(d) Opt']
    
    # insert energy:
    energy = read_energy(log_file_path)
    list_of_lines += [energy]

    
    #insert charge, multiplicity, and coordinates
    coordinates, charge = extract_coordinates(log_file_path)
    
    list_of_lines += [f'Charge, multiplicity: {charge}, 1']
    
    list_of_lines += ['Geometry:']
    
    list_of_lines += coordinates
    
    list_of_lines += [r'\end{lstlisting}']
    
    #replace all instances of \n
    lines = [x.replace('\n','') for x in list_of_lines]
    
    return lines



def generate_coordinates_file_for_SI(folder_path):
    """
    # Write a txt file which I can copy-paste into SI containing the 
    #  coordinates for all molecules for a mechanism

    gen_SI_coor: generate the coordinates for SI
    folder_path: the folder containing all the .log files
    molecule_class: 'Cox' or 'novel'
    """
    folders = folder_path.split('/')
    mechanism = folders[-1]
    filenames = next(walk(folder_path), (None, None, []))[2]  # [] if no file
    filenames.sort()

    if 'Cox' in folder_path:
        output_path = f'data/Cox-molecules/DFT/coordinates_for_SI/{mechanism}_coor.txt'
    elif 'novel' in folder_path:
        output_path = f'data/novel-molecules/DFT/coordinates_for_SI/{mechanism}_coor.txt'
    
    #write empty txt file:
    with open(output_path, 'w') as f:
    
        for filename in filenames:
            if filename[-4:] == '.log':
                if filename[2] == '_': #if this is true, then it's a 2 digit Cox molecule
                    mol_num = int(filename[:2])
                elif filename[3] == '_':
                    mol_num = int(filename[:3])

                if 'TS' in filename:
                    role = 'Transition State'
                elif 'int' in filename:
                    role = 'Intermediate'
                elif 'Ar-' in filename:
                    role = 'Intermediate, BA detached'
                elif 'reactant' in filename:
                    role = 'Reactant'
                elif '00' in filename:
                    role = 'Misc. small molecule'
                    
                lines = prepare_coordinates_for_SI(mol_num, mechanism, role, folder_path+'/'+filename)
                for line in lines:
                    f.write(line)
                    f.write('\n')
        f.close()
    print('Coordinates created successfully')
            
            
# Generate .com file given a smiles string
def generate_com_file_from_smiles(smiles, writePath, charge, TS, TS_key_atoms = None ):
    #TS is either true or false, depending on whether it is a transition state smiles
    m2 = Chem.MolFromSmiles(smiles)
    m3 = Chem.AddHs(m2)
    AllChem.EmbedMolecule(m3,randomSeed=0xf00d)   # optional random seed for reproducibility)
    molfile = Chem.MolToMolBlock(m3)
    df = pd.DataFrame([x.split('   ') for x in molfile.split('\n')])
    df = df.iloc[4: , 1:4]
    df = df.dropna()
    df = df.reset_index()
    
    #strip leading spaces
    df[3] = df[3].str.strip()
    
    #if there's a charge in the molecule, we need to remove the last line (CHG in the molfile)
#     if '+]' in smiles or '-]' in smiles:
#         #
#         print("There's a charge in ",smiles, ' so I deleted the following:') 
#         print(df.iloc[-1])
#         df.drop(df.tail(1).index,inplace=True) # drop last n rows

    x = df[1]
    y = df[2]
    
    z_and_atom_df = pd.DataFrame([x.split(' ') for x in df[3]])
    z = z_and_atom_df[0]
    atom = z_and_atom_df[1]
    
    new_df = pd.concat([atom, x,y,z],axis=1)
    # df is now done, it has the correct coordinates
    
    # now need to write the .com file
    # I will use B3LYP
    
    with open(writePath, 'w') as f:
        # header of the file
        if TS == True:
            header = f"""%NProcShared=2
#n M06L/6-311++G** SCRF=(Solvent=Water) Opt=(TS,CalcFC,noeigentest)

 Title

{charge} 1
"""
#         elif TS == False:
#             header = f"""%NProcShared=2
# #n B3LYP/6-31G(d) SCRF=(Solvent=Water) Opt

#  Title

# {charge} 1
# """

        elif TS == False:
            header = f"""%NProcShared=2
#n M06L/6-311++G** SCRF=(Solvent=Water) Opt

 Title

{charge} 1
"""
            
        f.write(header)
        
        # coordinates
        df_as_string = new_df.to_string(header=False, index=False)
        f.write(df_as_string)
        
        if TS_key_atoms != None:
            f.write(TS_key_atoms)
        f.write("""
        
""")
        f.close()
    

def generate_k2Ar_minus_com_files(df, mechanism_name, output_folder):
    """
    df: contains info on the relevant molecules
    mechanism_name = 'k2Ar'
    """
    list_of_molecules = []
    for i in range(len(df['mechanisms_active'])):
        if mechanism_name == df['mechanisms_active'][i]:

            smiles = df['smiles'][i]
            molecule_number = df['molecule_number'][i]
            
            # Change to smiles string due to mechanism
            if mechanism_name == 'k2Ar':
                #remove either OB(O) or B(O)O (which may have () around it)
                Ar = smiles.replace('cOB(O)','[c-]')
                if Ar == smiles:
                    Ar = smiles.replace('OB(O)c','[c-]')
                    if Ar == smiles:
                        Ar = smiles.replace('c(B(O)O)','[c-]')
                        if Ar == smiles:
                            Ar = smiles.replace('OB(O)C','[C-]')
                    


                list_of_molecules +=[Ar]
                write_path = output_folder + f'/{molecule_number}_k2Ar_Ar-.com'
                charge = -1
                TS = False
                generate_com_file_from_smiles(Ar, write_path, charge, TS)
    return list_of_molecules
                
                

#k2Ar TS
def generate_k2Ar_TS_com_files(df, mechanism_name, output_folder):
    """
    df: contains info on the relevant molecules
    mechanism_name = 'k2Ar'
    """
    list_of_molecules = []
    for i in range(len(df['mechanisms_active'])):
        if mechanism_name == df['mechanisms_active'][i]:

            smiles = df['smiles'][i]
            molecule_number = df['molecule_number'][i]
            
            # Change to smiles string due to mechanism
            if mechanism_name == 'k2Ar':
                #remove either OB(O) or B(O)O (which may have () around it)
                Ar = smiles.replace('B(O)','[B-](O)')
                Ar = Ar.replace('(O)', '(O)(O)')
            

                list_of_molecules +=[Ar]
                write_path = output_folder + f'/{molecule_number}_k2Ar_TS.com'
                charge = -1
                TS = True
                generate_com_file_from_smiles(Ar, write_path, charge, TS)

    return list_of_molecules
                
def read_energy(log_file_path):
    f = open(log_file_path, "r")
    output_file_lines = f.readlines()

    # Only read the energy if the file ends with 'Normal termination of Gaussian 16'
    if output_file_lines[-1][:34] == ' Normal termination of Gaussian 16':

        #check each line to see if it starts with 'SCF Done:'
        for i in range(len(output_file_lines)):
            if output_file_lines[i][:10] == ' SCF Done:':
                #read off energy from this line
                #energy = output_file_lines[i][24:38]
                energy = output_file_lines[i][12:38]
                # each time it finds a new 'SCF Done', the new energy measurement
                # will be overwritten by the new one
    return energy
        

def read_energy_in_whole_folder(folder_path):
    """
    Cox molecules are only 2 digits
    Novel molecules are 3 digits
    Cox_molecules = True or False
    """
    filenames = next(walk(folder_path), (None, None, []))[2]  # [] if no file
    
    #check if it's a cox molecule or a novel molecule
    if 'Cox' in folder_path:
        Cox_molecules = True
    else:
        Cox_molecules = False
    
    molecule_dict = {}
    
    number_of_errors = 0
    number_of_files = 0
    
    for filename in filenames:
        
        # We are only interested in the .log files
        if filename[-4:] == '.log':
            number_of_files +=1
            
            if Cox_molecules:
                molecule = filename[0:2]
                mechanism = filename[3:-4]
                
            elif not Cox_molecules:
        
                molecule = filename[0:3]
                mechanism = filename[4:-4]

            #open the file to find the energy
            f = open(f"{folder_path}/{filename}", "r")
            output_file_lines = f.readlines()

            # Only read the energy if the file ends with 'Normal termination of Gaussian 16'
            if output_file_lines[-1][:34] == ' Normal termination of Gaussian 16':

                #check each line to see if it starts with 'SCF Done:'
                for i in range(len(output_file_lines)):
                    if output_file_lines[i][:10] == ' SCF Done:':
                        #read off energy from this line
                        energy = output_file_lines[i][24:38]
                        # each time it finds a new 'SCF Done', the new energy measurement
                        # will be overwritten by the new one

            else:
                #probably timeout error
                energy = 'Error_occured'
                number_of_errors += 1
                
                # If gaussian didn't have 'normal termination', it was probably due to 
                # a timeout error.             
                

            #now add molecule, mechanism, and energy to list/dict

            # if it's not already a dict:
            try:
                molecule_dict[molecule][mechanism] = energy
            except KeyError:
                molecule_dict[molecule] = {}
                molecule_dict[molecule][mechanism] = energy                
        
    # Now create a df from molecule_dict
    molecule_df = pd.DataFrame.from_dict(molecule_dict, orient='index')
    molecule_df.index.name = 'molecule_number'
    
    # Now create a csv from this file
    molecule_df.to_csv(f'{folder_path}/energy.csv')
    print(number_of_errors, ' errors out of ', number_of_files, ' calculations.')