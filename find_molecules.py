from collections import deque
from typing import Iterable
import sys, copy

#Flattens nested lists like so: [1, [1, 2, [3], 4], 5, 6] -> [1, 1, 2, 3, 4, 5, 6]
def flatten(items):
    for x in items:
        if isinstance(x, Iterable) and not isinstance(x, (str, bytes)):
            for sub_x in flatten(x):
                yield sub_x
        else:
            yield x

#Analizes the input file and produces a dictionary with a list of atoms present in the system and the indecies of their specific properties like type and position
def get_data(filename):
    file = open(filename, "r")
    contents = file.read()
    que = deque(contents.split("\n"))
    results_dict = dict()

    line = ""
    header = ""
    done = False
    while not done:
        line = que.popleft()
        header += (line + '\n')
        if "ITEM: BOX BOUNDS" in line:
            results_dict["box"] = list()
            for i in range(3):
                line = que.popleft()
                header += (line + '\n')
                bounds = line.split()
                results_dict["box"].append(float(bounds[1]) - float(bounds[0]))
        if "ITEM: NUMBER OF ATOMS" in line:
            line = que.popleft()
            results_dict["numatoms"] = line
            header += line + '\n'
        if "ITEM: ATOMS" in line:
            split_header = line.split()[2:]
            results_dict["id"] = split_header.index("id")
            results_dict["type"] = split_header.index("type")
            results_dict["x"] = split_header.index("x")
            results_dict["y"] = split_header.index("y")
            results_dict["z"] = split_header.index("z")
            done = True

    results_dict["header"] = header

    results_dict["atoms"] = []
    for row in que:
        if "Velocities" in row:
            break
        split_row = row.split()
        if split_row:
            new_row = []
            for i, instance in enumerate(split_row):
                if i == results_dict["type"]:
                    new_row.append(instance)
                else:
                    new_row.append(float(instance))
            results_dict["atoms"].append(new_row)
#    results_dict["atoms"] = list(filter(None, [list(map(lambda x: float if split_row.index(x) != results_dict["type"] else int, split_row)) for row in que]))
    return results_dict

def dist(atom1, atom2, box):
#    print(atom1, atom2)
    x = abs(atom1[0] - atom2[0])
    x = (box[0] - x) if (x > box[0]/2) else x
    
    y = abs(atom1[1] - atom2[1])
    y = (box[1] - y) if (y > box[1]/2) else y
    
    z = abs(atom1[2] - atom2[2])
    z = (box[2] - z) if (z > box[2]/2) else z
#    print(x, y, z)

    return (x**2 + y**2 + z**2)**(0.5)
    
def read_input():
    mol_num = int(input("Enter the number of molecules you wish to highlight: "))
    molecules = list()
    cutoffs_list = list()
    print("\nA space-separated sequence of types on one line will be treated as a separate atom chain within the main one. To specify the anchor atom put the ampersand in front of the type (1 &2 1).")
    for i in range(mol_num):
        cutoffs = dict()
        molecules.append(list())
        print(f"\nEnter the types of atoms in the molecule {i+1}. When done enter 0")
        done = False
        while not done:
            str_input = input("Type: ")
            if " " in str_input:
                answ = str_input.split()
            else:
                answ = str_input
            if answ == "0":
                done = True
                print(f"Molecule {i+1}: {molecules[i]}")
            else:
                molecules[i].append(answ)

        molecules[i].append(True if input("Do you wish to highlight chained species?\nPlease answer Yes or No: ").lower() in ["yes", "y"] else False)

        print("Provide pair cutoff distances")

        flat_molecule = list(flatten(molecules[i]))
        flat_molecule = [el.replace("&", "") for el in flat_molecule[:-1]]
        molecule_types_num = len(flat_molecule) - 1 
        for k in range(molecule_types_num):
            if flat_molecule[k] not in cutoffs:
                cutoffs[flat_molecule[k]] = dict()
            for j in range(k + 1, molecule_types_num):
                if flat_molecule[k] != flat_molecule[j]:
                    if flat_molecule[j] not in cutoffs[flat_molecule[k]]:
                        cutoffs[flat_molecule[k]][flat_molecule[j]] = float(input(f"type pair {flat_molecule[k]}-{flat_molecule[j]}: "))
                    if flat_molecule[j] not in cutoffs:
                        cutoffs[flat_molecule[j]] = dict()
                    cutoffs[flat_molecule[j]][flat_molecule[k]] = cutoffs[flat_molecule[k]][flat_molecule[j]]
        cutoffs_list.append(cutoffs)


    return {"molecules": molecules, "cutoffs": cutoffs_list}

def get_molecules_ids(to_leave, atom_data):
    ids_toleave = []
    molecules = to_leave["molecules"]
    for i, molecule in enumerate(molecules):
        mol_counter = 0
        chain = molecule[-1]
        molecule = molecule[:-1]
        copy_atoms = atom_data["atoms"].copy()
        for line in atom_data["atoms"]:
            if molecule and line[atom_data["type"]] in molecule:
                molecule_found = find_molecule(copy.deepcopy(molecule), line, atom_data, to_leave["cutoffs"][i], copy_atoms, molecule.index(line[atom_data["type"]])) 
                ids_toleave += molecule_found
                if molecule_found:
                    mol_counter += 1
                    print(f"Found molecule #{mol_counter} of species {i+1}: {molecule_found}")
                if not chain:
                    copy_atoms = list(filter(lambda x: False if x[atom_data["id"]] in molecule_found else True, copy_atoms))

    return [*set(ids_toleave)]

def find_neighbor(atom_coords, source_type, target_type, atom_data, copy_atoms, exclude, cutoffs):
    done = False
    result = list()
    while not done:
        closest = find_closest(atom_coords, atom_data, target_type, copy_atoms, exclude)
        if closest:
            closest_coords = (closest[atom_data["x"]], closest[atom_data["y"]], closest[atom_data["z"]])
            distance_between = dist(atom_coords, closest_coords, atom_data["box"])
            if distance_between <= cutoffs[source_type][target_type]:
                result.append(closest[atom_data["id"]])
                exclude.append(closest[atom_data["id"]])
            else:
                done = True
        else:
            done = True

    for i in range(len(result)):
        exclude.pop()
    return result

def extract_anchor(molecule):
    pos = "default"
    for i, el in enumerate(molecule):
        if "&" in el:
            pos = i
    if pos == "default":
        return 0, molecule
    new_molecule = copy.deepcopy(molecule)
    new_molecule[pos] = new_molecule[pos].replace("&", "")
    return pos, new_molecule

def find_molecule(molecule, atom, atom_data, cutoffs, copy_atoms, rel_position, exclude = []):
#    print("Molecule: ", molecule)
#    print("Atom: ", atom)
#    print("Relative position: ", rel_position)
    molecule_ids = [atom[atom_data["id"]]]
    #Setup the list of atoms that have already been found
    if not exclude:
        exclude = [atom[atom_data["id"]]]
    #Check the basecase where only one atom in the requested molecule
    initial_flat = list(flatten(molecule))
    if len(initial_flat) == 1:
        return [atom[atom_data["id"]]]
    if type(molecule[rel_position]) == type(list()):
        anchor_position, clean_molecule = extract_anchor(molecule[rel_position])
        inner_molecule = find_molecule(clean_molecule, atom, atom_data, cutoffs, copy_atoms, anchor_position, exclude)
        if not inner_molecule:
            #Inner molecule not found so the molecule dne
            return []
        #Inner molecule found
        found_molecule = find_molecule(molecule[:rel_position] + [clean_molecule[anchor_position]] + molecule[rel_position + 1:], atom, atom_data, cutoffs, copy_atoms, rel_position, exclude)
        if len(molecule) > 1 and not found_molecule:
            return []
        else:
            #Return the ids of the inner molecule and the outer one
            return list(set(inner_molecule + found_molecule))

    else:
        #First atom position is ocupied by an atom
        atom_coords = (atom[atom_data["x"]], atom[atom_data["y"]], atom[atom_data["z"]])
        #Going forward
        if len(molecule) - 1 > rel_position:
            if type(molecule[rel_position + 1]) == type(list()):
                #Next spot to the right is a molecule
                anchor_position, clean_molecule = extract_anchor(molecule[rel_position + 1])
                neighbors = find_neighbor(atom_coords, molecule[rel_position], clean_molecule[anchor_position], atom_data, copy_atoms, exclude, cutoffs)
            else:
                #Next position is an atom
                neighbors = find_neighbor(atom_coords, molecule[rel_position], molecule[rel_position + 1], atom_data, copy_atoms, exclude, cutoffs)
            #For each possible neighbor at the current step recursively look for the molecule as if the neighbor was the correct choice
            for neighbor in neighbors:
                exclude.append(neighbor)
                found_molecule = find_molecule(molecule[rel_position + 1:], next((x for x in copy_atoms if x[atom_data["id"]] == neighbor), None), atom_data, cutoffs, copy_atoms, 0,  exclude)
                if found_molecule:
                    #Molecule on the right exists
                    molecule_ids += found_molecule
                    break
                else:
                    #Right part of the molecule not found, so the whole molecule dne
                    return []
        if rel_position > 0:
            if type(molecule[rel_position - 1]) == type(list()):
                #Next spot to the left is a molecule
                anchor_position, clean_molecule = extract_anchor(molecule[rel_position - 1])
                neighbors = find_neighbor(atom_coords, molecule[rel_position], clean_molecule[anchor_position], atom_data, copy_atoms, exclude, cutoffs)
            else:
                #Next position is an atom
                neighbors = find_neighbor(atom_coords, molecule[rel_position], molecule[rel_position - 1], atom_data, copy_atoms, exclude, cutoffs)
            #For each possible neighbor at the current step recursively look for the molecule as if the neighbor was the correct choice
            for neighbor in neighbors:
                exclude.append(neighbor)
                found_molecule = find_molecule(molecule[:rel_position], next((x for x in copy_atoms if x[atom_data["id"]] == neighbor), None), atom_data, cutoffs, copy_atoms, rel_position - 1,  exclude)
                if found_molecule:
                    #Molecule on the right exists
                    molecule_ids += found_molecule
                    break
                else:
                    #Left part of the molecule not found, so the whole molecule dne
                    return []
#        print("IDS: ", molecule_ids, "Molecule: ", molecule)
        if len(molecule_ids) != len(list(flatten(molecule))):
            return []
        else:
            return molecule_ids

        

def find_closest(atom, atom_data, atom_type, copy_atoms, exclude):
    filtered_by_type = list(filter(lambda x: True if x[atom_data["type"]] == atom_type else False, copy_atoms))
    filtered_without_excluded = []
    for line in filtered_by_type:
        if line[atom_data["id"]] not in exclude:
            filtered_without_excluded.append(line)
    sorted_atoms = sorted(filtered_without_excluded, key = lambda x: dist(atom, (x[atom_data["x"]], x[atom_data["y"]], x[atom_data["z"]]), atom_data["box"]))
    if sorted_atoms:
        return sorted_atoms[0]
    else:
        return 0

def create_new_file(filename, atom_data, ids):
    split_filename = filename.split('.')
    file = open('.'.join(split_filename[:-1]) + ".molecules." + split_filename[-1], "w")
    file.write(atom_data["header"].replace(f"\n{atom_data['numatoms']}\n", f"\n{len(ids)}\n"))
    for line in atom_data["atoms"]:
        if line[atom_data["id"]] in ids:
            file.write(' '.join([str(el) for el in line]) + '\n')

    file.close()


def main():
    filename = sys.argv[1]
    print("Analyzing the inputted file")
    atom_data = get_data(filename)
    type_column = atom_data["type"]
    to_leave = read_input()
    ids_to_leave = get_molecules_ids(to_leave, atom_data)
    create_new_file(filename, atom_data, ids_to_leave)

if __name__ == "__main__":
    main()
