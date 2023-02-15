from collections import deque
import numpy as np
import sys

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
        split_row = row.split()
        if split_row:
            new_row = []
            for i, instance in enumerate(split_row):
                if i == results_dict["type"]:
                    new_row.append(int(instance))
                else:
                    new_row.append(float(instance))
            results_dict["atoms"].append(new_row)
#    results_dict["atoms"] = list(filter(None, [list(map(lambda x: float if split_row.index(x) != results_dict["type"] else int, split_row)) for row in que]))
    return results_dict

def dist(atom1, atom2, box):
    x = abs(atom1[0] - atom2[0])
    x = (box[0] - x) if (x > box[0]/2) else x
    
    y = abs(atom1[1] - atom2[1])
    y = (box[1] - y) if (y > box[1]/2) else y
    
    z = abs(atom1[2] - atom2[2])
    z = (box[2] - z) if (z > box[2]/2) else z

    return (x**2 + y**2 + z**2)**(0.5)
    
def read_input():
    mol_num = int(input("Enter the number of molecules you wish to highlight: "))
    molecules = list()
    cutoffs = dict()
    for i in range(mol_num):
        molecules.append(list())
        print(f"Enter the types of atoms in the molecule {i+1}. When done enter 0")
        done = False
        while not done:
            answ = int(input("Type: "))
            if answ == 0:
                done = True
                print(f"Molecule {i+1}: {molecules[i]}")
            else:
                molecules[i].append(answ)

        molecules.append(True if input("Do you wish to highlight chained species?\nPlease answer Yes or No: ").lower() in ["yes", "y"] else False)

        print("Provide pair cutoff distances")

        molecule_types_num = len(molecules[i]) - 1 
        for k in range(molecule_types_num):
            if molecules[i][k] not in cutoffs:
                cutoffs[molecules[i][k]] = dict()
            for j in range(k + 1, molecule_types_num):
                if molecules[i][k] != molecules[i][j]:
                    if molecules[i][j] not in cutoffs[molecules[i][k]]:
                        cutoffs[molecules[i][k]][molecules[i][j]] = float(input(f"type pair {molecules[i][k]}-{molecules[i][j]}: "))
                    if molecules[i][j] not in cutoffs:
                        cutoffs[molecules[i][j]] = dict()
                    cutoffs[molecules[i][j]][molecules[i][k]] = cutoffs[molecules[i][k]][molecules[i][j]]

    return {"molecules": molecules, "cutoffs": cutoffs}

def get_molecules_ids(to_leave, atom_data):
    ids_toleave = []
    chain = to_leave["molecules"][-1]
    molecules = to_leave["molecules"][:-1]
    for molecule in molecules:
        copy_atoms = atom_data["atoms"].copy()
        for line in atom_data["atoms"]:
            if line[atom_data["type"]] in molecule:
                ids_toleave += find_molecule(molecule, line, atom_data, to_leave["cutoffs"], copy_atoms, chain)

    return [*set(ids_toleave)]
    
def find_molecule(molecule, atom, atom_data, cutoffs, copy_atoms, chain):
    ids = [atom[atom_data["id"]]]
    starting_position = molecule.index(atom[atom_data["type"]])
    to_compare_index = starting_position
    atom_coords = (atom[atom_data["x"]], atom[atom_data["y"]], atom[atom_data["z"]])
    exclude = [atom[atom_data["id"]]]
    #Going backwards
    for i in range(starting_position - 1, -1, -1):
        done = False
        while not done:
            closest = find_closest(atom_coords, atom_data, molecule[i], copy_atoms, exclude)
            if closest:
                closest_coords = (closest[atom_data["x"]], closest[atom_data["y"]], closest[atom_data["z"]])
                distance_between = dist(atom_coords, closest_coords, atom_data["box"])
                if distance_between > cutoffs[molecule[to_compare_index]][molecule[i]]:
                    done = True
                else:
                    ids.append(closest[atom_data["id"]])
                    exclude.append(closest[atom_data["id"]])
            else:
                done = True
        to_compare_index -= 1
        
    #Going forward
    to_compare_index = starting_position
    for i in range(starting_position + 1, len(molecule)):
        done = False
        while not done:
            closest = find_closest(atom_coords, atom_data, molecule[i], copy_atoms, exclude)
            if closest:
                closest_coords = (closest[atom_data["x"]], closest[atom_data["y"]], closest[atom_data["z"]])
                distance_between = dist(atom_coords, closest_coords, atom_data["box"])
                if distance_between > cutoffs[molecule[to_compare_index]][molecule[i]]:
                    done = True
                else:
                    ids.append(closest[atom_data["id"]])
                    exclude.append(closest[atom_data["id"]])
            else:
                done = True
        to_compare_index += 1


    if len(molecule) == len(ids):
        if not chain:
            for line in copy_atoms:
                if line[atom_data["id"]] in ids:
                    copy_atoms.remove(line)
        return ids
    else:
        return []


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
