from collections import deque
import sys, copy

def flatten(L):
    for item in L:
        try:
            yield from flatten(item)
        except TypeError:
            yield item

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
    for i in range(mol_num):
        cutoffs = dict()
        molecules.append(list())
        print(f"Enter the types of atoms in the molecule {i+1}. When done enter 0")
        done = False
        while not done:
            str_input = input("Type: ")
            try:
                answ = int(str_input)
            except:
                answ = list(map(int, str_input.split()))
            if answ == 0:
                done = True
                print(f"Molecule {i+1}: {molecules[i]}")
            else:
                molecules[i].append(answ)

        molecules.append(True if input("Do you wish to highlight chained species?\nPlease answer Yes or No: ").lower() in ["yes", "y"] else False)

        print("Provide pair cutoff distances")

        flat_molecule = list(flatten(molecules[i]))
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
    chain = to_leave["molecules"][-1]
    molecules = to_leave["molecules"][:-1]
    for i, molecule in enumerate(molecules):
        copy_atoms = atom_data["atoms"].copy()
        for line in atom_data["atoms"]:
            if line[atom_data["type"]] in molecule:
                molecule_found = find_molecule(copy.deepcopy(molecule), line, atom_data, to_leave["cutoffs"][i], copy_atoms, chain) 
                ids_toleave += molecule_found["ids"]
                if molecule_found["ids"]:
                    print(f"Found molecule: {molecule_found['ids']} !!!")
                copy_atoms = molecule_found["copy_atoms"].copy()

    return [*set(ids_toleave)]
"""
def find_neighbor(atom_coords, source_type, target_type, atom_data, copy_atoms, exclude, cutoffs):
        closest = find_closest(atom_coords, atom_data, target_type, copy_atoms, exclude)
#        print(f"Closest atom to type {source_type} is: ", closest)
        if closest:
            closest_coords = (closest[atom_data["x"]], closest[atom_data["y"]], closest[atom_data["z"]])
            distance_between = dist(atom_coords, closest_coords, atom_data["box"])
#            print("Distance between them: ", distance_between)
            if distance_between > cutoffs[source_type][target_type]:
                return False
            else:
                exclude.append(closest[atom_data["id"]])
                return {"id": closest[atom_data["id"]], "coords":closest_coords}
        else:
            return False
"""
def find_neighbor(atom_coords, source_type, target_type, atom_data, copy_atoms, exclude, cutoffs):
    done = False
    result = list()
    while not done:
        closest = find_closest(atom_coords, atom_data, target_type, copy_atoms, exclude)
        if closest:
            closest_coords = (closest[atom_data["x"]], closest[atom_data["y"]], closest[atom_data["z"]])
            distance_between = dist(atom_coords, closest_coords, atom_data["box"])
            if distance_between <= cutoffs[source_type][target_type]:
                result.append({"id": closest[atom_data["id"]], "coords":closest_coords})
        else:
            return []

    
def find_molecule(molecule, atom, atom_data, cutoffs, copy_atoms, chain, exclude = "def"):
    print_condition = True if atom[atom_data["id"]] in [2331.0, 432.0, 1553.0, 4061.0, 751.0, 2992.0] else False
    print("Investigating a molecule: ", molecule) if print_condition else None
    print("Looking at neighbors of atom: ", atom)  if print_condition else None
    initial_flat = list(flatten(molecule))
    ids = [atom[atom_data["id"]]]
    starting_position = molecule.index(atom[atom_data["type"]])
    to_compare_index = starting_position
    atom_coords = (atom[atom_data["x"]], atom[atom_data["y"]], atom[atom_data["z"]])
    if exclude == "def":
        exclude = [atom[atom_data["id"]]]
    #Going backwards
    print("Going backwards") if print_condition else None
    for i in range(starting_position - 1, -1, -1):
        print("Searching for neighbors of atom at index: ", to_compare_index)  if print_condition else None
        if type(molecule[i]) == type(list()):
            print("The neighbor from the left is a molecule: ", molecule[i]) if print_condition else None
            neighbor = find_neighbor(atom_coords, molecule[to_compare_index], molecule[i][0], atom_data, copy_atoms, exclude, cutoffs)
            if neighbor:
                print("Found the attaching atom: ", neighbor["id"]) if print_condition else None
                atom_coords = neighbor["coords"]
                found_molecule = find_molecule(molecule[i], next((x for x in copy_atoms if x[atom_data["id"]] == neighbor["id"]), None), atom_data, cutoffs, copy_atoms, chain, exclude)
                if found_molecule["ids"]:
                    print("Found the intersticial molecule with ids: ", found_molecule["ids"]) if print_condition else None
                    ids += found_molecule["ids"]
                    copy_atoms = found_molecule["copy_atoms"]
                    molecule[i] = molecule[i][0]
                else:
                    print("Molecule not found ;(") if print_condition else None
                    break
                to_compare_index -= 1
                continue
            else:
                print("Attaching atom not found ;(") if print_condition else None
                break
            

        print("(after moltest)Looking for a neighbor for atom of type: ", molecule[to_compare_index]) if print_condition else None
        neighbor = find_neighbor(atom_coords, molecule[to_compare_index], molecule[i], atom_data, copy_atoms, exclude, cutoffs)
        if neighbor:
            print("Found the neighbor: ", neighbor["id"])  if print_condition else None
            ids.append(neighbor["id"])
            atom_coords = neighbor["coords"]
            to_compare_index -= 1
        else:
            print("Cound not find the neighbor") if print_condition else None
            break
        
    #Going forward
    print("Going forward") if print_condition else None
    to_compare_index = starting_position
    atom_coords = (atom[atom_data["x"]], atom[atom_data["y"]], atom[atom_data["z"]])
    for i in range(starting_position + 1, len(molecule)):
        print("Searching for neighbors of atom at index: ", to_compare_index) if print_condition else None 
        if type(molecule[i]) == type(list()):
            print("The neighbor from the right is a molecule: ", molecule[i]) if print_condition else None
            neighbor = find_neighbor(atom_coords, molecule[to_compare_index], molecule[i][0], atom_data, copy_atoms, exclude, cutoffs)
            if neighbor:
                print("Found the attaching atom: ", neighbor["id"]) if print_condition else None
                atom_coords = neighbor["coords"]
                found_molecule = find_molecule(molecule[i], next((x for x in copy_atoms if x[atom_data["id"]] == neighbor["id"]), None), atom_data, cutoffs, copy_atoms, chain, exclude)
                if found_molecule["ids"]:
                    print("Found the intersticial molecule with ids: ", found_molecule["ids"]) if print_condition else None
                    ids += found_molecule["ids"]
                    copy_atoms = found_molecule["copy_atoms"]
                    molecule[i] = molecule[i][0]
                else:
                    print("Molecule not found ;(") if print_condition else None
                    break
                to_compare_index += 1
                continue
            else:
                print("Attaching atom not found ;(") if print_condition else None
                break
            

        print("(after moltest)Looking for a neighbor for atom of type: ", molecule[to_compare_index]) if print_condition else None
        neighbor = find_neighbor(atom_coords, molecule[to_compare_index], molecule[i], atom_data, copy_atoms, exclude, cutoffs)
        if neighbor:
            print("Found the neighbor: ", neighbor["id"])  if print_condition else None
            ids.append(neighbor["id"])
            atom_coords = neighbor["coords"]
            to_compare_index += 1
        else:
            print("Cound not find the neighbor") if print_condition else None
            break


    result = {"ids":[], "copy_atoms":copy_atoms}
    if len(initial_flat) == len(ids):
        result["ids"] = ids
        if not chain:
            result["copy_atoms"] = list(filter(lambda x: False if x[atom_data["id"]] in ids else False, copy_atoms))
    return result


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
