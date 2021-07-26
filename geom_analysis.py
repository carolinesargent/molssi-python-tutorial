import numpy as np
import os
import argparse

def calculate_distance(atom1_coord, atom2_coord):
    """Calculate the distance between two three-dimensional points."""
    
    x_distance = atom1_coord[0] - atom2_coord[0]
    y_distance = atom1_coord[1] - atom2_coord[1]
    z_distance = atom1_coord[2] - atom2_coord[2]
    bond_length_12 = np.sqrt(x_distance ** 2 + y_distance ** 2 + z_distance ** 2)
    return bond_length_12

def bond_check(atom_distance, minimum_length=0, maximum_length=1.5):
    """
    Check if a distance is a bond based on a minimum and maximum bond length.
    """
    if atom_distance < 0:
        raise ValueError(F'Invalid atom distance {atom_distance}. Distance can not be less than 0!')
    if atom_distance > minimum_length and atom_distance <= maximum_length:
        return True
    else:
        return False
    
def open_xyz(file_location):
    fpath, extension = os.path.splitext(file_location)

    if extension.lower() != '.xyz':
        raise ValueError("Incorrect file type! File must be type xyz")

    xyz_file = np.genfromtxt(fname = file_location, skip_header=2, dtype='unicode')
    atom_labels = xyz_file[:,0]
    coords = xyz_file[:,1:]
    coords = coords.astype(np.float64)
    return atom_labels, coords

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This script analyzes a user given xyz file and outputs the length of the bonds.")
    parser.add_argument("xyz_file", help="The filepath for the xyz file to analyze")
    parser.add_argument('-minimum_length', help='The minimum distance to consider atoms bonded.', type=float, default=0)
    parser.add_argument('-maximum_length', help='The maximum distance to consider atoms bonded.', type=float, default=1.5)
    args = parser.parse_args()

    atom_labels, coords = open_xyz(args.xyz_file)
    num_atoms = len(atom_labels)

    for num1 in range(0, num_atoms):
        for num2 in range(0, num_atoms):
            if num1 < num2:
                bond_length_12 = calculate_distance(coords[num1], coords[num2])
                if bond_check(bond_length_12, minimum_length=args.minimum_length, maximum_length=args.maximum_length) is True:
                    print(F'{atom_labels[num1]} to {atom_labels[num2]} : {bond_length_12:.3f}')
