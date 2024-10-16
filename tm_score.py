import numpy as np
import argparse

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--template_path', type=str, default='A572547.pdb', help='The template structure pdb file path.')
    parser.add_argument('--native_path', type=str, default='B572547.pdb', help='The native structure pdb file path.')

    return parser.parse_args()

args = parse_args()

def calculate_distance_matrix(native_coords, template_coords):
    """
    calculate the distance between corresponding residue pairs
    parameters:
    native_coords: (N, 3)
    template_coords: (N, 3)
    return:
    distances: (N,)
    """
    return np.sqrt(np.sum((native_coords - template_coords) ** 2, axis=1))

def calculate_rmsd(native_coords, template_coords):
    """
    calculate the RMSD between two sets of atomic coordinates
    parameters:
    native_coords: (N, 3)
    template_coords: (N, 3)
    return:
    rmsd: RMSD value
    """
    L = template_coords.shape[0]
    return np.sqrt(np.sum((native_coords - template_coords) ** 2) / L)

def kabsch(native_coords, template_coords):
    """
    kabsch superposition
    parameters:
    native_coords: (N, 3)
    template_coords: (N, 3)
    return:
    - R: rotation matrix (3, 3)
    - translation: translation vector (3,)
    """
    centroid_native_coords = np.mean(native_coords, axis=0)
    centroid_template_coords = np.mean(template_coords, axis=0)
    native_coords_centered = native_coords - centroid_native_coords
    template_coords_centered = template_coords - centroid_template_coords

    H = template_coords_centered.T @ native_coords_centered
    U, S, Vt = np.linalg.svd(H)

    R = U @ Vt
    if np.linalg.det(R) < 0:
        Vt[-1, :] = -Vt[-1, :]
        R = U @ Vt
    translation = centroid_native_coords - centroid_template_coords @ R

    return R, translation

def read_pdb(file_path):
    """
    read pdb file, return atomic coordinates of all CA
    parameters:
    file_path: pdb file path
    return:
    coords: atomic coordinates of all CA (N, 3)
    """
    coords = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                # read atomic coordinates (x, y, z) of all CA
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                coords.append([x, y, z])
    return np.array(coords)

# calculate TM-score
def calculate_tm_score(native_coords, template_coords, d0):
    """
    calculate the TM-score between two sets of atomic coordinates
    parameters:
    native_coords: (N, 3)
    template_coords: (N, 3)
    return:
    tm_score
    """
    L_native = len(native_coords) 
    
    # calculate the distance between corresponding residue pairs
    distances = calculate_distance_matrix(native_coords, template_coords)
    
    # calculate TM-score
    tm_score = np.sum(1 / (1 + (distances / d0) ** 2)) / L_native
    
    return tm_score

def iterative_alignment(native_coords, template_coords, L_values, d0, d0_search, max_iterations=20):
    """ 
    perform heuristic iteration to select the optimal TM-score
    parameters:
    - native_coords: (N, 3)
    - template_coords: (N, 3)
    - L_values: list of initial fragment lengths
    return:
    - rmsd: global RMSD
    - best_rotation: rotation matrix when get best TM-score
    - best_translation: translation vector when get best TM-score
    - best_tm_score: best TM-score
    - final_rmsd: RMSD when get best TM-score
    """
    best_tm_score = -np.inf
    best_rotation = None
    best_translation = None
    rmsd = None
    final_rmsd = None

    N = len(template_coords)
    for L in L_values:
        print(f"Running iterations with initial fragment length L={L}")
        # sliding fragment from N-terminus to C-terminus
        for start in range(0, N - L + 1): 
            # select initial fragment
            native_fragment = native_coords[start:start + L]
            template_fragment = template_coords[start:start + L]

            # calculate rotation matrix
            R, t = kabsch(native_fragment, template_fragment)
            rotated_coords = template_coords @ R + t

            # calculate global RMSD
            if L == L_values[0]:
                rmsd = calculate_rmsd(native_coords, rotated_coords)

            # calculate TM-score
            tm = calculate_tm_score(native_coords, rotated_coords, d0)

            # if get higher TM-score, update rotation matrix and translation vector
            if tm > best_tm_score:
                best_tm_score = tm
                best_rotation, best_translation = R, t
                final_rmsd = calculate_rmsd(native_coords, rotated_coords)

            # select residue pairs with distances less than d0_init
            d0_init = d0_search - 1
            distances = np.linalg.norm(rotated_coords - native_coords, axis=1)
            selected_indices = distances < d0_init

            # update selected residues
            template_fragment = template_coords[selected_indices]
            native_fragment = native_coords[selected_indices]

            while len(template_fragment) < 3 and N > 3:
                d0_init += 0.5
                selected_indices = distances < d0_init
                template_fragment = template_coords[selected_indices]
                native_fragment = native_coords[selected_indices]

            converged = False
            previous_selected_indices = selected_indices

            for it in range(max_iterations):
                # modify d0_init
                d0_init = d0_search + 1

                # recalculate rotation matrix
                R, t = kabsch(native_fragment, template_fragment)
                rotated_coords = template_coords @ R + t

                # calculate TM-score
                tm = calculate_tm_score(native_coords, rotated_coords, d0)

                # if get higher TM-score, update rotation matrix and translation vector
                if tm > best_tm_score:
                    best_tm_score = tm
                    best_rotation, best_translation = R, t
                    final_rmsd = calculate_rmsd(native_coords, rotated_coords)

                # select residue pairs with distances less than d0_init
                distances = np.linalg.norm(rotated_coords - native_coords, axis=1)
                selected_indices = distances < d0_init

                # update selected residues
                template_fragment = template_coords[selected_indices]
                native_fragment = native_coords[selected_indices]

                while len(template_fragment) < 3 and N > 3:
                    d0_init += 0.5
                    selected_indices = distances < d0_init
                    template_fragment = template_coords[selected_indices]
                    native_fragment = native_coords[selected_indices]

                # check whether converged
                if previous_selected_indices is not None and np.array_equal(selected_indices, previous_selected_indices):
                    converged = True
                    break

                previous_selected_indices = selected_indices

            if converged:
                print(f"Converged for L={L}, start={start}, iteration={it+1}")
            else:
                print(f"Max iterations reached for L={L}, start={start}")

    return rmsd, best_rotation, best_translation, best_tm_score, final_rmsd

if __name__ == '__main__':

    # load pdb files
    # template_structure --> native_structure
    args = parse_args()
    template_pdb = args.template_path
    native_pdb = args.native_path
    template_coords = read_pdb(template_pdb)
    native_coords = read_pdb(native_pdb)

    LT = template_coords.shape[0]
    LN = native_coords.shape[0]
    n_init_max = 6
    L_ini_min = min(4, LT)

    L_ini = [max(LT // (2**i), L_ini_min) for i in range(n_init_max - 1)]
    L_ini.append(L_ini_min) # make sure L_ini_min in list
    L_ini = sorted(list(set(L_ini)))[::-1] # unique

    # TM-score: d0
    d0_min = 0.5
    if LN > 15:
        d0 = max(1.24 * (LN - 15) ** (1/3) - 1.8, d0_min)
    else:
        d0 = d0_min
    
    d0_search = d0
    if d0_search > 8:
        d0_search = 8
    if d0_search < 4.5:
        d0_search = 4.5
    
    # calculate TM-score, rmsd, rotation matrix and translation vector
    rmsd, best_R, best_t, best_score, final_rmsd = iterative_alignment(native_coords, template_coords, L_ini, d0=d0, d0_search=d0_search)
    
    print("---------- Results ----------")
    print(f"RMSD: {rmsd}")
    print(f"Best TM-score: {best_score}")
    print(f"Final RMSD: {final_rmsd}")
    print("Best Rotation Matrix:\n", best_R)
    print("Best Translation Vector:\n", best_t)

    # ---------- save rotated coordinates ----------
    print("-----------------------------")
    output_file = f"{template_pdb.split('.')[0]}_superpos.pdb"

    with open(template_pdb, 'r') as f:
        lines = f.readlines()

    with open(output_file, 'w') as f:
        for line in lines:
            if line.startswith("ATOM"):
                # read coordinates
                coords = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                # make rotation and translation
                rotated_coords = coords @ best_R + best_t
                # update coordinates
                new_line = (
                    line[:30] +
                    f"{rotated_coords[0]:8.3f}" +
                    f"{rotated_coords[1]:8.3f}" +
                    f"{rotated_coords[2]:8.3f}" +
                    line[54:])
                f.write(new_line)
            else:
                f.write(line)

    print(f"Rotated coordinates saved to {output_file}.")
