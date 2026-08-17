import abtem
import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.cluster import wulff_construction
abtem.config.set({"precision": "float32"})

def Zoneaxis4(family_plane):
    '''
    Zoneaxis: A function that find the rotation angle around the x and y-axis to make an abitrary zoneaxis face positive z.
    
    Input = The zoneaxis/facet of intrest in the form family_plane = [h,k,l].
    
    Output = The rotationangle (in degree not radians) around the x and y-axis gives respectfully as theta and phi.

    By: Johannes Varnes (s214475@dtu.dk)

    Honorable mentions: David Hart and Aleksander Karoli Christensen assisted with math.
    '''    
    # It is Assumed that the unit cell is cubic and that the miller indicies are non negative.
    # Then the family plane vector is normalised to a unit vector.
    v1 = family_plane / np.linalg.norm(family_plane)
    
    if v1[2] < 0:
        
        # theta is the rotation angle around the positive x-axis and will rotate the 3D vector into the xz-plane
        if v1[1] == 0:
            theta = 0.0
        else:
            theta = round(np.degrees(np.arctan(v1[1]/v1[2])),3) + 180
    
        # phi is the rotation angle around the negative y-axis, hence the (-), and rotates the 2D vector into the z-axis. 
        phi = -round(np.degrees(np.arcsin(v1[0])),3)
    else:
        
        # theta is the rotation angle around the positive x-axis and will rotate the 3D vector into the xz-plane
        if v1[1] == 0:
            theta = 0.0
        else:
            theta = round(np.degrees(np.arctan(v1[1]/v1[2])),3)
    
        # phi is the rotation angle around the negative y-axis, hence the (-), and rotates the 2D vector into the z-axis. 
        phi = -round(np.degrees(np.arcsin(v1[0])),3)

    if v1[0] == 0 and v1[1] == 0 and v1[2] == -1:
        phi = -180

    return (theta, phi)

def find_ontop_index(atom):
    '''
    find_ontop_value: A function that takes an ase-model and returns the index og all surface-atom on the facet facing positive z.
    
    Input = ase model from wulff_construction().
    
    Output = nunpy array ( array[(index1, index2 ..)] ) of all surface-atom index [1...N_atoms] facing positive z.

    By: Johannes Varnes (s214475@dtu.dk)
    '''
    
    # find the max z-value of the array of coordinates.
    max_value = np.max(atom.get_positions()[:,2])
    
    # making a sorting array of atom coords in the wolf structure. 
    cords_arr = atom.get_positions()
    index_arr = []
    
    # go through each index in the array of atom coordinates.
    for n, cords in enumerate(cords_arr[:,2]):
        
        # if the value is equal to the higest z-value plus minus 1/4 of the latice constant, set the value to True, otherwise False.
        if cords <= max_value and cords > max_value - lc/4:
            index_arr.append(n)
        else:
            continue
        
    return (index_arr)

def facet_sort(facet_list):
    '''
    facet_sort: A function that takes an array of faccets and sorts them after the absolut sum of (h,k,l). Lowest first.
    
    Input = numpy array of faccets where a faccet = (h,k,l).
    
    Output = a sorted numpy array, where the first faccets are {100} and the last faccets are the ones with the higest absolut sum of (h,k,l).

    By: Johannes Varnes (s214475@dtu.dk)
    '''
    # Creating a list for the absolute norm of the faccets.
    abs_hkl = []

    # Go through the list of faccets and sort them based on absolute norm using the function "np.argsort()".
    for f in facet_list:
        abs_hkl.append(np.linalg.norm(f))
    i_sort = np.argsort(abs_hkl)
    
    return (facet_list[i_sort])

def find_surface_index(atom):
    '''
    find_surface_index: A function that takes an ase-model and finds all the surface atoms of an ase wulff-construction model.
    
    Input = ase model from wulff_construction().
    
    Output = numpy matrix of [zoneaxis1 ...etc , index_of_surface_atoms ...etc], zonesaxis = [h,k,l], index_of_surface_atoms = numpy array of index's from find_ontop_values().

    By: Johannes Varnes (s214475@dtu.dk)
    '''
    # Find all relevant faccets for Pt based off http://crystalium.materialsvirtuallab.org/ using the function "model.get_surfaces()" from ASE.
    s = [(1, 0, 0), (1, 1, 0), (1, 1, 1), (3,3,2), (3,2,2), (2,2,1), (3,3,1), (2,1,1), (3,2,1), (3,1,1), (3,1,0), (2,1,0), (3,2,0)]
    e = [0.116, 0.117, 0.093, 0.097, 0.099, 0.100, 0.106, 0.110, 0.110, 0.112, 0.117, 0.118, 0.118]
    small = 50
    model = wulff_construction('Pt', s, e, small, structure = 'fcc', rounding='closest', latticeconstant=lc)
    facets = model.get_surfaces()
    
    on_top_indes = []

    # Go through all faccets systematically in the ase-model using the function "facet_sort".
    for facet in facet_sort(facets):
        
        # Find rotation angles for the n'th facet using function "Zoneaxis".
        theta, phi = Zoneaxis4(facet)
        
        # Rotate the model facet to face positiv z axis 
        atom.rotate(theta, (1, 0, 0), center=(0, 0, 0), rotate_cell=False)
        atom.rotate(phi, (0, 1, 0), center=(0, 0, 0), rotate_cell=False)

        # append all the index's for on top position of this facet.
        new_cord = (facet, find_ontop_index(atom))
        on_top_indes.append(new_cord)

        # Rotate the model back to start position.
        atom.rotate(-phi, (0, 1, 0), center=(0, 0, 0), rotate_cell=False)
        atom.rotate(-theta, (1, 0, 0), center=(0, 0, 0), rotate_cell=False)

    return (on_top_indes)

def sort_surface_index(atom):
    '''
    sort_surface_index: A function that takes an ase-model and distribute surface index atoms to the respective facets.
    
    Observe: only use facets that are defined in the ase wulff_construction() model.
    
    Input = ase model from wulff_construction().
    
    Output = numpy matrix of [facet , facet_list], facet = [h,k,l], facet_list = facet sorted index atoms from find_surface_index() sorted by facet_sort().

    By: Johannes Varnes (s214475@dtu.dk)
    '''

    new_atom = atom.copy()
    
    # Find all the coordinates for on top positions of the atom.
    List_of_index = find_surface_index(new_atom)

    # Create a list of atom index's that is already occupied with a "CO" adsorbant.
    l = []

    # A new list of surface index's with the same format.
    new_list = []
    
    for f, facet in enumerate(List_of_index):
        # Find rotation angles for the n'th facet from the list using function "Zoneaxis".
        theta, phi = Zoneaxis4(facet[0])


        # Rotate the n'th facet to face positiv z axis.
        new_atom.rotate(theta, (1, 0, 0), center=(0, 0, 0), rotate_cell=False)
        new_atom.rotate(phi, (0, 1, 0), center=(0, 0, 0), rotate_cell=False)
        
        # Save all on top coordinates of the atom.
        facet_list = []
        
        for index in List_of_index[f][1]:
            
            # If the coords (atom index) already has a "CO" adsorbant skip to next cords.
            if index in l:
                continue

            # Add the new index to the list of index atoms in this facet.
            facet_list.append(index)
            
            # Add the new index to the list of occupied index atoms.
            l.append(index)

        # Append the new list of indes atoms for each facet. 
        new = (facet[0], facet_list)
        new_list.append(new)
        
        # Rotate the model back to start position.
        new_atom.rotate(-phi, (0, 1, 0), center=(0, 0, 0), rotate_cell=False)
        new_atom.rotate(-theta, (1, 0, 0), center=(0, 0, 0), rotate_cell=False)
        
    return (new_list)

def coverage_list(atom, coverage, facet_list=None, seed=None):
    '''
    coverage_list: A function that takes an ase model and compute the the surface index atoms that will contain an "On top" CO adsorbate from a specified coverage.
    
    Input = atom = ase model from wulff_construction(), coverage = an int64 betweem 0-1, seed = numpy seed randomiser which by default is random. 
    
    Output = a sorted numpy array, where the first faccets are {100} and the last faccets are the ones with the higest absolut sum of (h,k,l).

    By: Johannes Varnes (s214475@dtu.dk)
    '''
    if seed is not None: 
        np.random.seed(seed)
        
    new_atom = atom.copy()
    
    # Find all the coordinates for on top positions of the atom.
    full_list = sort_surface_index(new_atom)

    if facet_list is not None: 
        full_list = facet_list

    # New list for facets and index's with randome coverage.
    cov_list = []

    # Go through the list of faccets and sort them based on absolute norm using the function "np.argsort()".
    for i in range(len(full_list)):
        # N, list of all index accessible in this facets.
        N = full_list[i][1]

        # Randomly picks a procent of the accessible index's rounded down to whole integer based on coverage.
        cov = int(len(N)*coverage)
        #cov_N = np.random.permutation(len(N))[:c]
        
        cov_N = np.random.choice(N, size = cov, replace = False) 
        
        new = (full_list[i][0], cov_N)
        cov_list.append(new)
    
    return (cov_list)

def place_ontop(atom, coverage = 1, facet_list=None, seed=None):
    '''
    place_ontop: A function that takes an ase-model and place "CO" molecules on all "on top" locations of the models facets.
    
    Observe: only use facets that are defined in the ase wulff_construction() model, and only works for nonrotated ase models. Rotation can be performed afterwards.
    
    Input = ase model from wulff_construction(), coverage = an int64 betweem 0-1 by default is 1, seed = numpy seed randomiser by default is random. 
    
    Output = ase model from wulff_construction() with "CO" adsorbant on all "On top" positions.

    By: Johannes Varnes (s214475@dtu.dk)
    '''

    new_atom = atom.copy()
    
    # Find all the coordinates for on top positions of the atom.
    List_of_index = coverage_list(new_atom, facet_list=facet_list, coverage=coverage, seed=seed)
    
    d1 = (1.867+1.877+1.879)/3 # distance from atom to "C" in "CO". Data from Andreas DFT faccet 100, 110, 111 in Å.
    d2 = (1.152+1.153+1.152)/3 # distance from "C" to "O" in "CO". Data from Andreas DFT faccet 100, 110, 111 in Å.
    
    for f, facet in enumerate(List_of_index):
        # Find rotation angles for the n'th facet from the list using function "Zoneaxis".
        theta, phi = Zoneaxis4(facet[0])
        
        # Rotate the n'th facet to face positiv z axis.
        new_atom.rotate(theta, (1, 0, 0), center=(0, 0, 0), rotate_cell=False)
        new_atom.rotate(phi, (0, 1, 0), center=(0, 0, 0), rotate_cell=False)

        # Find all the current coords for the wullf constructions index atoms. (Be aware, the index coords change with rotation).
        index_coords = new_atom.get_positions()
        
        # Place "CO" molecule on top of all on top coordinates of the atom.
        for cords in List_of_index[f][1]:

            # Finds and place "CO" adsorbant on the "On Top" position of the n'th surface atom.
            C_cords = index_coords[cords] + np.array([0,0,d1])
            O_cords = index_coords[cords] + np.array([0,0,d1+d2])
            substrate = Atoms('CO', positions = [C_cords, O_cords])
            new_atom = new_atom + substrate
    
        # Rotate the model back to start position.
        new_atom.rotate(-phi, (0, 1, 0), center=(0, 0, 0), rotate_cell=False)
        new_atom.rotate(-theta, (1, 0, 0), center=(0, 0, 0), rotate_cell=False)
        
    return (new_atom)

def find_neighbours(atom, facet_list=None):
    new_atom = atom.copy()
    
    # Create a list of sorted facets and their surface atoms using the function" sort_surface_index"
    List_of_index = sort_surface_index(new_atom)

    if facet_list is not None: 
        List_of_index = facet_list

    # Create a list of atom index's that is already occupied with a "CO" adsorbant.
    l = []

    # Create a list for defining the neighbour distances for each facet.
    new_list = []
    
    for f, facet in enumerate(List_of_index):
        # Fail safe, Not all facets has surface atoms, in which case skip this facet.
        if len(List_of_index[f][1]) == 0:
            continue
            
        # Find rotation angles for the n'th facet from the list using function "Zoneaxis".
        theta, phi = Zoneaxis4(facet[0])
        
        # Rotate the n'th facet to face positiv z axis.
        new_atom.rotate(theta, (1, 0, 0), center=(0, 0, 0), rotate_cell=False)
        new_atom.rotate(phi, (0, 1, 0), center=(0, 0, 0), rotate_cell=False)

        # Find all the current coords for the wullf constructions index atoms. (Be aware, the index coords change with rotation).
        index_coords = new_atom.get_positions()

        # Create a refference atom using the first surface index atom in the list from surface f.
        ref = index_coords[List_of_index[f][1][0]]

        # Create a list atom distances from each surface atom in respect to the refference atom.
        d_list = []
        
        # Create a list of distances from refference atom, to determine bridg, long bridg and hollow site.
        for cords in List_of_index[f][1]:
            new_cord = index_coords[cords]
            dx = new_cord[0] - ref[0]  # x-coordinate
            dy = new_cord[1] - ref[1]  # y-coordinate
            d = np.sqrt(dx**2 + dy**2) # distance between the 2 atoms in 2D (No z-coordinate)
            d_list.append(np.round(d,6))

        # Sort the list of distances using np.argsort
        d_list = np.array(d_list) # from list to array
        i_sort = np.argsort(d_list) # Create mask for sorting low to high 
        d_list_sorted = d_list[i_sort] # Apply mask (both must be np.arrays)
        d_list_sorted = list(dict.fromkeys(d_list_sorted)) # Remove duplicates
        d_list_sorted = [d for d in d_list_sorted if d > 0 and d < (1.99*d_list_sorted[1])] # removes 0 and 2nd order neightbours.
        
        # define distances for, bridg, long bridg and hollow using d_list_sorted.
        def_list = np.zeros(3)
        def_list[:len(d_list_sorted)] = d_list_sorted # exchange 0 with definition
        bridge = def_list[0] # The closest neighbour is a bridge
        hollow = def_list[1] # The next closest should be hollow
        long_bridge = def_list[2] # The one furthest away is long bridge

        # Append the new list of indes atoms for each facet. 
        new = (facet[0], facet[1], def_list)
        new_list.append(new)
        
        # Rotate the model back to start position.
        new_atom.rotate(-phi, (0, 1, 0), center=(0, 0, 0), rotate_cell=False)
        new_atom.rotate(-theta, (1, 0, 0), center=(0, 0, 0), rotate_cell=False)
        
    return (new_list)

def place_bridge(atom, facet_list=None):
    new_atom = atom.copy()
    #List_of_index = test
    List_of_index = find_neighbours(new_atom, facet_list)

    # If the Co coord already exist, skip to next iterated surface atom
    def cord_pair_in_list(pair, lst, tol=1e-4):
        return any(all(np.allclose(p1, p2, atol=tol) for p1, p2 in zip(pair, item)) for item in lst)

    # Define adsorbate distances surface atom.
    d1 = 1.48 # 1.48 nm from Pt to C
    d2 = 1.17 # 1.17 nm from C to O 


    for f, facet in enumerate(List_of_index):
        # print(f, facet)
        # Fail safe, Not all facets has surface atoms, in which case skip this facet.
        if len(List_of_index[f][1]) == 0:
            continue
            
        # Create a list of atom index's that is already occupied with a "CO" adsorbant.
        l = []
            
        # Find rotation angles for the n'th facet from the list using function "Zoneaxis".
        theta, phi = Zoneaxis4(facet[0])
        
        # Rotate the n'th facet to face positiv z axis.
        new_atom.rotate(theta, (1, 0, 0), center=(0, 0, 0), rotate_cell=False)
        new_atom.rotate(phi, (0, 1, 0), center=(0, 0, 0), rotate_cell=False)
        
        # Find all the current coords for the wullf constructions index atoms. (Be aware, the index coords change with rotation).
        index_coords = new_atom.get_positions()
        
        # Place "CO" molecule on top of all on top coordinates of the atom.
        for surf_atom in List_of_index[f][1]:
            # Get coords of iterate atom used to check for bridg neighbours
            iter_atom = index_coords[surf_atom]
            
            # Loop over all other surface atoms to check for bridge neighbour 
            for neighbour in List_of_index[f][1]:
                new_cord = index_coords[neighbour]
                dx = new_cord[0] - iter_atom[0]  # x-coordinate
                dy = new_cord[1] - iter_atom[1]  # y-coordinate
                d = np.sqrt(dx**2 + dy**2) # distance between the 2 atoms in 2D (No z-coordinate)
                d = np.round(d,6) # Round to 6 ciffers
                #print(f, surf_atom, d)

                tol = 1e-4  # tolerance value for float comparison
                bridg = List_of_index[f][2][0]
                long_bridge = List_of_index[f][2][2]

                if (np.isclose(d, bridg, atol=tol) or np.isclose(d, long_bridge, atol=tol)) and d != 0:
                    # print(f' bridge:{List_of_index[f][2][0]}', f'lbridge:{List_of_index[f][2][2]}', f'iterated atom:{d}')
                    # Finds and place "CO" adsorbant on the "bridge" position of the n'th surface atom.
                    #x_shift = iter_atom[0] +  # bridg x shift from iterated atom
                    #y_shift = iter_atom[1] + dy # bridg y shift from iterated atom
                    C_cords = iter_atom + np.array([dx/2,dy/2,d1])
                    O_cords = iter_atom +  np.array([dx/2,dy/2,d1+d2])
                    new_cord = [C_cords, O_cords]
                    # print(f' Coordinates:{new_cord}', f'original{iter_atom}')
                    
                    # This checks if all arrays in new_cord are present in l
                    if cord_pair_in_list(new_cord, l):
                        continue

                    substrate = Atoms('CO', positions = new_cord)
                    new_atom = new_atom + substrate
                    l.append(new_cord)
                    # print('2', new_cord, surf_atom)
                    

                    # flip = new_atom.copy()
                    # flip.rotate(180, (1, 0, 0), center=(0, 0, 0), rotate_cell=True)
                    # abtem.show_atoms(flip, scale=1)
        
        # Rotate the model back to start position.
        new_atom.rotate(-phi, (0, 1, 0), center=(0, 0, 0), rotate_cell=False)
        new_atom.rotate(-theta, (1, 0, 0), center=(0, 0, 0), rotate_cell=False)
        # print(len(l))

    return (new_atom)

def stable_solution(atom):
    new_atom = atom.copy()
    facets = sort_surface_index(new_atom)

    bridg_list = []
    ontop_list = []
    B = []
    T = []

    for facet in facets:
        if np.sum(abs(facet[0])) == 1:
            bridg_list.append(facet)
            B.append(facet[1])
        else:
            ontop_list.append(facet)
            T.append(facet[1])
        
    new_atom = place_bridge(new_atom, facet_list=bridg_list)
    new_atom = place_ontop(new_atom, coverage=1, facet_list=ontop_list)
    
    return(new_atom)