#!/usr/bin/env python3

import numpy as np
import pandas as pd
from periodictable import elements
import argparse
from argparse import RawTextHelpFormatter
from gaussianutility.utilities import readinput

def parse_args():
    parser = argparse.ArgumentParser(
        description= "Sort geometry by an index; 'x', 'y', or 'z' coordinate, atomic number,\n"
                     "or ONIOM layer from low to high in the Gaussian input file (.com or .gjf)\n"
                     "Sorting by an atomic number is defualt"
                     "Sorting in ascending order is default\n\n"
                     "Return file_name.com: Gaussian input file with a sorted geometry",
        formatter_class=RawTextHelpFormatter)

    parser.add_argument('file_name', help='Gaussian input file (.com or .gjf)')
    parser.add_argument('-s', '--sort', nargs='?', const=1, \
    help="Index for sorting;\n"+\
         "'x' or 'y' or 'z' for coordinate, 'A' for atomic number, and 'L' for ONIOM layer\n"+\
         "'L' can be combined with another index\n"+\
         "Ex) 'Lx' or 'xL' combines sorting by ONIOM layer and x coordinate in each layer\n"+\
         "'A' is default for non-ONIOM input and 'AL' is default for ONIOM input", default='AL')

    parser.add_argument('-o', '--order', nargs='?', const=1, \
    help="Index for the order for sorting;\n"+\
         "'a' for ascending order (default) and 'r' for reversed (descending) order\n"+\
         "Only one index can be provided\n"+\
         "If two sorting indeices (after -s) are provided for an ONIOM input,\n"+\
         "the ONIOM layer will be sorted from low to high automatically,\n"+\
         "and the atoms in each layer will be sorted based on the other index", default='a')

    args = parser.parse_args()
    return args
    
def sort_input(file_name, sortIdx, orderIdx):
    # Read a Gaussian input file
    route, title, charge_mult, df_geom, connectivity = readinput(file_name)

    oniomIdx = "oniom" in route.lower()
    
    # Check sorting and ordering indecies
    sortOptions = ['x', 'y', 'z', 'A', 'L', 'xL', 'Lx', 'yL', 'Ly', 'zL',\
    'Lz', 'AL', 'LA']

    if sortIdx not in sortOptions:
        raise ValueError("Unrecognized sorting index is provided")

    sortIdx = list(sortIdx)

    if not oniomIdx and 'L' in sortIdx:
        raise ValueError("The provided input structure does not have ONIOM partitions")
        
    if orderIdx not in ('a', 'r'):
        raise ValueError("Order index must be 'a for ascending or r for reversed (descending) order")
        
    # Set ordering index
    if len(sortIdx) == 1:
        if orderIdx == 'a': orderIdx = [1]
        elif orderIdx == 'r': orderIdx = [0]
        
    elif len(sortIdx) == 2:
        if orderIdx == 'a': orderIdx = [1, 1]
        elif orderIdx == 'r': orderIdx = [1, 0]

    # Sorting atoms
    ## Sorting based on the ONIOM layer first
    if oniomIdx and 'L' in sortIdx:
        df_layer_sort = pd.DataFrame({'ONIOM_layer': ['L', 'M', 'H']})
        layer_sort_mapping = df_layer_sort.reset_index().set_index('ONIOM_layer')
        df_geom['ONIOM_layer_num'] = df_geom['ONIOM_layer'].map(layer_sort_mapping['index'])
    
    ## If sorting only based on the ONION later
    if sortIdx == ['L']:
        df_geom = df_geom.sort_values(by=['ONIOM_layer_num'], ascending = orderIdx)
        df_geom = df_geom.drop(columns=['ONIOM_layer_num'])

    ## If sorting based on the atomic number
    elif 'A' in sortIdx:
        elem = np.array([[el.symbol, el.number] for el in elements])
        elem = {elem[i,0]: int(elem[i,1]) for i in range(len(elem))}
        df_geom['Atomic_num'] = [elem[atomSb] for atomSb in np.array(df_geom['Atom'])]

        if 'L' in sortIdx:
            df_geom = df_geom.sort_values(by=['ONIOM_layer_num', 'Atomic_num'], ascending = orderIdx)
            df_geom = df_geom.drop(columns=['ONIOM_layer_num', 'Atomic_num'])
        else:
            df_geom = df_geom.sort_values(by=['Atomic_num'], ascending = orderIdx)
            df_geom = df_geom.drop(columns=['Atomic_num'])
            
    ## If sorting based on the cartesian coordinate
    else:
        if 'x' in sortIdx: sortCoord = 'x'
        elif 'y' in sortIdx: sortCoord = 'y'
        elif 'z' in sortIdx: sortCoord = 'z'
            
        df_geom[sortCoord] = [float(val) for val in df_geom[sortCoord]]
             
        if 'L' in sortIdx:
            df_geom = df_geom.sort_values(by=['ONIOM_layer_num', sortCoord], ascending = orderIdx)
            df_geom = df_geom.drop(columns=['ONIOM_layer_num'])
        else:
            df_geom = df_geom.sort_values(by=sortCoord, ascending = orderIdx)
            
        df_geom[sortCoord] = ["{:.8f}".format(val) for val in df_geom[sortCoord]]

    # Modify connectivity data
    if 'connectivity' in route:
        atomIdxSwitch = np.vstack([df_geom.index.values, np.array(range(len(df_geom)))]).T
        atomIdxSwitch = atomIdxSwitch + np.ones(atomIdxSwitch.shape, dtype=int)
        atomIdxSwitch = {atomIdxSwitch[i,0]: atomIdxSwitch[i,1] for i in range(len(df_geom))}

        newConnect = []
        for line in connectivity:
            vals = line.split()
            for idx, val in enumerate(vals):
                try: vals[idx] = str(atomIdxSwitch[int(val)])
                except ValueError:
                    pass

            newConnect.append(vals)

        idx2remove = []
        connect2add = []
        for idx1, line in enumerate(newConnect):            
            if len(line) > 1:
                #1, 3, 5, ...
                for idx2 in range(len(line)):
                    if idx2 % 2 == 1:
                        connect2add.append([line[0], line[idx2], line[idx2+1]])
                                          
                idx2remove.append(idx1)
                                          
        for idx in sorted(idx2remove, reverse=True):
            del newConnect[idx]
                    
        newConnect = newConnect + connect2add

        for idx, line in enumerate(newConnect):
            if len(line) > 1 and int(line[0]) > int(line[1]):
                newConnect[idx] = [line[1], line[0], line[2]]
                
        missingAtom = np.setdiff1d(np.array(range(len(df_geom)))+1, np.unique([int(elem[0]) for elem in newConnect]))
        for elem in missingAtom:
            newConnect.append([str(elem)])
                
        newConnect = sorted(newConnect, key=lambda x:int(x[0]))
        newConnectSorted = []

        for i in range(len(df_geom)):
            sublist = [elem for elem in newConnect if int(elem[0]) == i+1]
            if len(sublist) > 1:
                starting = sublist[0][0]
                connected = []
                for elem in sublist:
                    connected.append(elem[1:])
                
                connected = [x for x in connected if x]
                
                if len(connected) > 1:
                    connected = [[starting]] + sorted(connected, key=lambda x:int(x[0]))
                else:
                    connected = [[starting]] + connected
            
                connected = [item for sublist in connected for item in sublist]
                newConnectSorted.append(connected)
                
            else:
                newConnectSorted.append(sublist[0])

    # Write Gaussian input file
    with open(file_name, 'w') as output:
        output.write(f"{route}\n")
        output.write(f"{title}\n")
        output.write(charge_mult)
        output.write(df_geom.to_csv(index=False, header=False, sep='\t'))
        output.write('\n\n')
        if 'newConnectSorted' in locals():
            for line in newConnectSorted:
                output.write(' '.join(line) + '\n')
            output.write('\n')

def main():
    args = parse_args()
    sort_input(args.file_name, args.sort, args.order)

if __name__ == "__main__":
    main()
    
