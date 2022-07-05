#=======================================================================
def freeze_low(file_name):
    """Freeze Cartesian coordinates of atoms in the low ONIOM layer.
    Only .com and .gjf file is supported.
    """
    informat = file_name.rsplit(".",1)[-1]
    if not informat == 'com' or informat == 'gjf':
        raise TypeError('The input file format must be Gaussian input (com or gjf)')
        
    lines = open(file_name).readlines()

    for idx, line in enumerate(lines):
        if "#" in line:
            idx_route = idx
            route = line
            break

    if len(lines[idx_route+1]) >=2:
        route = route.replace('\n', ' ')
        route += lines[idx_route+1]
        idx_route += 1

    if not ("oniom" or "ONIOM" or "Oniom") in route:
        raise AttributeError("Must have ONIOM formulation")

    if not "modredundant" in route:
        lines[idx_route] = lines[idx_route].replace('opt', 'opt=modredundant')

    connect = 0
    if "geom=connectivity" in route:
        connect += 1

    geom = []

    for idx2, line in enumerate(lines[idx_route+5:]):
        geom.append(line.split())
        geomEndIdx = idx_route+idx2+5
        if len(line) <= 2: break

    # if connect == 1:
    #     for id3, line in enumerate(lines[])
    geom = list(filter(None, geom))

    lineLen = 0
    for line in geom:
        if len(line) > lineLen: lineLen = len(line)

    clmnName = []
    for i in range(lineLen):
        name = "Column {}".format(i)
        clmnName.append(name)

    clmnName[5] = 'ONIOM layer'
    df_geom = pd.DataFrame(geom, columns= clmnName)
    oniom_layer = pd.DataFrame(df_geom, columns = ['ONIOM layer'])

    Low_level = np.array(oniom_layer.index[oniom_layer['ONIOM layer']=='L'])+1

    if connect == 0:
        editfile = open(file_name, 'w')
        for line in lines[:geomEndIdx]:
            editfile.write(line)
        editfile.write('\n')
        for elem in Low_level:
            editfile.write('X ' + str(elem) + ' F\n')
        editfile.write('\n')
        editfile.close()

    if connect == 1:
        for idx, line in enumerate(lines[geomEndIdx+1:]):
            if len(line) <= 2:
                connectEndIdx = geomEndIdx+idx+1
                break

        editfile = open(file_name, 'w')
        for line in lines[:connectEndIdx+1]:
            editfile.write(line)
        for elem in Low_level:
            editfile.write('X ' + str(elem) + ' F\n')
        editfile.write('\n')
        editfile.close()

