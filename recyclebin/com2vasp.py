import click
from gaussianutility.utils import com2VASP

@click.command(name = 'com2vasp',
               context_settings = {'help_option_names': ['-h', '--help'],
                                   'show_default': True})
@click.argument('com_file', type=str)
@click.option('-t', '--title', required=False, help='Title of the input file')
@click.option('-l', '--lattice_param', required=False, nargs = 3, help='List of lattice parameters, a, b, and c')

def main(com_file, title, lattice_param):
    """
    Convert .com file to .vasp file.
    Currently, only support cubin structure, where a,b,c represent lattice parameters.
    """

    com2VASP(com_file, title, lattice_param)
    
## Core part that was in utils is below
#=======================================================================
# def com2VASP(com_file, title = 0, lattice_param = 0):
#     '''
#     Convert .com file to .vasp file.
#     Currently, only support cubin structure, where a,b,c represent lattice parameters.
#     '''
#     if not title:
#         title = com_file
#     if not lattice_param:
#         a, b, c = 20.0168991089, 19.8754005432, 13.3344001770
#     else:
#         [a, b, c] = lattice_param
        
#     informat, route, title_and_spin, df_geom, rest = read_input(com_file)

#     elem = []
#     for el in elements:
#         elem.append(el.symbol)

#     df_atom_sort = pd.DataFrame(elem[1:], columns = ['symbol'])
#     atom_sort_mapping = df_atom_sort.reset_index().set_index('symbol')
#     df_geom['Atomic_num'] = df_geom['Atom'].map(atom_sort_mapping['index'])

#     df_geom = df_geom.sort_values(by=['Atomic_num'], ascending = 0)
#     df_geom = df_geom.drop(columns=['Atomic_num'])
#     df_geom = df_geom.drop(df_geom[df_geom['Atom'] == 'Tv'].index)
#     df_geom = df_geom.drop(df_geom[df_geom['Atom'] == '?'].index)
    
#     atomlist = np.array(df_geom['Atom'])
#     coord = df_geom[['x','y','z']]

#     uniqatomlist = []
#     uniqatomcount = []
#     for elem in atomlist:
#         if elem not in uniqatomlist:
#             uniqatomlist.append(elem)
#             uniqatomcount.append(np.count_nonzero(atomlist == elem))

#     uniqatomlist = np.array(uniqatomlist)
#     uniqatomcount= np.array(uniqatomcount, dtype=str)
    
#     outfilename = title.rsplit(".",1)[0] + ".vasp"
#     outfile = open(outfilename, 'w')
#     outfile.write(title+'\n')
#     outfile.write('1.0\n')
#     outfile.write('  '+str(a)+'  0.0000000000  0.0000000000\n')
#     outfile.write('  0.0000000000  '+str(b)+'  0.0000000000\n')
#     outfile.write('  0.0000000000  0.0000000000  '+str(c)+'\n')
#     outfile.write('   '+' '.join(uniqatomlist)+'\n')
#     outfile.write('   '+' '.join(uniqatomcount)+'\n')
#     outfile.write('Cartesian\n')
#     outfile.write(coord.to_csv(index=False, header=False, sep='\t'))
#     outfile.write('\n')
#     outfile.close()
    
