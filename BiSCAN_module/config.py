from ase.io import read
import pandas as pd
from ase.neighborlist import NeighborList, natural_cutoffs
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import re
import pathlib
import dash
from dash import dcc, html, Input, Output

## Default properties for each atom.
default_properties = {
    "H": {"colour": "white", "radius": 10.2652},
    "He": {"colour": "blue", "radius": 10.2652},
    "Li": {"colour": "yellow", "radius": 30.2652},
    "Be": {"colour": "green", "radius": 30.2652},
    "B": {"colour": "green", "radius": 30.2652},
    "C": {"colour": "grey", "radius": 30.4652},
    "N": {"colour": "blue", "radius": 30.4652},
    "O": {"colour": "red", "radius": 30.4652}
    "F": {"colour": "green", "radius": 30.2652},
    "Ne": {"colour": "yellow", "radius": 30.2652},
    "Na": {"colour": "orange", "radius": 30.2652},
    "Mg": {"colour": "orange", "radius": 30.2652},
    "Al": {"colour": "grey", "radius": 30.2652},    
}

## Auxiliary function I used to get the geometries.
def is_number(i):
    try:
        float(i)
        return True
    except ValueError:
        return False




class BiSCAN:   
    '''
    Class defined upon the feeding of an output file from Orca5 or Gaussian16 of a relaxed bidimensional scan.
    Instancing an object of this class will automatically generate the bidimensional analysis.
    
    Example: biscan = BiSCAN('output_file.log')
    '''
    def __init__(self,file_string):
        with open(file_string, 'r') as f:
            file = f.read()
        if self.is_orca(file) is True:
            fig1 = self.get_scan_orca(file)
            mol = self.get_coords_orca(file)
        elif self.is_gaussian(file) is True:
            fig1 = self.get_scan_gaussian(file)
            mol = self.get_coords_gaussian(file)
        else:
            print('Either bad termination or not an Orca/Gaussian scan output file.')
            SystemExit
        lower_x, upper_x, lower_y, upper_y, lower_z, upper_z = self.create_box(mol)
        init_fig2 = self.draw_structure(mol[0], lower_x, upper_x, lower_y, upper_y, lower_z, upper_z)
        app = self.start_dash(fig1,init_fig2, mol,lower_x, upper_x, lower_y, upper_y, lower_z, upper_z)
        app.run_server(debug=True)
            
        
    #######################################
     
     
    def get_scan_gaussian(self,file):
        '''
        Method that generates the figure showing the result of the bidimensional scan from Gaussian.
        
        Example: fig = biscan.get_scan_gaussian('output_file.log')
        '''
        ## First find the summary of the scan:
        ini = re.search('Summary of Optimized Potential Surface Scan', file)
        fin = re.search('GradGradGrad', file[ini.start():])
        real_fin = int(ini.start()) + int(fin.start())
        summary = file[ini.start():real_fin]

        ### Find the root of the energy
        pattern_energy = r'-\d+.\d+'
        extra_energy = re.search(pattern_energy,summary)


        data1 = pd.DataFrame(columns=['energies'])
        data2 = pd.DataFrame(columns=['d1'])
        data3 = pd.DataFrame(columns=['d2'])

        # print(summary.strip().splitlines()[2].split())
        for line in summary.strip().splitlines():
            if 'Eigenvalues' in line:
                for string in line.split():
                    if is_number(string) is True:
                        new_data = pd.DataFrame({"energies" : [float(extra_energy.group())-float(string)]})
                        data1 = pd.concat([data1,new_data],ignore_index=True)
            elif 'R1' in line:
                for string in line.split():
                    if is_number(string) is True:
                        new_data = pd.DataFrame({"d1" : [string]})
                        data2 = pd.concat([data2,new_data],ignore_index=True)
            elif 'R2' in line:
                for string in line.split():
                    if is_number(string) is True:
                        new_data = pd.DataFrame({"d2" : [string]})
                        data3 = pd.concat([data3,new_data],ignore_index=True)
        
        data_scan = pd.concat([data1, data2, data3], axis=1)               
        fig1 = go.Scatter3d(
                        x=data_scan['d1'], y=data_scan['d2'], z=data_scan['energies'],
                        mode='markers',
                        marker=dict(size=8, color=data_scan['energies'],colorscale='spectral',),
                        name='Main Points',
                        # opacity=0.7
                    )

        return fig1
     
     
     ########################################
    
    def get_scan_orca(self,file):
        '''
        Method that generates the figure showing the result of the bidimensional scan from Orca.
        
        Example: atoms_object = biscan.get_scan_orca('output_file.log')
        '''

        ini = re.search("The Calculated Surface using the 'Actual Energy'", file)
        end = re.search('The Calculated Surface using the SCF energy', file)

        text = file[ini.end():end.start()]

        d1 = []
        d2 = []
        energies = []
        for line in text.strip().splitlines():
            for index,string in enumerate(line.split()):
                if index == 0:
                    d1.append(float(string))
                if index == 1:
                    d2.append(float(string))
                if index == 2:
                    energies.append(float(string))

        data1 = pd.DataFrame({'d1' : d1})
        data2 = pd.DataFrame({'d2' : d2})
        data3 = pd.DataFrame({'energies' : energies})


        data_scan1 = pd.concat([data1, data2, data3], axis=1)

        data_scan = data_scan1.sort_values('energies')

        fig1 = go.Scatter3d(
                x=data_scan['d1'], y=data_scan['d2'], z=data_scan['energies'],
                mode='markers',
                marker=dict(size=8, color=data_scan['energies'],colorscale='spectral',),
                name='Main Points',
                # opacity=0.7
            )
        return fig1
    
    
    ############################################
    
    
    
    def get_coords_orca(self,file):
        '''
        Method that returns an Atom object (class from ASE) containing the converged structures from an Orca output file.
        
        Example: biscan.get_coord_orca('output_file.log')
        '''
        output_file = 'temporary_file_coord.xyz'
        
        
        ### number of atoms
        match = re.search('Number of atoms', file)
        match_n_at = re.search(r'\d+', file[match.start():])
        n_at = int(match_n_at.group())


        matches_ini = [ match.start() for match in re.finditer('FINAL ENERGY EVALUATION AT THE STATIONARY POINT', file)]
        for match in matches_ini:
            match_end = re.search('CARTESIAN COORDINATES \\(A.U.\\)', file[match:])
            text = file[match:match+match_end.start()]
            len_text = len(text.splitlines())
            with open(output_file, 'a') as f:
                f.write(f'{n_at}\n\n')
                for i in range(6,len_text-2):
                    f.write(f'{text.splitlines()[i]}\n')


        mol = read(output_file, index=":")
        pathlib.Path.unlink(output_file)
        return mol
    #############################################

    def is_orca(self,file):
        '''
        Method that returns True if the output file is from Orca.
        
        Example: is_file_orca = biscan.is_orca('output_file.log')
        '''
        string_end = 'ORCA TERMINATED NORMALLY'
        match = re.search(string_end[::-1], file[::-1])
        if match:
            return True
        

    #############################################


    def get_coords_gaussian(self,file):
        '''
        Method that returns an Atom object (class from ASE) containing the converged structures from a Gaussian output file.
        
        Example: atoms_object = biscan.get_coord_gaussian('output_file.log')
        '''
        output_file = 'temporary_file_coord.xyz'
        ### Number of atoms

        re_charge=re.search('Multiplicity', file)
        re_modre= re.search('The following ModRedundant', file)

        temp = file[re_charge.start():re_modre.start()]
        n_at = len(temp.splitlines())-3

        ### List of atomic symbols

        list_symbols = []
        pattern= r'\w+'
        for i in range(1,len(temp.splitlines())-2):
            match = re.search(pattern,temp.splitlines()[i])
            list_symbols.append(match.group())


        ### Get coords
        matches_conv = [match.start() for match in re.finditer('Stationary point found', file)]

        string_match ='Coordinates'

        matches_coord = []
        for match_conv in matches_conv:
            match = re.search(string_match[::-1], file[:match_conv][::-1])
            match_new = match_conv - match.start()
            matches_coord.append(match_new)



        matches_rot = []
        for match_coord in matches_coord:
            match_rot = re.search('Rotational constants', file[match_coord:])
            match_new = match_rot.start() + match_coord
            matches_rot.append(match_new)
            

        pattern = r'-*\d+.\d+'
        for i,j in zip(matches_coord, matches_rot):
            temp = file[i:j]
            list_x = []
            list_y = []
            list_z = []
            for linea in temp.strip().splitlines()[3:-1]:
                matches = [match for match in re.finditer(pattern, linea)]
                list_x.append(float(matches[0].group()))
                list_y.append(float(matches[1].group()))
                list_z.append(float(matches[2].group()))
            with open(output_file, 'a') as out:
                out.write(f'{n_at}\n')
                out.write(f'\n')
                for i in range(len(list_symbols)):
                    out.write(f'{list_symbols[i]}   {list_x[i]}     {list_y[i]}     {list_z[0]}\n')

        mol = read(output_file, index=":")
        pathlib.Path.unlink(output_file)
        return mol


    ##########################################



    def is_gaussian(self,file):
        '''
        Method that returns True if the output file is from Gaussian.
        
        Example: is_file_gaussian = biscan.is_gaussian('output_file.log')
        '''
        string_end = 'Entering Gaussian System'
        match = re.search(string_end, file)
        if match:
            return True

    #############################################
    def create_box(self,mol):
        '''
        Method that returns suitable figure limits for the structure representation from an Atoms object.
        
        Example: lower_x, upper_x, lower_y, upper_y, lower_z, upper_z = biscan.create_box(Atoms_object)
        '''
        list_x =[]
        list_y =[]
        list_z =[]
        for image in mol:
            for atom in range(len(image.get_positions())):
                list_x.append(image.get_positions()[atom][0])
                list_y.append(image.get_positions()[atom][1])
                list_z.append(image.get_positions()[atom][2])
        
        min_x = np.min(list_x)
        min_y = np.min(list_y)
        min_z = np.min(list_z)


        max_x = np.max(list_x)
        max_y = np.max(list_y)
        max_z = np.max(list_z)


        dist_x = max_x - min_x
        dist_y = max_y - min_y
        dist_z = max_z - min_z

        max_dist = np.max([dist_x,dist_y,dist_z])

        x_center = sum(list_x)/len(list_x)
        y_center = sum(list_y)/len(list_y)
        z_center = sum(list_z)/len(list_z)

        upper_x = x_center + max_dist
        lower_x = x_center - max_dist

        upper_y = y_center + max_dist
        lower_y = y_center - max_dist

        upper_z = z_center + max_dist
        lower_z = z_center - max_dist
        
        return lower_x, upper_x, lower_y, upper_y, lower_z, upper_z
    
    def draw_structure(self, mol, lower_x, upper_x, lower_y, upper_y, lower_z, upper_z):
        '''
        Method that generates a 3D plot (returns a figure) to represent a molecular structures. 
        It needs as input the Atoms object and the figure limits created by the create_box() method.
        
        Example: fig = biscan.draw_structure(atoms_object, lower_x, upper_x, lower_y, upper_y, lower_z, upper_z)
        '''
        data_x = pd.DataFrame(columns=['x'])
        data_y = pd.DataFrame(columns=['y'])
        data_z = pd.DataFrame(columns=['z'])
        data_symbol = pd.DataFrame(columns=['symbol'])
        data_color = pd.DataFrame(columns=['color'])
        data_size = pd.DataFrame(columns=['size'])

        ### COORDS
        for atom in mol.get_positions():
            for index,coord in enumerate(atom):
                if index == 0:
                    if len(data_x) == 0:
                        data_x = pd.DataFrame({'x' : [coord]})
                    else:
                        dato = pd.DataFrame({'x' : [coord]})
                        data_x = pd.concat([data_x,dato], ignore_index=True)
                if index == 1:
                    if len(data_y) == 0:
                        data_y = pd.DataFrame({'y' : [coord]})
                    else:
                        dato = pd.DataFrame({'y' : [coord]})
                        data_y = pd.concat([data_y,dato], ignore_index=True)
                if index == 2:
                    if len(data_z) == 0:
                        data_z = pd.DataFrame({'z' : [coord]})
                    else:
                        dato = pd.DataFrame({'z' : [coord]})
                        data_z = pd.concat([data_z,dato], ignore_index=True)
        ### SYMBOLS
        for atom in mol.get_chemical_symbols():
            dato = pd.DataFrame({'symbol' : [atom]})
            data_symbol = pd.concat([data_symbol,dato], ignore_index=True)
        ### COLORS ACCORDING TO SYMBOL
        lista_colores = []
        for symbol in data_symbol['symbol']:
            if symbol in default_properties:
                lista_colores.append(default_properties[symbol]['colour'])
            else:
                lista_colores.append('orange')
        
        dato = pd.DataFrame({'color' : lista_colores})    
        data_color = pd.concat([data_color,dato], ignore_index=True)


        ### SIZE
        lista_size = []
        for symbol in data_symbol['symbol']:
            if symbol in default_properties:
                lista_size.append(default_properties[symbol]['radius'])
            else:
                lista_size.append(30.4652)
        
        data_size = pd.DataFrame({'size' : lista_size})    
        

        data_coord = pd.concat([data_x,data_y,data_z,data_symbol,data_color,data_size],axis=1)
        fig1 = px.scatter_3d(data_coord, x='x', y='y', z='z', color='color',color_discrete_map='identity', size='size', size_max=10)

        ### BONDS 
        ntrl_ctff = natural_cutoffs(mol)
        nl = NeighborList(cutoffs=ntrl_ctff)
        nl.update(mol)
        con_mat=nl.get_connectivity_matrix().asformat('array')
        pos1 = list(zip(*np.where(con_mat == 1)))


        ## pos contains the list with the bonded atoms
        pos = []
        for i,j in pos1:
            if i == j:
                continue
            pos.append([i,j])


        list_lines = []
        for i,j in pos:
            list_x = []
            list_y = []
            list_z = []
            list_color_bond = []
            list_color_bond.append('black')
            list_color_bond.append('black')
            list_x.append(mol.get_positions()[i][0])
            list_x.append(mol.get_positions()[j][0])
            list_y.append(mol.get_positions()[i][1])
            list_y.append(mol.get_positions()[j][1])
            list_z.append(mol.get_positions()[i][2])
            list_z.append(mol.get_positions()[j][2])
            data_bond = pd.DataFrame({'x' : list_x, 'y' : list_y, 'z' : list_z, 'color' : list_color_bond})
            linea = go.Scatter3d(x=data_bond['x'], y=data_bond['y'], z=data_bond['z'], mode='lines', line=dict(color='black', width=5))
            list_lines.append(linea)

        fig3 = go.Figure(data=fig1.data + tuple(list_lines))


        fig3.update_traces(
        marker=dict(
            size=20),
        line=dict(
            color='black',
            width=15
        )
        )
        fig3.update_layout(showlegend=False)
        fig3.update_layout(
            scene=dict(
                xaxis=dict(range=[lower_x, upper_x]),
                yaxis=dict(range=[lower_y, upper_y]),
                zaxis=dict(range=[lower_z, upper_z])
            )
        )
        
        return fig3
    
    def start_dash(self,fig1,init_fig2, mol,lower_x, upper_x, lower_y, upper_y, lower_z, upper_z):
        '''
        Method that creates an app that turns both figures interactive though a callback function. The input of the callback function is the point you click on,
        while the output is a figure product of the draw_structure() method with the updated geometry of the system corresponding to that point of the scan.
        This method receives as input the figures to be displayed (the scan from get_scan_gaussian() or get_scan_orca()), the geometry representation (from draw_structure()),
        the Atoms object (from get_coords_gaussian() or get_coords_orca()) and the figure limits (from create_box()).
        
        Example: app = biscan.start_dash(fig1, fig2, atoms_object, lower_x, upper_x, lower_y, upper_y, lower_z, upper_z)
        '''
        app = dash.Dash(__name__)

        app.layout = html.Div(children=[
        dcc.Graph(id='scan', figure={'data': [fig1], 'layout' : go.Layout(title='Bidimensional Scan', scene=dict(xaxis=dict(title='d1 / Angstrom'),yaxis=dict(title='d2 / Angstrom'),zaxis=dict(title='Energy / Ha')))}, style={'height' : '950px', 'width' : '50%', 'display' : 'inline-block'}),
                        dcc.Graph(id='geometry', figure={'data': [init_fig2], 'layout' : go.Layout(title='Geometry of selected point', scene=dict(xaxis=dict(title='x / Angstrom'),yaxis=dict(title='y / Angstrom'),zaxis=dict(title='z / Angstrom')))}, style={'height' : '950px', 'width' : '50%', 'display' : 'inline-block'})])
        @app.callback(
        Output('geometry', 'figure'),
        [Input('scan', 'clickData')]
        )
        def update_fig2(clickData):
            if clickData is None:
                return init_fig2

            point_index = clickData['points'][0]['pointNumber']
            image = mol[point_index]
            fig2_new = self.draw_structure(image, lower_x, upper_x, lower_y, upper_y, lower_z, upper_z)
            return fig2_new
        
        return app
    
