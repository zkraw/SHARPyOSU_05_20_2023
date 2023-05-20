import h5py as h5
import numpy as np
import pandas as pd
import subprocess
import os

class Aero():
    def __init__(self, route, case, data, airfoil_data, pass_aero, xfoil, excel):
        self.route = route
        self.case = case
        self.airfoil_data = airfoil_data
        self.num_elem = pass_aero.num_elem
        self.num_node = pass_aero.num_node
        self.boundary = pass_aero.boundary
        self.x = pass_aero.x
        self.y = pass_aero.y
        self.z = pass_aero.z
        self.conn = pass_aero.conn
        self.xfoil = xfoil
        
        if excel == True:
            self.excel(data)
        
        if excel == False:    
            self.Aero_Calc(data)
            
        
    def airfoil_inter(self, value, air_data):
        """Genearates new airfoil that is the result of linearly interpolating between two airfoils

        Args:
            value (float): airfoil ID number EX: 1.2 will calculate the airfoil 20% between airfoils 1 and 2 
            air_data (list[float]): the known airfoils used on the planform. n row, 2 columns: x/c then y/c

        Returns:
            list(float): two column list that has x/c and y/c values of the new interpolated value
        """
        begin = int(value)
        end = int(value) + 1

        def interp(middle, x, y):
            y_interp = ((middle-x[0])*(y[1]-y[0])) / (x[1] - x[0]) + y[0]
            return y_interp

        air_begin = air_data[begin]
        air_end = air_data[end]
        
        new_air = []
        for i in range(len(air_begin)):
            for j in range(len(air_end)):
                try:
                    if air_begin[i][0] > air_end[j][0] and air_begin[i][0] < air_end[j+1][0]:
                        y_interp = interp(air_begin[i][0], [air_end[j][0], air_end[j+1][0]], [air_end[j][1], air_end[j+1][1]])
                        nu = (y_interp - air_begin[i][1]) * (value - begin) + air_begin[i][1]
                        new_air.append([air_begin[i][0], nu])
                        break
                    elif air_begin[i][0] == air_end[j][0]:
                        nu = (air_end[j][1] - air_begin[i][1]) * (value - begin) + air_begin[i][1]
                        new_air.append([air_begin[i][0], nu])
                        break
                except:
                    pass
        return new_air
    
    # used to get drag polars for interpolated airfoils
    def xfoil_oper(self, chord):
        x = self.xfoil
        Re = (x['rho'] * x['u_inf'] * chord) / x['u']
        
        if os.path.exists("polar_file.txt"):
            os.remove("polar_file.txt")
    
        XFOIL_PATH = str(os.getcwd()) + '/OSU_Contribution/XFOIL6.99/xfoil.exe'

        xfoil = subprocess.Popen(XFOIL_PATH, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            universal_newlines=True)
        
        actions =  ['load Temp_Interp.dat\n',
        # ["LOAD {0} \n ".format('Temp_Interp.dat'),
                                        "\n",
                                        "PANE \n ",
                                        "OPER \n ",
                                        "Visc {0} \n ".format(Re),
                                        "Mach {0} \n ".format(x['Mach']),
                                        "ITER {0} \n ".format(x['n_iter']),
                                        "PACC \n ",
                                        "polar_file.txt \n\n ",
                                        "aseq {0} {1} {2} \n ". format(x['alfa_1'], x['alfa_2'], x['alfa_step']),
                                        "quit"]
        command = ''.join(actions)
            
        xfoil.communicate(input=command)
    
        # Updating the drag polar for the change in chord
        nu_pol = np.loadtxt("polar_file.txt", skiprows=12)
        new_polar = []
        for i in range(len(nu_pol)):
            new_polar.append([nu_pol[i][0]*np.pi / 180, nu_pol[i][1], nu_pol[i][2], nu_pol[i][4]])
        
        return new_polar
    
    def interp_check(self, check_name, first_dat):
        with open(check_name) as check, open(first_dat) as dat:
            check_lines = check.readlines()
            dat_lines = dat.readlines()
        if dat_lines == check_lines:
            return True
        else:
            return False 
        
    def Aero_Calc(self, data):
        case_name = self.case
        route = self.route
        num_elem = self.num_elem
        num_node = self.num_node
        num_node_elem = data['num_node_elem']
        x = self.x
        y = self.y
        z = self.z
        conn = self.conn
        air_data = self.airfoil_data
        n_airfoils = len(air_data)
        # chords
        chords = np.zeros((num_elem, num_node_elem))

        # twist
        twist = np.zeros((num_elem, num_node_elem))

        # sweep
        sweep = np.zeros((num_elem, num_node_elem))

        # airfoil distribution
        airfoil_distribution = np.zeros((num_elem, num_node_elem), dtype=int)

        # aero_node, aero node definition
        aero_node = np.zeros((num_node,), dtype=bool)

        # elastic_axis
        elastic_axis = np.zeros((num_elem, num_node_elem))

        element_counter = 0
        node_counter = 0
        airfoil_hold = n_airfoils
        added_airfoils = []
        for k in range(len(data['origins'])):
            for i in range(len(data['origins'][k]['branch'])):
                for j in range(len(data['origins'][k]['branch'][i]['section'])):
                    try:
                        aero1 = data['origins'][k]['branch'][i]['section'][j]['aero']
                        try:
                            aero2 = data['origins'][k]['branch'][i]['section'][j]['aero2']
                        except:
                            aero2 = aero1

                        n_elems = data['origins'][k]['branch'][i]['section'][j]['n_elements']

                        chord_temp = np.linspace(aero1[0], aero2[0], n_elems * 2 + 1)
                        twist_temp = np.linspace(aero1[1], aero2[1], n_elems * 2 + 1)
                        sweep_temp = np.linspace(aero1[2], aero2[2], n_elems * 2 + 1)
                        airfoil_temp = np.linspace(aero1[3], aero2[3], n_elems * 2 + 1)
                        elastic_temp = np.linspace(aero1[4], aero2[4], n_elems * 2 + 1)

                        if element_counter == 0:
                            chord_temp = np.linspace(aero1[0], aero2[0], (n_elems * 2) + 2)
                            twist_temp = np.linspace(aero1[1], aero2[1], (n_elems * 2) + 2)
                            sweep_temp = np.linspace(aero1[2], aero2[2], (n_elems * 2) + 2)
                            airfoil_temp = np.linspace(aero1[3], aero2[3], (n_elems * num_node_elem + 2))
                            elastic_temp = np.linspace(aero1[4], aero2[4], (n_elems * 2) + 2)

                            chords[element_counter] = [chord_temp[0], chord_temp[1], chord_temp[2]]
                            twist[element_counter] = [twist_temp[0], twist_temp[1], twist_temp[2]]
                            sweep[element_counter] = [sweep_temp[0], sweep_temp[1], sweep_temp[2]]
                            elastic_axis[element_counter] = [elastic_temp[0], elastic_temp[1], elastic_temp[2]]

                            count_count = 3
                            for a in range(1, n_elems, 1):
                                three = 0
                                for aa in range(3):
                                    if aa == 2:
                                        chords[element_counter][three] = chords[element_counter - 1][2]
                                        twist[element_counter][three] = twist[element_counter - 1][2]
                                        sweep[element_counter][three] = sweep[element_counter - 1][2]
                                        elastic_axis[element_counter][three] = elastic_temp[element_counter - 1][2]
                                        continue
                                    else:
                                        chords[element_counter][three] = chord_temp[count_count]
                                        twist[element_counter][three] = twist_temp[count_count]
                                        sweep[element_counter][three] = sweep_temp[count_count]
                                        elastic_axis[element_counter][three] = elastic_temp[count_count]
                                    three += 1
                                    count_count += 1
                                    if k == 0 and j == 0 and i == 0 and aa == 0:
                                        break
                                    elif aa == 2:
                                        element_counter += 1

                            count_count = 0
                            for a in range(n_elems):
                                three = 0
                                for aa in range(3):
                                    if aero1[3] != aero2[3] and a == 0 and aa == 0:
                                        airfoil_distribution[element_counter][three] = airfoil_temp[count_count]
                                        pass
                                    elif aero1[3] != aero2[3] and a == n_elems - 1 and aa == 2:
                                        airfoil_distribution[element_counter][three] = airfoil_temp[count_count]
                                        pass
                                    elif aero1[3] != aero2[3]:
                                        airfoil_distribution[element_counter][three] = airfoil_hold
                                        added_airfoils.append(self.airfoil_inter(airfoil_temp[count_count], air_data))
                                        airfoil_hold += 1
                                    else:
                                        airfoil_distribution[element_counter][three] = airfoil_temp[count_count]
                                    three += 1
                                    count_count += 1
                                    if k == 0 and j == 0 and i == 0 and aa == 0:
                                        break
                                    elif aa == 2:
                                        element_counter += 1
                            continue

                        count_count = 0
                        for a in range(n_elems):
                            three = 0
                            for aa in range(3):
                                if aa == 0 and count_count != 0:
                                    chords[element_counter][three] = chords[element_counter - 1][1]
                                    twist[element_counter][three] = twist[element_counter - 1][1]
                                    sweep[element_counter][three] = sweep[element_counter - 1][1]
                                    elastic_axis[element_counter][three] = elastic_axis[element_counter - 1][1]
                                    airfoil_distribution[element_counter][three] = \
                                    airfoil_distribution[element_counter - 1][1]
                                    three += 1
                                    continue

                                elif aa == 0:
                                    chords[element_counter][three] = chord_temp[count_count]
                                    twist[element_counter][three] = twist_temp[count_count]
                                    sweep[element_counter][three] = sweep_temp[count_count]
                                    elastic_axis[element_counter][three] = elastic_temp[count_count]
                                    if aero1[3] != aero2[3] and a == 0:
                                        airfoil_distribution[element_counter][three] = airfoil_temp[count_count]
                                        pass

                                elif aa == 1:
                                    chords[element_counter][three] = chord_temp[count_count + 1]
                                    twist[element_counter][three] = twist_temp[count_count + 1]
                                    sweep[element_counter][three] = sweep_temp[count_count + 1]
                                    elastic_axis[element_counter][three] = elastic_temp[count_count + 1]
                                    if aero1[3] != aero2[3] and a == n_elems - 1:
                                        airfoil_distribution[element_counter][three] = airfoil_temp[-1]
                                        airfoil_hold += 1
                                        pass
                                    elif aero1[3] != aero2[3]:
                                        airfoil_distribution[element_counter][three] = airfoil_hold
                                        added_airfoils.append(self.airfoil_inter(airfoil_temp[count_count + 1], air_data))
                                        airfoil_hold += 1
                                    else:
                                        airfoil_distribution[element_counter][three] = airfoil_temp[count_count]

                                elif aa == 2:
                                    chords[element_counter][three] = chord_temp[count_count - 1]
                                    twist[element_counter][three] = twist_temp[count_count - 1]
                                    sweep[element_counter][three] = sweep_temp[count_count - 1]
                                    elastic_axis[element_counter][three] = elastic_temp[count_count - 1]
                                    if aero1[3] != aero2[3] and a == n_elems - 1:
                                        airfoil_distribution[element_counter][three] = airfoil_hold - 1
                                        added_airfoils.append(self.airfoil_inter(airfoil_temp[count_count - 1], air_data))
                                    elif aero1[3] != aero2[3]:
                                        airfoil_distribution[element_counter][three] = airfoil_hold
                                        added_airfoils.append(self.airfoil_inter(airfoil_temp[count_count - 1], air_data))
                                        airfoil_hold += 1
                                    else:
                                        airfoil_distribution[element_counter][three] = airfoil_temp[count_count]

                                three += 1
                                count_count += 1
                                if k == 0 and j == 0 and i == 0 and aa == 0:
                                    break
                                elif aa == 2:
                                    element_counter += 1

                        if node_counter == 0:
                            aero_node[node_counter] = True
                            for aaa in range(n_elems * 2):
                                node_counter += 1
                                aero_node[node_counter] = True
                        else:
                            aero_node[node_counter] = True
                            for aaa in range(n_elems * 2):
                                node_counter += 1
                                aero_node[node_counter] = True
                    except:
                        n_elems = data['origins'][k]['branch'][i]['section'][j]['n_elements']
                        for aaaa in range(n_elems):
                            if k == 0 and j == 0 and i == 0 and a == 0:
                                pass
                            else:
                                element_counter += 1
                        for aaa in range(n_elems * 2):
                            node_counter += 1
                            aero_node[node_counter] = False

        # surface_distribution
        surface_distribution = []

        num_surfaces = 0

        # surface_m chordwise panelling
        surface_m = []

        counter_surf = 0
        for k in range(len(data['origins'])):
            for i in range(len(data['origins'][k]['branch'])):
                try:
                    chordwise = data['origins'][k]['branch'][i]['chordwise']
                    surface_m.append(chordwise)
                    # counter_surf += 1
                except:
                    pass
                for j in range(len(data['origins'][k]['branch'][i]['section'])):
                    n_elems = data['origins'][k]['branch'][i]['section'][j]["n_elements"]
                    for a in range(n_elems):
                        try:
                            if isinstance(data['origins'][k]['branch'][i]['section'][j]["aero"][0], float) == True:
                                surface_distribution.append(counter_surf)
                        except:
                            surface_distribution.append(-1)

                    try:
                        if surface_distribution[-1] == -1 and surface_distribution[-2] >= 0:
                            counter_surf += 1
                            continue
                    except:
                        pass

                    if j == len(data['origins'][k]['branch'][i]['section']) - 1 and surface_distribution[-1] != -1:
                        counter_surf += 1

        # m_distribution, discretisation method

        m_distribution = "uniform"

# -----------------------------------------------------------------------------------------------------
#----------------------CONTROL SURFACES----------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------

        # control surfaces
        control_surface = np.ones((num_elem, num_node_elem), dtype=int)
        control_surface = control_surface * -1

        # control surface type
        control_surface_type = []

        # control surface chord
        control_surface_chord = []

        # control surface hinge coord
        control_surface_hinge_coord = []

        # control surface deflection
        control_surface_deflection = []

        control_start = []
        control_end = []

        # airfoil efficiency
        # airfoil_efficiency = np.zeros((num_elem, num_node_elem, 2, 3))

        element_counter = 0
        control_count = 0
        surface_interp = []
        for k in range(len(data['origins'])):
            for i in range(len(data['origins'][k]['branch'])):
                for j in range(len(data['origins'][k]['branch'][i]['section'])):
                    try:
                        control = data['origins'][k]['branch'][i]['section'][j]['control']
                        control2 = data['origins'][k]['branch'][i]['section'][j]['control2']
                        coord = data['origins'][k]['branch'][i]['section'][j]['coord']
                        n_elems = data['origins'][k]['branch'][i]['section'][j]['n_elements']
                        surface_type_temp = np.linspace(int(control[0]), int(control2[0]), n_elems * num_node_elem)
                        surface_chord_temp = np.linspace(int(control[1]), int(control2[1]), n_elems * num_node_elem)
                        surface_deflection_temp = np.linspace(control[2], control2[2], n_elems * num_node_elem)

                        surface_interp.append(1)

                        control_start.append(coord[1])
                        control_end.append(coord[2])

                        for a in range(len(surface_deflection_temp)):
                            control_surface_type.append(int(surface_type_temp[a]))
                            control_surface_chord.append(int(surface_chord_temp[a]))
                            control_surface_deflection.append(surface_deflection_temp[a])
                            control_surface_hinge_coord.append(coord[0])

                    except:
                        try:
                            control = data['origins'][k]['branch'][i]['section'][j]['control']
                            coord = data['origins'][k]['branch'][i]['section'][j]['coord']
                            n_elems = data['origins'][k]['branch'][i]['section'][j]['n_elements']
                            surface_interp.append(0)

                            control_surface_type.append((control[0]))
                            control_surface_chord.append(int(control[1]))
                            control_surface_deflection.append(control[2])
                            control_surface_hinge_coord.append(coord[0])
                            control_start.append(coord[1])
                            control_end.append(coord[2])

                        except:
                            pass

        # Finds the nodes closest to the control surfaces
        x_distance = np.zeros(len(x))
        y_distance = np.zeros(len(x))
        z_distance = np.zeros(len(x))
        distance = np.zeros(len(x))

        position = np.column_stack((x, y, z))

        control_nodes_start = []
        control_nodes_end = []
        for i in range(len(control_start)):
            for j in range(len(x)):
                x_distance[j] = control_start[i][0] - x[j]
                y_distance[j] = control_start[i][1] - y[j]
                z_distance[j] = control_start[i][2] - z[j]
                distance[j] = np.sqrt(x_distance[j] ** 2 + y_distance[j] ** 2 + z_distance[j] ** 2)
            minimal = min(distance)
            min_index = [i for i, j in enumerate(distance) if j == minimal]  # finds the minimum index to make th
            min_index = str(min_index)[1:-1]
            min_index = min_index.split(", ")
            min_index = int(min_index[0])
            control_nodes_start.append(min_index)

        for i in range(len(control_end)):
            for j in range(len(x)):
                x_distance[j] = control_end[i][0] - x[j]
                y_distance[j] = control_end[i][1] - y[j]
                z_distance[j] = control_end[i][2] - z[j]
                distance[j] = np.sqrt(x_distance[j] ** 2 + y_distance[j] ** 2 + z_distance[j] ** 2)
            minimal = min(distance)
            min_index = [i for i, j in enumerate(distance) if j == minimal]  # finds the minimum index to make th
            min_index = str(min_index)[1:-1]
            min_index = min_index.split(", ")
            min_index = int(min_index[0])
            control_nodes_end.append(min_index)
            if control_nodes_end[i] < control_nodes_start[i]:
                control_nodes_end[i] = control_nodes_start[i]

                # finding duplicate points
                conn_one = list(conn[:, 0])
                dupps = [x for x in conn_one if conn_one.count(x) > 1]
                duppss = list(set(dupps))
                for aa in range(0, len(duppss) - 1):
                    if duppss[aa] > duppss[aa + 1]:
                        small = duppss[aa + 1]
                        large = duppss[aa]
                        duppss[aa] = small
                        duppss[aa + 1] = large

                for aa in duppss:
                    if aa < control_nodes_start[i]:
                        hold = aa

                control_nodes_start[i] = aa

        surface_counter = 0
        interp_count = 0
        hold = -1
        for i in range(num_elem):
            try:
                if conn[i][0] == control_nodes_start[surface_counter]:
                    hold = interp_count
            except:
                pass

            control_surface[i][0] = hold
            control_surface[i][1] = hold
            control_surface[i][2] = hold

            if hold != -1 and surface_interp[surface_counter] == 1:
                control_surface[i][0] = hold
                control_surface[i][1] = hold + 2
                control_surface[i][2] = hold + 1
                hold = hold + 3
                interp_count = hold

            try:
                if conn[i][2] == control_nodes_end[surface_counter] - 1:
                    hold = -1
                    if surface_interp[surface_counter] == 0:
                        interp_count += 1

                    surface_counter += 1
            except:
                pass

        for i in range(len(added_airfoils)):
            air_data.append(added_airfoils[i])

        control_surface_deflection = np.array(control_surface_deflection)

# ------------------------------------------------------------------------------------------------------
# ----------------------AERO WRITE----------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------

        with h5.File(route + '/' + case_name + '.aero.h5', 'a') as h5file:
            airfoils_group = h5file.create_group('airfoils')

            for i in range(len(air_data)):  # len(air_data)
                airfoil = airfoils_group.create_dataset(str(i), data=air_data[i])
            # chord
            chord_input = h5file.create_dataset('chord', data=chords)
            dim_attr = chord_input.attrs['units'] = 'm'

            # twist
            twist_input = h5file.create_dataset('twist', data=twist)
            dim_attr = twist_input.attrs['units'] = 'rad'

            # sweep
            sweep_input = h5file.create_dataset('sweep', data=sweep)
            dim_attr = sweep_input.attrs['units'] = 'rad'

            # airfoil distribution
            airfoil_distribution_input = h5file.create_dataset('airfoil_distribution', data=airfoil_distribution)

            surface_distribution_input = h5file.create_dataset('surface_distribution', data=surface_distribution)
            surface_m_input = h5file.create_dataset('surface_m', data=surface_m)
            m_distribution_input = h5file.create_dataset('m_distribution',
                                                         data=m_distribution.encode('ascii', 'ignore'))

            aero_node_input = h5file.create_dataset('aero_node', data=aero_node)
            elastic_axis_input = h5file.create_dataset('elastic_axis', data=elastic_axis)

            control_surface_input = h5file.create_dataset('control_surface', data=control_surface)
            control_surface_deflection_input = h5file.create_dataset('control_surface_deflection',
                                                                     data=control_surface_deflection)
            control_surface_chord_input = h5file.create_dataset('control_surface_chord', data=control_surface_chord)
            control_surface_hinge_coord_input = h5file.create_dataset('control_surface_hinge_coord',
                                                                      data=control_surface_hinge_coord)
            control_surface_types_input = h5file.create_dataset('control_surface_type', data=control_surface_type)
    
    def excel(self, data):
        
        
        def linspace(start, end, n, flag):
            """function used to linearly space the various paramaters

            Args:
                start (float): starting value
                end (float): ending value
                n_nodes (float): number of nodes of the given section
                flag (0 or 1): conditional statement for starting section and subsequent sections 

            Returns:
                list[float]: list of neccesary values at each node or element 
            """
            if flag == 0:
                val = np.linspace(start, end, n)
            elif flag == 1:
                start = (end - start) / n + start
                val = np.linspace(start, end, n)
            return val
        
        df = pd.read_excel(data, sheet_name = 'Master')
        rows = df.shape[0]

        empt = []
        chords = np.array(empt) 
        twist =  np.array(empt) 
        sweep = np.array(empt)  
        airfoil_distribution = np.array(empt)   
        elastic_axis = np.array(empt)  
        surface_distribution  = np.array(empt, dtype = int)    
        aero_node  = np.array(empt, dtype= bool)    

        num_node = 0
        num_elem = 0
        node_hold = 1
        elem_hold = 0
        surface_num = 0
        m_panel = []
        start = 1
        aero_start = 1
        sym_before = 0

        for i in range(1, rows):
            chord_1 = df.values[i, 36]
            twist_1 = df.values[i, 37]
            sweep_1 = df.values[i, 38]
            airfoil_1 = df.values[i, 39]
            ea_loc = df.values[i, 40]

            chord_2 = df.values[i, 41]
            twist_2 = df.values[i, 42]
            sweep_2 = df.values[i, 43]
            airfoil_2 = df.values[i, 44]
            ea_loc_2 = df.values[i, 45]

            x_2 = df.values[i, 5]
            y_2 = df.values[i, 6]
            z_2 = df.values[i, 7]
                
            
            if str(df.values[i, 34]).lower() == "yes"  and start == 1:
                start = 0
                if i == 1:
                    node_hold = 1
                    elem_hold = num_elem
                
                else:
                    node_hold = num_node
                    elem_hold = num_elem 
                
                if i != 1:
                    sym_before = 1
                
            n_elem = df.values[i, 1]
            if str(n_elem) == 'nan':
                break  
            num_elem += n_elem 

            if i == 1:
                n_nodes = df.values[i, 1] * 2 + 1
                num_node += n_nodes
                flag = 0
                
            else:
                n_nodes = df.values[i, 1] * 2 
                num_node += n_nodes
                flag = 1
            
            # indicates that there is a new aero section
            if i != 1 and  df.values[i-1, 36] <= 0 and chord_1 > 0:
                aero_start = 1
            
            # adds another surface if there is a discontinuity in an aero surface      
            if [df.values[i, 2], df.values[i, 3], df.values[i, 4]] != [df.values[i-1, 5], df.values[i-1, 6], df.values[i-1, 7]]:
                aero_start = 1
                if chord_1 > 0 and i != 1 and aero_start == 1:
                    surface_num += 1

            # continuous section, but it is an aero surface and the previous section was not an aero surface
            elif chord_1 > 0 and surface_distribution[-1] == -1:
                surface_num += 1
                
            # adds surface if the preious section was a vertical fin leading to h-tail
            elif i!= 1 and df.values[i-1, 3] == df.values[i-1, 6]:
                if y_2 != df.values[i-1, 3]:
                    surface_num += 1
                    aero_start = 1
                
            # signifies that it is an aero section
            if chord_1 > 0:
                      
                if aero_start == 1:
                    m_panel.append(df.values[i, 46])
                    aero_start = 0
                    
                for ii in range(n_elem):
                    surface_distribution = np.append(surface_distribution, surface_num)
                
                # catches the ending of the previous section which is now an aero node
                if surface_num !=0:
                    aero_node[-1] = True

                for ii in range(n_nodes):
                    aero_node = np.append(aero_node, True)
            
            # identifies that a section is not an aero section  
            if chord_1 <= 0:
                for ii in range(n_elem):
                    surface_distribution = np.append(surface_distribution, -1)

                for ii in range(n_nodes):
                    aero_node = np.append(aero_node, False)
            
            chords = np.append(chords, linspace(chord_1, chord_2, n_nodes, flag))      
            twist = np.append(twist, linspace(twist_1, twist_2, n_nodes, flag))      
            sweep = np.append(sweep, linspace(sweep_1, sweep_2, n_nodes, flag))  
            airfoil_distribution = np.append(airfoil_distribution, linspace(airfoil_1, airfoil_2, n_nodes, flag))  
            elastic_axis = np.append(elastic_axis, linspace(ea_loc, ea_loc_2, n_nodes, flag))  

            # symmetric data values added here
            # checks that the ending pionts are not the same as the starting points of the next section. This signifies and change in section
            if i != rows - 1 and [x_2, y_2, z_2] != [df.values[i + 1 , 2], df.values[i + 1, 3], df.values[i + 1, 4]]:
                start = 1
                # checks to see that the symmetric flag is on. If it is, the other wing section paramaters are added 
                if str(df.values[i, 34]).lower() == "yes":
                    #surface_num += 1
                    for ii in range(node_hold, len(chords)):
                        chords = np.append(chords, chords[ii])
                        twist = np.append(twist, twist[ii])
                        sweep = np.append(sweep, sweep[ii])
                        airfoil_distribution = np.append(airfoil_distribution, airfoil_distribution[ii])  
                        elastic_axis = np.append(elastic_axis, elastic_axis[ii])  
                        aero_node = np.append(aero_node, aero_node[ii])
                        num_node += 1
                        
                    for ii in range(elem_hold, len(surface_distribution)):
                        # updating surface # for other side 
                        value = surface_distribution[ii]
                        if value != -1:
                            value += 1
                        surface_distribution = np.append(surface_distribution, value)
                        num_elem += 1
                    
                    if chord_1 > 0:
                        m_panel.append(m_panel[-1])
                        
                    surface_num = max(surface_distribution)
                        
            # if the last line is symmetric, then it adds the other side                         
            elif i == rows -1 and str(df.values[i, 34]).lower() == "yes":
                start = 1
                for ii in range(node_hold, len(chords)):
                    chords = np.append(chords, chords[ii])
                    twist = np.append(twist, twist[ii])
                    sweep = np.append(sweep, sweep[ii])
                    airfoil_distribution = np.append(airfoil_distribution, airfoil_distribution[ii])  
                    elastic_axis = np.append(elastic_axis, elastic_axis[ii])  
                    aero_node = np.append(aero_node, True)
                    num_node += 1
        
                for ii in range(elem_hold, len(surface_distribution)):
                    # updating surface # for other side 
                    value = surface_distribution[ii]
                    if value != -1:
                        value += 1
                    surface_distribution = np.append(surface_distribution, value)
                    num_elem += 1

                if chord_1 > 0:
                    m_panel.append(m_panel[-1])
                        
                surface_num = max(surface_distribution) + 1
                
        # finalizing matrices
        chord = np.zeros((self.num_elem, 3))
        twists = np.zeros((self.num_elem, 3))
        sweeps = np.zeros((self.num_elem, 3))
        elastic_axes = np.zeros((self.num_elem, 3))
        airfoil_temp = np.zeros((self.num_elem, 3))

        count = 2
        for i in range(len(self.conn)):
            con = self.conn 
            chord[i, ...] = chords[int(con[i, 0])], chords[int(con[i, 1])], chords[int(con[i, 2])]
            twists[i, ...] = twist[int(con[i, 0])], twist[int(con[i, 1])], twist[int(con[i, 2])]
            sweeps[i, ...] = sweep[int(con[i, 0])], sweep[int(con[i, 1])], sweep[int(con[i, 2])]
            elastic_axes[i, ...] = elastic_axis[int(con[i, 0])], elastic_axis[int(con[i, 1])], elastic_axis[int(con[i, 2])]
            airfoil_temp[i, ...] = airfoil_distribution[int(con[i, 0])], airfoil_distribution[int(con[i, 1])], airfoil_distribution[int(con[i, 2])]
            
            # fixes nodes that are shared by multiple sections
            if aero_node[int(con[i,1])] == 0:
                first = 0
                chord[i, 0] = first
                twists[i, 0] = first
                sweeps[i, 0] = first
                elastic_axes[i, 0] = first
                airfoil_temp[i, 0] = first

            if aero_node[int(con[i,1])] == 1:
                chord[i, 0] = (chord[i, 1] - chord[i, 2]) * -1 + chord[i, 2]
                twists[i, 0] = (twists[i, 1] - twists[i, 2]) * -1 + twists[i, 2]
                sweeps[i, 0] = (sweeps[i, 1] - sweeps[i, 2]) * -1 + sweeps[i, 2]
                elastic_axes[i, 0] = (elastic_axes[i, 1] - elastic_axes[i, 2]) * -1 + elastic_axes[i, 2]
                airfoil_temp[i, 0] = (airfoil_temp[i, 1] - airfoil_temp[i, 2]) * -1 + airfoil_temp[i, 2]
                
#----------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------POLARS---------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
        dfP = pd.read_excel(data, sheet_name = 'Polars')
        rows = dfP.shape[0]
        column = dfP.shape[1]
        polar = []
        
        col_start = 0
        for i in range(int(column/4)):
            pol_temp = []
            if i != 0:
                col_start += 4
            for j in range(2, rows):
                check = [str(dfP.values[j, col_start]), str(dfP.values[j, col_start+1]), str(dfP.values[j, col_start+2]), str(dfP.values[j, col_start+3])]
                if 'nan' in check:
                    pol_temp = np.array(pol_temp)
                    polar.append(pol_temp)
                    break 
                pol_temp.append([dfP.values[j, col_start], dfP.values[j, col_start+1], dfP.values[j, col_start+2], dfP.values[j, col_start+3]])
                if j == rows-1:
                    pol_temp = np.array(pol_temp)
                    polar.append(pol_temp) 
        
        # Calculating new airfoils from iterpolation
        airfoil_dis_final = np.zeros((self.num_elem, 3), dtype = int)
        for i in range(self.num_elem):
            airfoil_dis_final[i, ...] = airfoil_temp[i, 0], airfoil_temp[i, 1], airfoil_temp[i, 2]
                 
        new_airfoils = []
        interp_polars = []
        value_history = []
        count_history = []
        interp_write_flag = False
        count = int(len(self.airfoil_data))
        
        # INPUTS to drag polar interpolation
        cl_in = []
        cd_in = []
        cm_in = []
        for z in range(len(polar)):
            in_1 = []
            in_2 = []
            in_3 = []
            for zz in range(len(polar[z])):
                in_1.append([polar[z][zz][0], polar[z][zz][1]])
                in_2.append([polar[z][zz][0], polar[z][zz][2]])
                in_3.append([polar[z][zz][0], polar[z][zz][3]])
                if zz == len(polar[z]) - 1:
                    cl_in.append(in_1)
                    cd_in.append(in_2)
                    cm_in.append(in_3)
        
        for i in range(self.num_elem):
            for j in range(1,3):
                if airfoil_temp[i, j] % 1 != 0:
                    if i!= 0 and j == 2:
                        airfoil_dis_final[i, 0] = airfoil_dis_final[i-1, 1]

                # calculating new airfoil sections from interpolating and preventing duplicate airfoil production
                    if len(value_history) == 0:
                        new_airfoils.append(self.airfoil_inter(airfoil_temp[i, j], self.airfoil_data))
                        
                        # INTERPOLATION OF POLARS
                        new_cl = self.airfoil_inter(airfoil_temp[i,j], cl_in)
                        new_cl = np.array(new_cl)
                        new_cd = self.airfoil_inter(airfoil_temp[i,j], cd_in)
                        new_cd = np.array(new_cd)
                        new_cm = self.airfoil_inter(airfoil_temp[i,j], cm_in)
                        new_cm = np.array(new_cm)
                        interp_polars.append(np.column_stack((new_cl[:, 0], new_cl[:, 1], new_cd[:, 1], new_cm[:, 1])))
                            
                        value_history.append(airfoil_temp[i,j])
                        count_history.append(count)
                        airfoil_dis_final[i, j] = count
                        count += int(1)

                    else:
                        for ii in range(len(value_history)):
                            if value_history[ii] == airfoil_temp[i,j]:
                                airfoil_dis_final[i, j] = count_history[ii]
                                break
                            
                            elif ii == len(value_history) - 1:
                                new_airfoils.append(self.airfoil_inter(airfoil_temp[i, j], self.airfoil_data))
                                
                                # INTERPOLATION OF POLARS
                                new_cl = self.airfoil_inter(airfoil_temp[i,j], cl_in)
                                new_cl = np.array(new_cl)
                                new_cd = self.airfoil_inter(airfoil_temp[i,j], cd_in)
                                new_cd = np.array(new_cd)
                                new_cm = self.airfoil_inter(airfoil_temp[i,j], cm_in)
                                new_cm = np.array(new_cm)
                                interp_polars.append(np.column_stack((new_cl[:, 0], new_cl[:, 1], new_cd[:, 1], new_cm[:, 1])))
                                
                                value_history.append(airfoil_temp[i,j]) 
                                count_history.append(count)
                                airfoil_dis_final[i, j] = count
                                count += int(1)
                                break
                            
        # appending interpolated airfoils
        for i in range(len(new_airfoils)):
            self.airfoil_data.append(new_airfoils[i])
            polar.append(interp_polars[i])
        

        surface_m = []
        for i in range(len(m_panel)):
            if m_panel[i] > 0:
                surface_m.append(m_panel[i])
                
        m_distribution = 'uniform'
              
# -----------------------------------------------------------------------------------------------------
#----------------------CONTROL SURFACES----------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------
        # control surfaces
        control_read = pd.read_excel(data, sheet_name = 'Control_Surfaces')
        rows = control_read.shape[0]

        # collecting inputs from input file
        ref_num = []
        surface_type_temp = []
        control_chord_temp = []
        deflect_start = []
        deflect_end = []
        control_x_start = []
        control_y_start = []
        control_z_start= []
        control_x_end= []
        control_y_end= []
        control_z_end= []
        surface_interp = []
        sym = []
        
        # setting default tail_index value
        self.tail_index = 0 

        for i in range(1,rows):
            ref_num.append(control_read.values[i, 0])
            surface_type_temp.append(control_read.values[i, 2])
            control_chord_temp.append(control_read.values[i, 3])
            deflect_start.append(control_read.values[i, 4])
            control_x_start.append(control_read.values[i, 5])
            control_y_start.append(control_read.values[i, 6])
            control_z_start.append(control_read.values[i, 7])
            deflect_end.append(control_read.values[i, 8])
            control_x_end.append(control_read.values[i, 9])
            control_y_end.append(control_read.values[i, 10])
            control_z_end.append(control_read.values[i, 11])
            sym.append(0)
                                  
            # flag to mark which control surface will be used in static trim           
            if str(control_read.values[i, 13]).lower() == 'yes':
                self.tail_index = ref_num[-1]
                
            # flaging if the section has interpolation 
            if control_read.values[i,4] != control_read.values[i,8]:
                surface_interp.append(1)
            else: 
                surface_interp.append(0)
                
            if str(control_read.values[i,12]).lower() == 'yes':
                control_x_start.append(control_read.values[i, 5])
                control_y_start.append(-control_read.values[i, 6])
                control_z_start.append(control_read.values[i, 7])
                control_x_end.append(control_read.values[i, 9])
                control_y_end.append(-control_read.values[i, 10])
                control_z_end.append(control_read.values[i, 11])
                surface_type_temp.append(control_read.values[i, 2])
                control_chord_temp.append(control_read.values[i, 3])
                deflect_start.append(control_read.values[i, 4])
                deflect_end.append(control_read.values[i, 8])
                sym[-1] = 1
                sym.append(1)
                ref_num.append(ref_num[-1])
                surface_interp.append(surface_interp[-1])

            
        # finding the starting and ending nodes 
        pos = np.column_stack((self.x, self.y, self.z))        
        distance = np.zeros(len(self.x))
        control_nodes_start = []
        control_nodes_end = []
        con = self.conn

        # finding the starting nodes of the control surfaces
        for i in range(len(control_x_start)):
            for j in range(len(self.x)):
                x_distance = control_x_start[i] - self.x[j]
                y_distance = control_y_start[i] - self.y[j]
                z_distance = control_z_start[i] - self.z[j]
                distance[j] = np.sqrt(x_distance ** 2 + y_distance ** 2 + z_distance ** 2)
            minimal = min(distance)
            min_index = [i for i, j in enumerate(distance) if j == minimal]  # finds the minimum index
            min_index = str(min_index)[1:-1]
            min_index = min_index.split(", ")
            min_index = int(min_index[0])
            control_nodes_start.append(min_index)

        # finding the ending nodes of the control surfaces
        for i in range(len(control_x_start)):
            for j in range(len(self.x)):
                x_distance = control_x_end[i] - self.x[j]
                y_distance = control_y_end[i] - self.y[j]
                z_distance = control_z_end[i] - self.z[j]
                distance[j] = np.sqrt(x_distance ** 2 + y_distance ** 2 + z_distance ** 2)
            minimal = min(distance)
            min_index = [i for i, j in enumerate(distance) if j == minimal]  # finds the minimum index
            min_index = str(min_index)[1:-1]
            min_index = min_index.split(", ")
            min_index = int(min_index[0])
            control_nodes_end.append(min_index)
            if control_nodes_end[i] < control_nodes_start[i]:
                control_nodes_end[i] = control_nodes_start[i]

                # finding duplicate points
                conn_one = list(con[:, 0])
                dupps = [x for x in conn_one if conn_one.count(x) > 1]
                duppss = list(set(dupps))
                for aa in range(0, len(duppss) - 1):
                    if duppss[aa] > duppss[aa + 1]:
                        small = duppss[aa + 1]
                        large = duppss[aa]
                        duppss[aa] = small
                        duppss[aa + 1] = large

                for aa in duppss:
                    if aa < control_nodes_start[i]:
                        hold = aa

                control_nodes_start[i] = aa    
        
        control_surface = np.zeros((num_elem, 3), dtype=int) - 1
        con = self.conn
        
        surface_counter = 0
        hold = -1
        stop = 0
        start_con = []
        end_con = []

        # Combine the original list and the reference numbers into tuples
        combined_list = [(value, ref_num[i]) for i, value in enumerate(control_nodes_start)]

        # making sure that the reference number is changed in the same way that the start locations change
        sorted_list = sorted(combined_list, key=lambda x: x[0])

        # Extract the sorted values and their associated reference numbers into separate lists
        control_nodes_start = [t[0] for t in sorted_list]
        ref_num = [t[1] for t in sorted_list]
        
        control_nodes_end.sort()
        # If there are control surfaces
        if control_nodes_start != []:
            for i in range(num_elem):
                if stop == 1:
                    break
                for j in [0, 2, 1]:  
                    # starting control point cannot start on end point of con (column 1)
                    if con[i][j] == control_nodes_start[surface_counter] and j!= 1:
                        hold = ref_num[surface_counter]
                        start_con.append([i,j]) 
                        
                    control_surface[i][j] = hold
                        
                    # ending control point cannot start on start0 point of con (column 0)
                    if con[i][j] == control_nodes_end[surface_counter] and j!= 0:
                        surface_counter += 1
                        end_con.append([i,j]) 
                        hold = -1
                    
                        if surface_counter > len(control_nodes_end)-1 or surface_counter > len(control_nodes_start) - 1:
                            stop = 1
                            break
        
        
            control_surface_type = np.zeros(max(ref_num) + 1, dtype= int) # needs to be dtype int 
            control_surface_deflection = np.zeros(max(ref_num) + 1)  
            control_surface_chord = np.zeros(max(ref_num) + 1, dtype= int)  # needs to be dtype int 

            # getting the m_panel for each control surface 
            m_panel_save = np.zeros(max(ref_num) + 1) 
            for i in range(len(control_nodes_end)):
                go = False 
                for j in range(len(con)):
                    for k in [0, 2, 1]:
                        if con[j, k] == control_nodes_end[i]:
                            m_panel_save[ref_num[i]] = m_panel[surface_distribution[j]]
                            go = True
                            break
                    if go == True:
                        break
                    
            # adding control surface variables based off ref_num 
            for i in range(len(ref_num)):
                if i == 0:
                    control_surface_type[ref_num[i]] = int(surface_type_temp[i])
                    control_surface_deflection[ref_num[i]] = deflect_start[i]
                    control_surface_chord[ref_num[i]] = int(control_chord_temp[i]*m_panel_save[ref_num[i]])
                
                if i != 0 and ref_num[i] != ref_num[i-1]:
                    control_surface_type[ref_num[i]] = int(surface_type_temp[i])
                    control_surface_deflection[ref_num[i]] = deflect_start[i]
                    control_surface_chord[ref_num[i]] = int(control_chord_temp[i]*m_panel_save[ref_num[i]])
            
            # Adding interpolated control surfaces
            interp_count = len(control_surface_deflection) + 1
            lin_space_count = 0
            go = False
            start_con = np.array(start_con)
            end_con = np.array(end_con)
            for i in range(len(surface_interp)):
                if surface_interp[i] == 0:
                    continue
                else:
                    for j in range(start_con[i, 0], end_con[i, 0] + 1):
                        for k in [0, 2, 1]:
                            # if it is the first row, checks to see where to start applying interpt sections
                            if j == start_con[i, 0] and k == start_con[i, 1]:
                                go = True
                                    
                            if j == end_con[i, 0] and k == end_con[i, 1]:
                                control_surface[j, k] = interp_count
                                interp_count += 1
                                lin_space_count += 1
                                # appending interp control surface data to appopriate variable
                                deflect_nu = np.linspace(deflect_start[i], deflect_end[i], lin_space_count + 1)
                                type_nu = np.linspace(surface_type_temp[i], surface_type_temp[i], lin_space_count + 1, dtype= int)
                                chord_nu = np.linspace(control_chord_temp[i] * m_panel[surface_distribution[j]] , control_chord_temp[i] * m_panel[surface_distribution[j]], lin_space_count + 1, dtype = int)
                                
                                for kk in range(len(deflect_nu)):
                                    control_surface_type = np.append(control_surface_type, type_nu[kk])
                                    control_surface_deflection = np.append(control_surface_deflection, deflect_nu[kk])
                                    control_surface_chord = np.append(control_surface_chord, chord_nu[kk])
                                
                                lin_space_count = 0
                                go = False
                            
                            if go == True and k != start_con[i, 1]:
                                control_surface[j, k] = interp_count
                                interp_count += 1 
                                lin_space_count += 1
                            
            control_surface_hinge_coord = np.zeros((len(control_surface_type),), dtype=float)
            for i in range(1,rows):
                index = int(control_read.values[i, 0])
                control_surface_hinge_coord[index] = -control_read.values[i, 1]
        else:
            control_surface_type = np.zeros(1, dtype= int) # needs to be dtype int 
            control_surface_deflection = np.zeros(1)  
            control_surface_chord = np.zeros(1, dtype= int)  # needs to be dtype int 
            control_surface_hinge_coord = np.zeros((len(control_surface_type),), dtype=float)
            
        with h5.File(self.route + '/' + self.case + '.aero.h5', 'a') as h5file:
            airfoils_group = h5file.create_group('airfoils')
            polars_group = h5file.create_group('polars')

            for i in range(len(self.airfoil_data)):  # len(air_data)
                airfoil = airfoils_group.create_dataset(str(i), data=self.airfoil_data[i])
                  # New Polar Variable
                polars_input = polars_group.create_dataset(str(i), data=polar[i])
            
            # chord
            chord_input = h5file.create_dataset('chord', data=chord)
            dim_attr = chord_input.attrs['units'] = 'm'

            # twist
            twist_input = h5file.create_dataset('twist', data=twists)
            dim_attr = twist_input.attrs['units'] = 'rad'

            # sweep
            sweep_input = h5file.create_dataset('sweep', data=sweeps)
            dim_attr = sweep_input.attrs['units'] = 'rad'

            # airfoil distribution
            airfoil_distribution_input = h5file.create_dataset('airfoil_distribution', data=airfoil_dis_final)
            surface_distribution_input = h5file.create_dataset('surface_distribution', data=surface_distribution)
            surface_m_input = h5file.create_dataset('surface_m', data=surface_m)
            m_distribution_input = h5file.create_dataset('m_distribution', data=m_distribution.encode('ascii', 'ignore'))
            aero_node_input = h5file.create_dataset('aero_node', data=aero_node)
            elastic_axis_input = h5file.create_dataset('elastic_axis', data=elastic_axes)
            control_surface_input = h5file.create_dataset('control_surface', data=control_surface)
            control_surface_deflection_input = h5file.create_dataset('control_surface_deflection', data=control_surface_deflection)
            control_surface_chord_input = h5file.create_dataset('control_surface_chord', data=control_surface_chord)
            control_surface_hinge_coord_input = h5file.create_dataset('control_surface_hinge_coord', data=control_surface_hinge_coord)
            control_surface_types_input = h5file.create_dataset('control_surface_type', data=control_surface_type)
          
            
        