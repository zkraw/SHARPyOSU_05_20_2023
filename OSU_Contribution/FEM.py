from bdb import Breakpoint
from re import L
import h5py as h5
import numpy as np
import pandas as pd
from sharpy.solvers.beam_opt import Beam_Opt
import sharpy.utils.model_utils as matrix
import matplotlib.pyplot as plt


class FEM():
    def __init__(self, data, route, case_name, excel):
        self.data = data
        self.route = route
        self.case_name = case_name
        self.excel = excel
        self.beam_opt_nodes = []
        self.beam_opt_elems = []
        self.beam_opt_sym = []
        
        if self.excel == True:
           self.Excel_Calc()
            
        if self.excel == False:
            self.FEM_calc(self.data)
            
    
    def FEM_calc(self, data):
        num_node_elem = data['num_node_elem']
        # total number of elements
        num_elem = 0
# -----------------------------------------------------------------------------------------------------
#----------------------COORDINATES---------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------
        elems = []
        origins_x = []
        origins_y = []
        origins_z = []
        branch_end = []
        origin_end = []
        origin_elem = []
        piv_temp = []
        piv_log = []
        for k in range(len(data['origins'])):
            origins_x.append(data['origins'][k]['position'][0])
            origins_y.append(data['origins'][k]['position'][1])
            origins_z.append(data['origins'][k]['position'][2])
            origin_end.append(sum(elems))
            origin_elem.append(2 * sum(elems))
            for i in range(len(data['origins'][k]['branch'])):
                for j in range(len(data['origins'][k]['branch'][i]['section'])):
                    num_elem += data['origins'][k]['branch'][i]['section'][j]['n_elements']
                    elems.append(data['origins'][k]['branch'][i]['section'][j]['n_elements'])
                    piv_temp.append(data['origins'][k]['branch'][i]['section'][j]['pivot'])
                    if j == len(data['origins'][k]['branch'][i]['section']) - 1:
                        branch_end.append(sum(elems))

        for i in range(len(elems)):
            for j in range(elems[i]):
                if piv_temp[i] != "continue":
                    piv_log.append(piv_temp[i])
                elif piv_temp[i] == "continue":
                    piv_log.append("false")

        # total number of nodes
        num_node = 1
        for i in range(len(elems)):
            num_node += 2 * elems[i]

        # coordinates of nodes
        hold_x = 0
        hold_y = 0
        hold_z = 0

        x_start = []
        y_start = []
        z_start = []

        x_end = []
        y_end = []
        z_end = []

        # making start and end points
        for k in range(len(data['origins'])):
            for i in range(len(data['origins'][k]['branch'])):
                hold_x = origins_x[k]
                hold_y = origins_y[k]
                hold_z = origins_z[k]

                x_start.append(hold_x); y_start.append(hold_y); z_start.append(hold_z)

                for j in range(len(data['origins'][k]['branch'][i]['section'])):
                    x_end.append(data['origins'][k]['branch'][i]['section'][j]['position'][0])
                    y_end.append(data['origins'][k]['branch'][i]['section'][j]['position'][1])
                    z_end.append(data['origins'][k]['branch'][i]['section'][j]['position'][2])

                    if data['origins'][k]['branch'][i]['section'][j]['pivot'] == "true":
                        hold_x = data['origins'][k]['branch'][i]['section'][j]['position'][0]
                        hold_y = data['origins'][k]['branch'][i]['section'][j]['position'][1]
                        hold_z = data['origins'][k]['branch'][i]['section'][j]['position'][2]
                        if j != len(data['origins'][k]['branch'][i]['section']) - 1:
                            x_start.append(hold_x); y_start.append(hold_y); z_start.append(hold_z)

                    elif data['origins'][k]['branch'][i]['section'][j]['pivot'] == "continue":
                        if j != len(data['origins'][k]['branch'][i]['section']) - 1:
                            x_start.append(x_end[-1]);
                            y_start.append(y_end[-1]);
                            z_start.append(z_end[-1])

                    elif data['origins'][k]['branch'][i]['section'][j]['pivot'] == "false":
                        if j != len(data['origins'][k]['branch'][i]['section']) - 1:
                            x_start.append(hold_x);
                            y_start.append(hold_y);
                            z_start.append(hold_z)

        x = []
        y = []
        z = []

        # setting the coordinates
        x.append(x_start[0])
        y.append(y_start[0])
        z.append(z_start[0])

        for i in range(len(x_start)):
            x_temp = np.linspace(x_start[i], x_end[i], 2 * elems[i] + 1)
            y_temp = np.linspace(y_start[i], y_end[i], 2 * elems[i] + 1)
            z_temp = np.linspace(z_start[i], z_end[i], 2 * elems[i] + 1)
            for j in range(1, len(x_temp)):
                x.append(x_temp[j])
                y.append(y_temp[j])
                z.append(z_temp[j])

#-----------------------------------------------------------------------------------------------
#-------------------- CONNECTIVITY--------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

        conn = np.zeros((num_elem, num_node_elem), dtype=int)
        addition = [2, 2, 2]

        # initial start
        conn[0][0] = 0;
        conn[0][1] = 2;
        conn[0][2] = 1;

        end_sections = []
        end_count = 0
        for i in range(len(elems)):
            end_count += elems[i]
            end_sections.append(end_count)

        holder = 0
        conn_counter = 0
        origin_counter = 0
        for i in range(num_elem):
            try:
                if i == origin_end[origin_counter + 1]:
                    origin_counter += 1
            except:
                pass
            if i == 0:  # if not i        a = [];   if not a    if len(a) == 0
                holder = 2
                continue
            if i == branch_end[conn_counter]:
                holder = origin_elem[origin_counter]
                conn_counter += 1
            conn[i, :] = [conn[i - 1, 1], conn[i - 1, 1] + 2, conn[i - 1, 2] + 2]
            if piv_log[i] == 'false' and piv_log[i - 1] == 'true':
                holder = conn[i - 1][1]
            for j in range(len(end_sections)):
                if piv_log[i - 1] == 'false' and end_sections[j] == i and piv_temp[j] != 'continue':
                    conn[i][0] = holder
                    break

#-----------------------------------------------------------------------------------------------
#--------------------STIFFNESS------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

        # stiffness matrices
        stiff_temp = []
        mass_temp = []
        lin_mass = []
        lin_stiff = []

        mass_counter = 0
        section_counter = 0
        for k in range(len(data['origins'])):
            for i in range(len(data['origins'][k]['branch'])):
                for j in range(len(data['origins'][k]['branch'][i]['section'])):
                    try:  # linear interpolation stiffness
                        stiff_temp1 = data['origins'][k]['branch'][i]['section'][j]['stiffness']
                        stiff_temp2 = data['origins'][k]['branch'][i]['section'][j]['stiffness2']
                        lin_stiff.append(1)

                        EA = np.linspace(stiff_temp1[0], stiff_temp2[0],
                                         data['origins'][k]['branch'][i]['section'][j]['n_elements'])
                        GAy = np.linspace(stiff_temp1[1], stiff_temp2[1],
                                          data['origins'][k]['branch'][i]['section'][j]['n_elements'])
                        GAz = np.linspace(stiff_temp1[2], stiff_temp2[2],
                                          data['origins'][k]['branch'][i]['section'][j]['n_elements'])
                        GJ = np.linspace(stiff_temp1[3], stiff_temp2[3],
                                         data['origins'][k]['branch'][i]['section'][j]['n_elements'])
                        EIy = np.linspace(stiff_temp1[4], stiff_temp2[4],
                                          data['origins'][k]['branch'][i]['section'][j]['n_elements'])
                        EIz = np.linspace(stiff_temp1[5], stiff_temp2[5],
                                          data['origins'][k]['branch'][i]['section'][j]['n_elements'])

                        for a in range(len(EA)):
                            stiff_temp.append([EA[a], GAy[a], GAz[a], GJ[a], EIy[a], EIz[a]])


                    except:
                        stiff_temp.append(data['origins'][k]['branch'][i]['section'][j]['stiffness'])
                        lin_stiff.append(0)
                        #elems = data['origins'][k]['branch'][i]['section'][j]['n_elements']
                        #for a in range(elems):


#-----------------------------------------------------------------------------------------------
#--------------------MASS PROPERTIES------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
        for k in range(len(data['origins'])):
            for i in range(len(data['origins'][k]['branch'])):
                for j in range(len(data['origins'][k]['branch'][i]['section'])):
                    m_bar1 = data['origins'][k]['branch'][i]['section'][j]['mass']
                    xcg1 = np.array(data['origins'][k]['branch'][i]['section'][j]['xcg'])
                    mass_temp1 = data['origins'][k]['branch'][i]['section'][j]['inertia']
                    try:
                        m_bar2 = data['origins'][k]['branch'][i]['section'][j]['mass2']
                    except:
                        m_bar2 = m_bar1
                    try:
                        xcg2 = np.array(data['origins'][k]['branch'][i]['section'][j]['xcg2'])
                        xcg_same = False
                    except:
                        xcg2 = xcg1
                        xcg_same = True
                    try:
                        mass_temp2 = data['origins'][k]['branch'][i]['section'][j]['inertia2']
                    except:
                        mass_temp2 = mass_temp1

                    interp_hold = True
                    if m_bar2 == m_bar1 and xcg_same == True:
                        for a in range(len(mass_temp2)):
                            if mass_temp2[i] == mass_temp1[a]:
                                if a == len(mass_temp2) - 1:
                                    interp_hold = False
                                    lin_mass.append(0)


                    if interp_hold == True:
                        lin_mass.append(1)

                    m_bar_temp = np.linspace(m_bar1, m_bar2,
                                             data['origins'][k]['branch'][i]['section'][j]['n_elements'])
                    Ixx = np.linspace(mass_temp1[0][0], mass_temp2[0][0],
                                      data['origins'][k]['branch'][i]['section'][j]['n_elements'])
                    Iyy = np.linspace(mass_temp1[1][1], mass_temp2[1][1],
                                      data['origins'][k]['branch'][i]['section'][j]['n_elements'])
                    Izz = np.linspace(mass_temp1[2][2], mass_temp2[2][2],
                                      data['origins'][k]['branch'][i]['section'][j]['n_elements'])
                    xcg_temp = np.linspace(xcg1, xcg2, data['origins'][k]['branch'][i]['section'][j]['n_elements'])

                    for a in range(len(Ixx)):
                        inertia_f = np.diag([Ixx[a], Iyy[a], Izz[a]])
                        m_bar = m_bar_temp[a]
                        xcg = xcg_temp[a]
                        mass_temp.append(matrix.mass_matrix_generator(m_bar, xcg, inertia_f))
                    section_counter += 1

#-----------------------------------------------------------------------------------------------
#--------------------MASS & STIFF MATRICES------------------------------------------------------
#-----------------------------------------------------------------------------------------------
        stiffness = np.zeros((len(stiff_temp), 6, 6))

        stiff_mult = 1  # .36 for agreement with the lumped mass and dynamic test = .63, .45 for frequency =.7
        if stiff_mult != 1:
            print('SITFFNESS MULT IN USE')

        for i in range(len(stiffness)):
            # stiff_temp[i][0] = stiff_temp[i][0] * stiff_mult
            # stiff_temp[i][1] = stiff_temp[i][1] * stiff_mult
            # stiff_temp[i][2] = stiff_temp[i][2] * stiff_mult
            # stiff_temp[i][3] = stiff_temp[i][3] * stiff_mult
            # stiff_temp[i][4] = stiff_temp[i][4] * stiff_mult
            # stiff_temp[i][5] = stiff_temp[i][5] * stiff_mult
            
            stiff_temp[i] *= stiff_mult*np.diag(stiff_temp[i])
            # stiffness[i, ...] = np.diag(stiff_temp[i])

        cutoff = 3
        if cutoff > 0:
            print('CUTOFF IN USE')

        mass = np.zeros((len(mass_temp), 6, 6))
        for i in range(len(mass)):
            if mass_temp[i][3][3] < cutoff:
                mass_temp[i][3][3] = cutoff
                mass_temp[i][4][4] = cutoff * .5
                mass_temp[i][5][5] = cutoff * .5

            mass[i, ...] = mass_temp[i]

#-----------------------------------------------------------------------------------------------
#--------------------ELEM STIFFNESS & MASS------------------------------------------------------
#-----------------------------------------------------------------------------------------------
        elem_stiffness = np.zeros(num_elem, dtype=int)

        elem_counter = 0
        lin_counter = 0
        # elem_stiffness writing
        for i in range(len(elem_stiffness)):
            elem_stiffness[i] = elem_counter
            if lin_stiff[lin_counter] == 1 and i < end_sections[lin_counter] - 1:  # if linear interpolation is used
                elem_counter += 1

            elif i == end_sections[lin_counter] - 1:
                lin_counter += 1

            for j in range(len(end_sections)):
                if i == end_sections[j] - 1:
                    elem_counter += 1
                    break

        # elem_mass
        elem_mass = np.zeros(num_elem, dtype=int)

        elem_counter = 0
        lin_counter = 0
        for i in range(len(elem_mass)):
            elem_mass[i] = elem_counter
            if lin_mass[lin_counter] == 1 and i < end_sections[lin_counter] - 1:  # if linear interpolation is used
                elem_counter += 1

            elif i == end_sections[lin_counter] - 1:
                lin_counter += 1

            for j in range(len(end_sections)):
                if i == end_sections[j] - 1:
                    elem_counter += 1
                    break

#-----------------------------------------------------------------------------------------------
#--------------------FRAME OF REFERENCE---------------------------------------------------------
#-----------------------------------------------------------------------------------------------
        frame_reference_delta = np.zeros((num_elem, num_node_elem, 3))
        frame = []
        branch_end = []
        branch_count = 0

        for k in range(len(data['origins'])):
            for i in range(len(data['origins'][k]['branch'])):
                for j in range(len(data['origins'][k]['branch'][i]['section'])):
                    frame.append(data['origins'][k]['branch'][i]['section'][j]['frame'])
                    branch_count += data['origins'][k]['branch'][i]['section'][j]['n_elements']

        frame_final = []
        for i in range(len(frame)):
            for j in range(elems[i]):
                frame_final.append(frame[i])

        for i in range(num_elem):
            # for j in range(len(branch_end)):
            #     if i < branch_end[j]:
            if frame_final[i] == 'right':
                frame_reference_delta[i][0] = [-1, 0, 0]
                frame_reference_delta[i][1] = [-1, 0, 0]
                frame_reference_delta[i][2] = [-1, 0, 0]
                continue
            if frame_final[i] == 'left':
                frame_reference_delta[i][0] = [1, 0, 0]
                frame_reference_delta[i][1] = [1, 0, 0]
                frame_reference_delta[i][2] = [1, 0, 0]
                continue
            if frame_final[i] == 'fuselage':
                frame_reference_delta[i][0] = [0, 1, 0]
                frame_reference_delta[i][1] = [0, 1, 0]
                frame_reference_delta[i][2] = [0, 1, 0]
                continue
            if frame_final[i] == 'v_fin':
                frame_reference_delta[i][0] = [-1, 0, 0]
                frame_reference_delta[i][1] = [-1, 0, 0]
                frame_reference_delta[i][2] = [-1, 0, 0]
                continue
            if frame_final[i] == 'fuselage_neg':
                frame_reference_delta[i][0] = [0, -1, 0]
                frame_reference_delta[i][1] = [0, -1, 0]
                frame_reference_delta[i][2] = [0, -1, 0]
                continue

#-----------------------------------------------------------------------------------------------
#--------------------STRUCTURAL TWIST & BOUNDARY CONDITIONS-------------------------------------
#-----------------------------------------------------------------------------------------------
        # structural twist
        structural_twist = np.zeros((num_elem, num_node_elem))

        # boundary condition
        boundary_conditions = np.zeros(num_node, dtype=int)
        boundary_conditions[0] = 1
        bound_temp = []
        for k in range(len(data['origins'])):
            for i in range(len(data['origins'][k]['branch'])):
                for j in range(len(data['origins'][k]['branch'][i]['section'])):
                    bound_temp.append(data['origins'][k]['branch'][i]['section'][j]["open_end"])

        counter_bound_index = 0
        counter_end = 0
        for i in range(num_node):
            if i == counter_end + elems[counter_bound_index] * 2:
                counter_end += elems[counter_bound_index] * 2
                if bound_temp[counter_bound_index] == "true":
                    boundary_conditions[i] = -1
                counter_bound_index += 1

#-----------------------------------------------------------------------------------------------
#--------------------BEAM NUMBER----------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

        beam_number = np.zeros(num_elem, dtype=int)

#-----------------------------------------------------------------------------------------------
#--------------------APPLIED FORCES-------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

        app_forces = np.zeros((num_node, 6))

        def eliptical(y):
            F = np.sqrt((1 - (y ** 2 / 36.175 ** 2)) * 99.55 ** 2)
            return F

        def triangle(y):
            F = -4.312 * y + 156
            return F

        def modified_eliptical(y):
            y = y / 36.175
            F = 132.5 * ((1 - y ** 2) ** (1.5))
            return F

        def bell(y):
            scale = 6033.28
            sigma = 13
            mean = 10
            F = (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-((y ** 2 - mean) / (2 * sigma ** 2))) * scale
            return F

        thrust = 191.2
        nodes_app = [12, 46]  # nacelle nodes 12, 46 for other input file [102, 106 , 136, 140]

        for i in range(len(nodes_app)):
            if nodes_app[i] < 20:
                app_forces[nodes_app[i]][1] = thrust * np.cos(4.2 * np.pi / 180)
                print('APPLIED FORCES IN USE')
                # app_forces[nodes_app[i]][2] = -thrust * np.sin(4.2 * np.pi / 180)

            else:
                app_forces[nodes_app[i]][1] = -thrust * np.cos(4.2 * np.pi / 180)
                # app_forces[nodes_app[i]][2] = -thrust * np.sin(4.2 * np.pi / 180)

#-----------------------------------------------------------------------------------------------
#--------------------LUMPED MASSES--------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
        lumped_mass_nodes = []
        lumped_mass = []

        lumped_mass_position_temp = []
        lump_sec = []
        for k in range(len(data['origins'])):
            for i in range(len(data['origins'][k]['branch'])):
                for j in range(len(data['origins'][k]['branch'][i]['section'])):
                    try:
                        lumped_mass.append(data['origins'][k]['branch'][i]['section'][j]['lumped_mass'])
                        lumped_mass_position_temp.append(
                            data['origins'][k]['branch'][i]['section'][j]['lumped_position'])
                        lump_sec.append(1)
                    except:
                        lump_sec.append(0)

        n_lumped_mass = len(lumped_mass)
        lumped_mass_inertia = np.zeros((n_lumped_mass, 3, 3))
        lumped_mass_position = np.zeros((n_lumped_mass, 3))

        x_distance = np.zeros(len(x))
        y_distance = np.zeros(len(x))
        z_distance = np.zeros(len(x))

        distance = np.zeros(len(x))

        for i in range(len(lumped_mass)):
            for j in range(len(x)):
                x_distance[j] = lumped_mass_position_temp[i][0] - x[j]
                y_distance[j] = lumped_mass_position_temp[i][1] - y[j]
                z_distance[j] = lumped_mass_position_temp[i][2] - z[j]
                distance[j] = np.sqrt(x_distance[j] ** 2 + y_distance[j] ** 2 + z_distance[j] ** 2)
            minimal = min(distance)
            min_index = [i for i, j in enumerate(distance) if j == minimal]  # finds the minimum index to make th
            min_index = str(min_index)[1:-1]
            min_index = min_index.split(", ")
            min_index = int(min_index[0])
            lumped_mass_nodes.append(min_index)

            # finds the element where the node is
            save = []
            for aa in range(num_elem):
                if save != []:
                    break
                for aaa in range(3):
                    if min_index == conn[aa][aaa]:
                        save = aa
                        break

            if frame_final[save] == 'right':
                lumped_mass_position[i, ...] = [y_distance[min_index], -x_distance[min_index], z_distance[min_index]]

            elif frame_final[save] == 'left':
                lumped_mass_position[i, ...] = [-y_distance[min_index], x_distance[min_index], z_distance[min_index]]

            elif frame_final[save] == 'fuselage':
                lumped_mass_position[i, ...] = [x_distance[min_index], y_distance[min_index], z_distance[min_index]]

            elif frame_final[save] == 'fuselage_neg':
                lumped_mass_position[i, ...] = [-x_distance[min_index], -y_distance[min_index], z_distance[min_index]]

            elif frame_final[save] == 'v_fin':
                lumped_mass_position[i, ...] = [-z_distance[min_index], -x_distance[min_index], y_distance[min_index]]

# --------------------------------------------------------------------------------------------------------------------------
# --------------------INPUT WRITE-------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------

        with h5.File(self.route + '/' + self.case_name + '.fem.h5', 'a') as h5file:
            coordinates = h5file.create_dataset('coordinates', data=np.column_stack((x, y, z)))
            conectivities = h5file.create_dataset('connectivities', data=conn)
            num_nodes_elem_handle = h5file.create_dataset(
                'num_node_elem', data=num_node_elem)
            num_nodes_handle = h5file.create_dataset(
                'num_node', data=num_node)
            num_elem_handle = h5file.create_dataset(
                'num_elem', data=num_elem)
            stiffness_db_handle = h5file.create_dataset(
                'stiffness_db', data=stiffness)
            stiffness_handle = h5file.create_dataset(
                'elem_stiffness', data=elem_stiffness)
            mass_db_handle = h5file.create_dataset(
                'mass_db', data=mass)
            mass_handle = h5file.create_dataset(
                'elem_mass', data=elem_mass)
            frame_of_reference_delta_handle = h5file.create_dataset(
                'frame_of_reference_delta', data=frame_reference_delta)
            structural_twist_handle = h5file.create_dataset(
                'structural_twist', data=structural_twist)
            bocos_handle = h5file.create_dataset(
                'boundary_conditions', data=boundary_conditions)
            beam_handle = h5file.create_dataset(
                'beam_number', data=beam_number)
            app_forces_handle = h5file.create_dataset(
                'app_forces', data=app_forces)
            lumped_mass_nodes_handle = h5file.create_dataset(
                'lumped_mass_nodes', data=lumped_mass_nodes)
            lumped_mass_handle = h5file.create_dataset(
                'lumped_mass', data=lumped_mass)
            lumped_mass_inertia_handle = h5file.create_dataset(
                'lumped_mass_inertia', data=lumped_mass_inertia)
            lumped_mass_position_handle = h5file.create_dataset('lumped_mass_position', data=lumped_mass_position)

        self.num_elem = num_elem
        self.num_node = num_node
        self.x = x
        self.y = y
        self.z = z
        self.conn = conn

    def Excel_Calc(self):
        df = pd.read_excel(self.data, sheet_name= 'Master')
        rows = len(df.values)

# -----------------------------------------------------------------------------------------------------
# ---------------------COORDINATES/ STIFFNESS_DB/ MASS_DB----------------------------------------------
# -----------------------------------------------------------------------------------------------------

        # function used to linearly interpolate the various paramaters
        def interp(start, end, n, flag):
            """function used to linearly space the various paramaters

            Args:
                start (float): starting value
                end (float): ending value
                n_nodes (float): number of nodes of the given section
                flag (0 or 1): conditional statement for starting section and subsequent sections. 1 is needed for n_node properties 

            Returns:
                list[float]: list of neccesary values at each node or element 
            """
            if flag == 0:
                val = np.linspace(start, end, n)
            elif flag == 1:
                start = (end - start) / n + start
                val = np.linspace(start, end, n)
            return val

        # defining empty lists to append to
        empt = []
        x = np.array(empt)
        y = np.array(empt)
        z = np.array(empt)
        

# calculating the coordinates of the nodes  
        num_node_elem = 3
        num_elem = 0
        num_node = 0 
        
        mass_temp = []
        stiff_temp = []
        mass_sym_hold = [] 
        
        start = 1 # flag if we are at the start of a section

        for i in range(1, rows):
            # defining all of the input values from the excel sheet
            
            x_1 = df.values[i, 2]
            if str(x_1) == 'nan':
                break
            x_2 = df.values[i, 5]
            y_1 = df.values[i, 3]
            y_2 = df.values[i, 6]
            z_1 = df.values[i, 4]
            z_2 = df.values[i, 7]
            
            EA_1 = df.values[i, 8]
            EA_2 = df.values[i, 14]
            GAy_1 = df.values[i, 9]
            GAy_2 = df.values[i, 15]
            GAz_1 = df.values[i, 10]
            GAz_2 = df.values[i, 16]
            GJ_1 = df.values[i, 11]
            GJ_2 = df.values[i, 17]
            EIy_1 = df.values[i, 12]
            EIy_2 = df.values[i, 18]
            EIz_1 = df.values[i, 13]
            EIz_2 = df.values[i, 19]

            mass_1 = df.values[i, 20]
            mass_2 = df.values[i, 21]
            
            xcg1 = df.values[i, 22]
            ycg1 = df.values[i, 23]
            zcg1 = df.values[i, 24]
            xcg2 = df.values[i, 25]
            ycg2 = df.values[i, 26]
            zcg2 = df.values[i, 27]

            Ixx_1 = df.values[i, 28]
            Iyy_1 = df.values[i, 29]
            Izz_1 = df.values[i, 30]
            Ixx_2 = df.values[i, 31]
            Iyy_2 = df.values[i, 32]
            Izz_2 = df.values[i, 33]
            
            try:
                beam_opt = str(df.values[i, 47]).lower()
            except:
                beam_opt = 'no'
                
            # checks to see if current row has symmetric flag and is the start of a new section
            if str(df.values[i, 34]).lower() == "yes"  and start == 1:
                start = 0
                if i == 1:
                    node_hold = 1
                    elem_hold = num_elem
                
                else: 
                    node_hold = num_node
                    elem_hold = num_elem
                    
            n_elem = df.values[i, 1]
            if beam_opt == 'yes':
                for ii in range(num_elem, num_elem + n_elem):
                    self.beam_opt_elems.append(int(ii))

            num_elem += n_elem 

            if i == 1:
                n_nodes = df.values[i, 1] * 2 + 1
                flag = 0
            
            else:
                n_nodes = df.values[i, 1] * 2 
                flag = 1
                

            if beam_opt == 'yes':
                for ii in range(num_node, num_node + n_nodes):
                        self.beam_opt_nodes.append(int(ii))

            num_node += n_nodes

            
            x = np.append(x, interp(x_1, x_2, n_nodes, flag))      
            y = np.append(y, interp(y_1, y_2, n_nodes, flag))      
            z = np.append(z, interp(z_1, z_2, n_nodes, flag))  
            
            flag = 0 # flag is - as all of these properties are for n_elem not n_nodes
            EA = interp(EA_1, EA_2, n_elem, flag)      
            GAy = interp(GAy_1, GAy_2, n_elem, flag)      
            GAz = interp(GAz_1, GAz_2, n_elem, flag)  
            GJ = interp(GJ_1, GJ_2, n_elem, flag)      
            EIy = interp(EIy_1, EIy_2, n_elem, flag)      
            EIz = interp(EIz_1, EIz_2, n_elem, flag)  

            mass =  interp(mass_1, mass_2, n_elem, flag)  
            xcg_temp = interp(xcg1, xcg2, n_elem, flag)
            ycg_temp = interp(ycg1, ycg2, n_elem, flag)
            zcg_temp = interp(zcg1, zcg2, n_elem, flag)
            Ixx =  interp(Ixx_1, Ixx_2, n_elem, flag)  
            Iyy =  interp(Iyy_1, Iyy_2, n_elem, flag)  
            Izz =  interp(Izz_1, Izz_2, n_elem, flag)  

            for qq in range(len(EA)):
                stiff_temp.append(np.diag([EA[qq], GAy[qq], GAz[qq], GJ[qq], EIy[qq], EIz[qq]]))
                mass_temp.append(matrix.mass_matrix_generator(mass[qq], np.array([xcg_temp[qq],ycg_temp[qq],zcg_temp[qq]]), np.diag([Ixx[qq], Iyy[qq], Izz[qq]])))
                # holds the values of the mass matrix if the section was in the'left' frame of reference
                mass_sym_hold.append(matrix.mass_matrix_generator(mass[qq], np.array([xcg_temp[qq], -ycg_temp[qq], zcg_temp[qq]]), np.diag([Ixx[qq], Iyy[qq], Izz[qq]])))

            # symmetric data values added here 
            if i != rows - 1:
                # checks for discountinous and that sym flag is yes in order to add sym values    
                if str(df.values[i, 34]).lower() == 'yes': 
                    if beam_opt == 'yes':
                        self.beam_opt_sym = True
                    
                    if [df.values[i+1, 2], df.values[i+1, 3], df.values[i+1, 4]] != [df.values[i, 5], df.values[i, 6], df.values[i, 7]]:
                        for ii in range(node_hold, len(x)):
                            x = np.append(x, x[ii])
                            y = np.append(y, -y[ii])
                            z = np.append(z, z[ii])
                            num_node += 1
                        
                        for ii in range(elem_hold, len(stiff_temp)):
                            stiff_temp.append(stiff_temp[ii])
                            mass_temp.append(mass_sym_hold[ii])
                            mass_sym_hold.append(0)
                            num_elem += 1
                        start = 1
                    
            # if the last section has the sym flag        
            if i == rows -1 and str(df.values[i, 34]).lower() == 'yes': 
                if beam_opt == 'yes':
                    self.beam_opt_sym = True
                
                for ii in range(node_hold, len(x)):
                    x = np.append(x, x[ii])
                    y = np.append(y, -y[ii])
                    z = np.append(z, z[ii])
                    num_node += 1
                
                for ii in range(elem_hold, len(stiff_temp)):
                    stiff_temp.append(stiff_temp[ii])
                    mass_temp.append(mass_sym_hold[ii])
                    mass_sym_hold.append(0)
                    num_elem += 1
                start = 1
                    
        # finalizing stiffness and mass matrices to be sent to .h5
        stiffness = np.zeros((len(stiff_temp), 6, 6))
        for i in range(len(stiff_temp)):
            stiffness[i, ...] = stiff_temp[i]
            
        mass = np.zeros((len(mass_temp), 6, 6))
        for i in range(len(mass)):
            mass[i, ...] = mass_temp[i]
            
# -----------------------------------------------------------------------------------------------------
#----------------------ELEM STIFFNESS AND MASS---------------------------------------------------------
# -----------------------------------------------------------------------------------------------------

        elem_stiffness = np.linspace(0, num_elem - 1, num_elem, dtype= int)
        elem_mass = np.linspace(0, num_elem - 1, num_elem, dtype = int)       
            
# -----------------------------------------------------------------------------------------------------
#----------------------CONNECTIVITY--------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------
        position = np.column_stack((x, y, z))
        
        connectivites = np.zeros((num_elem, num_node_elem), dtype = int)
        connectivites[0, ...] = [0, 2, 1]
        
        # index that keeps track of the current element
        count = 0
        # node that serves as the pivot point when there is a discontinuity
        hold = 0
        total = 0
        hold_sym = 0
        start = 1
        for i in range(1, rows):
            if i == 6:
                print('Hi') 
            # saving the first point of sym section
            if start == 1 and str(df.values[i, 34]).lower() == 'yes':
                init_loc = [df.values[i, 2], df.values[i, 3], df.values[i, 4]]
                start = 0
            # keeping track of how many sym elements there are    
            if str(df.values[i, 34]).lower() == 'yes':
                hold_sym += df.values[i, 1]

            # finding the initial connectivity point of each section
            for j in range(len(position)):
                if df.values[i, 2] == position[j, 0] and df.values[i, 3] == position[j, 1] and df.values[i, 4] == position[j, 2]:
                    hold = j
                    break
            n_elements = df.values[i,1]
            
            if str(n_elements) == 'nan':
                break
            
            for j in range(n_elements):
                # keeping track of number of elements for sym
                total += 1
                count += 1 
                
                # skip initial point as it is already defined
                if i == 1 and j == 0:
                    count -= 1
                    continue
            
                # defines the first point as hold as this is the only point that will break the pattern 
                if j == 0:    
                    connectivites[count, ...] = [hold, connectivites[count-1, 1] + 2, connectivites[count-1, 2] + 2]
                
                if j != 0:
                    connectivites[count, ...] = [connectivites[count-1, 1], connectivites[count-1, 1] + 2, connectivites[count-1, 2] + 2]
                    
                if j == n_elements - 1:
                    # checks for disontinuity and that there is a symmetric flag
                    if i != rows - 1 and str(df.values[i, 34]).lower() == "yes":
                        if [df.values[i+1, 2], df.values[i+1, 3], df.values[i+1, 4]] != [df.values[i, 5], df.values[i, 6], df.values[i, 7]]:
                            # start flag indicating that the next sym section is the start
                            start = 1
                            # makes the hold for the sym case by setting the init_loc of the section
                            for j in range(len(position)):
                                if init_loc[0] == position[j, 0] and -init_loc[1] == position[j, 1] and init_loc[2] == position[j, 2]:
                                    hold = j
                                    break
                            for j in range(hold_sym):
                                count += 1
                                if j == 0:    
                                    connectivites[count, ...] = [hold, connectivites[count-1, 1] + 2, connectivites[count-1, 2] + 2]
                                
                                if j != 0:
                                    connectivites[count, ...] = [connectivites[count-1, 1], connectivites[count-1, 1] + 2, connectivites[count-1, 2] + 2]
                            # reset sym_count for next sym operation
                            hold_sym = 0
                    
                    # triggers sym operation for the last section
                    if i == rows - 1 and str(df.values[i, 34]).lower() == "yes":
                            # start flag indicating that the next sym section is the start
                            start = 1
                            # makes the hold for the sym case by setting the init_loc of the section
                            for j in range(len(position)):
                                if init_loc[0] == position[j, 0] and -init_loc[1] == position[j, 1] and init_loc[2] == position[j, 2]:
                                    hold = j
                                    break
                            for j in range(hold_sym):
                                count += 1
                                if j == 0:    
                                    connectivites[count, ...] = [hold, connectivites[count-1, 1] + 2, connectivites[count-1, 2] + 2]
                                
                                if j != 0:
                                    connectivites[count, ...] = [connectivites[count-1, 1], connectivites[count-1, 1] + 2, connectivites[count-1, 2] + 2]
                            # reset sym_count for next sym operation
                            hold_sym = 0
                                
# -----------------------------------------------------------------------------------------------------
#----------------------FRAME OF REFERENCE--------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------
        
        frame_reference_delta = np.zeros((num_elem, num_node_elem, 3))
        frame_temp = []
        ref_look_up = []
        sym_count = 0
        origin_flag = 0 # if the starting coordinate of a section has a y component of 0, flag is (1) so that it is known to change frame from 'right' to 'left'
        start = 1 # flag if we are at the start of a section
        for i in range(1, rows):
            n_elem = df.values[i, 1]
            if str(n_elem) == 'nan':
                break
            frame = str(df.values[i, 35]).lower()                  
            
            # raise origin flag if the starting coordinate has a y component of 0 and it is the start of a new section 
            if str(df.values[i, 34]).lower() == 'yes' and start == 1: 
                if df.values[i, 3] ==  0:
                    origin_flag = 1
                else:
                    origin_flag = 0
                start = 0
                sym_count = 0
            
            for j in range(n_elem):
    
                if frame == 'right':
                    right = [[-1, 0, 0], [-1, 0, 0], [-1, 0, 0]]
                    frame_temp.append(right)
                    sym_count += 1
                    ref_look_up.append('right')
                    
                if frame == 'left':
                    left = [[1, 0, 0], [1, 0, 0], [1, 0, 0]]
                    frame_temp.append(left)
                    sym_count += 1
                    ref_look_up.append('left')

                if frame == 'fuselage':
                    fuse = [[0, 1, 0], [0, 1, 0], [0, 1, 0]]
                    frame_temp.append(fuse)
                    sym_count += 1
                    ref_look_up.append('fuselage')
                    
                if frame == 'v_fin':
                    v_fin = [[-1, 0, 0], [-1, 0, 0], [-1, 0, 0]]
                    frame_temp.append(v_fin)
                    sym_count += 1
                    ref_look_up.append('v_fin')
                    
                if frame == 'fuselage_neg':
                    fuse_neg = [[0, -1, 0], [0, -1, 0], [0, -1, 0]]
                    frame_temp.append(fuse_neg)
                    sym_count += 1
                    ref_look_up.append('fuselage_neg')
             
            if origin_flag == 1:
                hold = ref_look_up[-1]
                frame_hold = frame_temp[-1]
                temp = np.zeros((3,3))
                temp[0, ...] = frame_hold[0]
                temp[1, ...] = frame_hold[1]
                temp[2, ...] = frame_hold[2]
                
            # checks for discountinuity and sym flag to add sym section     
            if i != rows -1 and str(df.values[i, 34]).lower() == 'yes': 
                if [df.values[i+1, 2], df.values[i+1, 3], df.values[i+1, 4]] != [df.values[i, 5], df.values[i, 6], df.values[i, 7]]:
                    start = 1
                    if hold == 'right' and origin_flag == 1:
                        hold = 'left'
                        temp *= -1
                            
                    if origin_flag == 1:
                        for ii in range(sym_count):
                            frame_temp.append(temp)
                            ref_look_up.append(hold)
                            
                    else:
                        count = 0
                        for ii in range(int(sym_count), 0, -1):
                            
                            if ref_look_up[-ii-count] == 'left' :
                                ref_look_up.append('right')
                                cool = [[-1,0,0], [-1,0,0], [-1,0,0]]
                                frame_temp.append(cool)
                                
                            elif ref_look_up[-ii-count] == 'right' :
                                ref_look_up.append('left')
                                cool = [[1,0,0], [1,0,0], [1,0,0]]
                                frame_temp.append(cool)
                            else:
                                ref_look_up.append(ref_look_up[-ii-count])
                                frame_temp.append(frame_temp[-ii-count])
                            count += 1
                        
                            
            # if the last section has the sym flag        
            if i == rows - 1 and str(df.values[i, 34]).lower() == 'yes': 
                start = 1
                if hold == 'right' and origin_flag == 1:
                    hold = 'left'
                    temp *= -1
                        
                if origin_flag == 1:
                    for ii in range(sym_count):
                        frame_temp.append(temp)
                        ref_look_up.append(hold)
                            
                else:
                    count = 0
                    for ii in range(int(sym_count), 0, -1):
                        
                        if ref_look_up[-ii-count] == 'left' :
                            ref_look_up.append('right')
                            cool = [[-1,0,0], [-1,0,0], [-1,0,0]]
                            frame_temp.append(cool)
                                
                        elif ref_look_up[-ii-count] == 'right' :
                            ref_look_up.append('left')
                            cool = [[1,0,0], [1,0,0], [1,0,0]]
                            frame_temp.append(cool)
                        else:
                            ref_look_up.append(ref_look_up[-ii-count])
                            frame_temp.append(frame_temp[-ii-count])
                        count += 1
                
        # finalizing frame of reference matrix 
        for i in range(len(frame_temp)):
            frame_reference_delta[i, ...] = frame_temp[i]
        
#-----------------------------------------------------------------------------------------------
#--------------------STRUCTURAL TWIST & BOUNDARY CONDITIONS-------------------------------------
#-----------------------------------------------------------------------------------------------

        # stuctural twist
        structural_twist = np.zeros((num_elem, num_node_elem))
        
        # boundary conditions
        boundary_conditions = np.zeros(num_node, dtype=int)
        boundary_conditions[0] = 1
        
        start = []
        end = []
        elem_end = []
        for i in range(len(connectivites)):
            start.append(connectivites[i][0])
            end.append(connectivites[i][1])
            elem_end.append(i)
        
        boundary_conditions[0] = 1
        for i in range(len(end)):
            for j in range(len(start)):
                if end[i] != start[j]:
                    if j == len(start)-1:
                        boundary_conditions[end[i]] = -1
                else:
                    break
    
#-----------------------------------------------------------------------------------------------
#--------------------BEAM NUMBER----------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

       
        beam_number = np.zeros(num_elem, dtype=int)
        
#-----------------------------------------------------------------------------------------------
#--------------------APPLIED FORCES-------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
        
        df_2 = pd.read_excel(self.data, sheet_name= 'Applied_Forces_Lumped_Masses')
        rows = len(df_2.values)
        
        # make dictionary
        force = []
        app_loc = []
        
        for i in range(1, rows):
            if str(df_2.values[i,0]).lower() != 'nan':
                try:
                    force.append([df_2.values[i, 0], df_2.values[i, 1], df_2.values[i, 2], df_2.values[i, 3], df_2.values[i, 4], df_2.values[i, 5]])
                    app_loc.append([df_2.values[i, 6], df_2.values[i, 7], df_2.values[i, 8]])
                except:
                    continue

        app_forces = np.zeros((num_node, 6))
        self.app_nodes = []
        
        distance = np.zeros(len(x))

        for i in range(len(force)):
            for j in range(len(x)):
                x_distance = app_loc[i][0] - x[j]
                y_distance = app_loc[i][1] - y[j]
                z_distance = app_loc[i][2] - z[j]
                distance[j] = np.sqrt(x_distance ** 2 + y_distance ** 2 + z_distance ** 2)
            minimal = min(distance)
            min_index = [i for i, j in enumerate(distance) if j == minimal]  # finds the minimum index to make th
            min_index = str(min_index)[1:-1]
            min_index = min_index.split(", ")
            min_index = int(min_index[0])
            app_forces[min_index, ...] = force[i]
        
        # saves nodes where there is a y force as these nodes are needed for static trim
        for i in range(len(app_forces)):
            if app_forces[i, 1] != 0:
                self.app_nodes.append(i)

                 
#-----------------------------------------------------------------------------------------------
#--------------------LUMPED MASSES--------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
        
        # lumped masses input
        lumped_mass_nodes = []
        
        lumped_mass = []
        lumped_loc = []
        lumped_inertia = []
        
        for i in range(1, rows):
            if str(df_2.values[i,12]).lower() != 'nan':
                lumped_mass.append([df_2.values[i, 12]])
                lumped_loc.append([df_2.values[i, 13], df_2.values[i, 14], df_2.values[i, 15]])
                lumped_inertia.append([df_2.values[i, 16], df_2.values[i, 17], df_2.values[i, 18]])
        
        n_lumped_mass = len(lumped_mass)
        lumped_mass_inertia = np.zeros((n_lumped_mass, 3, 3))
        lumped_mass_position = np.zeros((n_lumped_mass, 3))
        
        for i in range(len(lumped_inertia)):
            lumped_mass_inertia[i, ...] = lumped_inertia[i]

        x_distance = np.zeros(len(x))
        y_distance = np.zeros(len(x))
        z_distance = np.zeros(len(x))

        distance = np.zeros(len(x))

        for i in range(len(lumped_mass)):
            for j in range(len(x)):
                x_distance[j] = lumped_loc[i][0] - x[j]
                y_distance[j] = lumped_loc[i][1] - y[j]
                z_distance[j] = lumped_loc[i][2] - z[j]
                distance[j]  = np.sqrt(x_distance[j]**2 + y_distance[j]**2 + z_distance[j]**2)
            minimal = min(distance)
            min_index = [i for i, j in enumerate(distance) if j == minimal] # finds the minimum index
            min_index = str(min_index)[1:-1]
            min_index = min_index.split(", ")
            min_index = int(min_index[0])
            lumped_mass_nodes.append(min_index)

            # finds the element where the node is
            save = [] 
            for aa in range(num_elem):
                if save != []:
                    break
                for aaa in range(3):
                    if min_index == connectivites[aa][aaa]:
                        save = aa
                        break

            # right = [-1, 0, 0]
            if ref_look_up[save] == 'right':
                lumped_mass_position[i, ...] = [y_distance[min_index], -x_distance[min_index] , z_distance[min_index]]

            # left = [1, 0, 0]
            elif ref_look_up[save] == 'left':
                lumped_mass_position[i, ...] = [-y_distance[min_index], x_distance[min_index] , z_distance[min_index]]

            # fuselage = [0, 1, 0]
            elif ref_look_up[save] == 'fuselage':
                lumped_mass_position[i, ...] = [x_distance[min_index] , y_distance[min_index], z_distance[min_index]]

            # fuselage_neg = [0, -1, 0]
            elif ref_look_up[save] == 'fuselage_neg':
                lumped_mass_position[i, ...] = [-x_distance[min_index] , -y_distance[min_index], z_distance[min_index]]

            # right = [-1, 0, 0]
            elif ref_look_up[save] == 'v_fin':
                lumped_mass_position[i, ...] = [-z_distance[min_index], -x_distance[min_index], y_distance[min_index]]

# --------------------------------------------------------------------------------------------------------------------------
# --------------------INPUT WRITE-------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------

        with h5.File(self.route + '/' + self.case_name + '.fem.h5', 'a') as h5file:
            coordinates = h5file.create_dataset('coordinates', data = position)
            conectivities = h5file.create_dataset('connectivities', data = connectivites)
            num_nodes_elem_handle = h5file.create_dataset('num_node_elem', data = num_node_elem)
            num_nodes_handle = h5file.create_dataset('num_node', data = num_node)
            num_elem_handle = h5file.create_dataset('num_elem', data = num_elem)
            stiffness_db_handle = h5file.create_dataset('stiffness_db', data = stiffness)
            stiffness_handle = h5file.create_dataset('elem_stiffness', data = elem_stiffness)
            mass_db_handle = h5file.create_dataset('mass_db', data = mass)
            mass_handle = h5file.create_dataset('elem_mass', data = elem_mass)
            frame_of_reference_delta_handle = h5file.create_dataset('frame_of_reference_delta', data=frame_reference_delta)
            structural_twist_handle = h5file.create_dataset('structural_twist', data = structural_twist)
            bocos_handle = h5file.create_dataset('boundary_conditions', data=boundary_conditions)
            beam_handle = h5file.create_dataset('beam_number', data=beam_number)
            app_forces_handle = h5file.create_dataset('app_forces', data=app_forces)
            lumped_mass_nodes_handle = h5file.create_dataset('lumped_mass_nodes', data=lumped_mass_nodes)
            lumped_mass_handle = h5file.create_dataset('lumped_mass', data=lumped_mass)
            lumped_mass_inertia_handle = h5file.create_dataset('lumped_mass_inertia', data=lumped_mass_inertia)
            lumped_mass_position_handle = h5file.create_dataset('lumped_mass_position', data=lumped_mass_position)
            
        self.num_elem = num_elem
        self.num_node = num_node
        self.x = x
        self.y = y
        self.z = z
        self.conn = connectivites
        self.ref_lookup = ref_look_up
        self.boundary = boundary_conditions
        self.lumped_loc = lumped_loc