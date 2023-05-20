import ctypes as ct
import sys
import numpy as np

import sharpy.solvers.staticcoupled as staticcoupled
import sharpy.aero.utils.mapping as mapping
import sharpy.utils.cout_utils as cout
import sharpy.utils.solver_interface as solver_interface
from sharpy.utils.solver_interface import solver, BaseSolver
import sharpy.utils.settings as settings
import sharpy.utils.algebra as algebra
import sharpy.postproc.beamloads as beamloads
import matplotlib.pyplot as plt
import sharpy.structure.models.beam as beam
import sympy as sp
import sharpy

@solver
class Beam_Opt(BaseSolver):
    """
    Beam_Opt is cool
    """
    solver_id = 'BeamOpt'
    solver_classification = 'structural'

    settings_types = dict()
    settings_default = dict()
    settings_description = dict()

    # width of the spar 
    settings_types['B'] = 'float'
    settings_default['B'] =  0.2999980281284993 

    # height of the spar
    settings_types['H'] = 'float'
    settings_default['H'] = 0.11706065387250503

    settings_types['density_spar'] = 'float'
    settings_default['density_spar'] = 2720

    settings_types['min_thickness'] = 'float'
    settings_default['min_thickness'] = .00001
    
    settings_types['og_thickness'] = 'float'
    settings_default['og_thickness'] = 0.000132331767098875

    settings_types['E'] = 'float'
    settings_default['E'] = 8817987907.14019
    settings_description['E'] = 'Youngs Modulus'

    settings_types['G'] = 'float'
    settings_default['G'] = 17264944659.8733
    settings_description['G'] = 'Shear modulus of spar material'

    # EI of everything else other than the spar
    settings_types['EI_everything'] = 'float'
    settings_default['EI_everything'] = 2999.99999999999994

    # GJ of everything else other than the spar
    settings_types['GJ_everything'] = 'float'
    settings_default['GJ_everything'] = 14999.9999999999998

    # GJ of everything else other than the spar
    settings_types['Mass_Everything'] = 'float'
    settings_default['Mass_Everything'] = .45 
    
    settings_types['poissons_ratio'] = 'float'
    settings_default['poissons_ratio'] = .19
    
    # ultimate shear strength of material. Value used from Solar Impulse model. Need to change if serious about looking at check
    settings_types['Ultimate_Shear'] = 'float'
    settings_default['Ultimate_Shear'] = .6 * 80e6

    
#--------------------------------------------------------------------------
#--------------------END OF SETTINGS SUBJECT TO CHANGE---------------------    
#--------------------------------------------------------------------------

    settings_types['Spar couple H'] = 'float'
    settings_default['Spar couple H'] = .86
    settings_description['Spar couple H'] = 'Height of Centroid of cap relative to thickness'

    settings_types['Max Strain'] = 'float'
    settings_default['Max Strain'] = .0001
    
    settings_types['Beam_Opt_Nodes'] = 'list(float)'
    settings_default['Beam_Opt_Nodes'] = []

    settings_types['Beam_Opt_Elems'] = 'list(float)'
    settings_default['Beam_Opt_Elems'] = []
    
    settings_types['Sym'] = 'bool'
    settings_default['Sym'] = False

    settings_types['sections'] = 'list(float)'
    settings_default['sections'] = [0, 24.824, 24.824, 36.174]
    settings_description['sections'] = 'Start and End points of each section'

    settings_types['spar_width'] = 'list(float)'
    settings_default['spar_width'] = [.986, .986, .986, .508]

    settings_types['y_spar_lookup'] = 'list(float)'
    settings_default['y_spar_lookup'] = [0, .507, 1.621, 2.735, 4.536, 7.239, 9.940, 12.363, 14.897, 16.219, 18.355, 20.980, 23.543, 24.824, 26.108, 28.675, 30.114, 31.553, 33.094, 34.634, 35.404, 36.174]

    settings_types['up_low_thick'] = 'list(float)'
    settings_default['up_low_thick'] = [1.59,1.17, 0.95, 0.73, 0.73, 0.73, 0.73, 0.73, 0.73, 0.73, 0.64, 0.64, 0.52, 0.52, 0.52, 0.52, 0.52, 0.43, 0.43, 0.31, 0.31, 0.31 ]

    settings_types['fwd_aft_thick'] = 'list(float)'
    settings_default['fwd_aft_thick'] = [1.59, 1.17, 0.95, 0.73, 0.73, 0.73, 0.73, 0.73, 0.73, 0.73, 0.64, 0.64, 0.52,
                                        0.52, 0.52, 0.52, 0.52, 0.43, 0.43, 0.31, 0.31, 0.31]

    settings_types['Solar_Mass_Dist'] = 'float'
    settings_default['Solar_Mass_Dist'] = .8095

    settings_types['Cap_Area'] = 'list(float)'
    settings_default['Cap_Area'] = [944, 944, 944, 944, 840, 791, 784, 650, 575, 420, 330, 274, 200, 185, 170, 122, 98, 73, 73, 73, 61, 50 ]

    settings_types['Solar_Coverage'] = 'list(float)'
    settings_default['Solar_Coverage'] = [1, 1, 1, 1, 1, 1, 1, 1, .94, .94, .94, .94, .94, .94, .77, .77, .77, .77, .77, .77, .77, .77]

    settings_types['t/c'] = 'list(float)'
    settings_default['t/c'] = [.16, .16, .16, .164]
    settings_description['t/c'] = 'Thickness per chord '

    settings_types['Density_Laminate'] = 'float'
    settings_default['Density_Laminate'] = 1540
    settings_description['Density_Laminate'] = 'Density of Spar material'

    settings_types['Unit_Chord_Sect_Area'] = 'float'
    settings_default['Unit_Chord_Sect_Area'] = .1043
    settings_description['Unit_Chord_Sect_Area'] = 'Area that scales with chord'

    settings_types['LE_TE_Density'] = 'float'
    settings_default['LE_TE_Density'] = 3.722
    settings_description['LE_TE_Density'] = 'Density of LE and TE material'

    settings_types['LE_TE_Percent_Area'] = 'float'
    settings_default['LE_TE_Percent_Area'] = .596
    settings_description['LE_TE_Percent_Area'] = 'Percentage of the area of the cross section that is not spar'

    settings_types['Misc_Spar_Factor'] = 'float'
    settings_default['Misc_Spar_Factor'] = .59
    settings_description['Misc_Spar_Factor'] = 'Additional mass of the spar [fudge factor]'

    settings_types['print_info'] = 'bool'
    settings_default['print_info'] = True
    settings_description['print_info'] = 'Print info to screen'

    settings_types['solver'] = 'str'
    settings_default['solver'] = ''
    settings_description['solver'] = 'Solver to run in trim routine'

    settings_types['solver_settings'] = 'dict'
    settings_default['solver_settings'] = dict()
    settings_description['solver_settings'] = 'Solver settings dictionary'

    settings_types['Max_iter'] = 'int'
    settings_default['Max_iter'] = 5
    settings_description['Max_iter'] = 'Maximum number of iterations of routine'

    settings_types['fz_tolerance'] = 'float'
    settings_default['fz_tolerance'] = 0.01
    settings_description['fz_tolerance'] = 'Tolerance in vertical force'

    settings_types['fx_tolerance'] = 'float'
    settings_default['fx_tolerance'] = 0.01
    settings_description['fx_tolerance'] = 'Tolerance in horizontal force'

    settings_types['m_tolerance'] = 'float'
    settings_default['m_tolerance'] = 0.01
    settings_description['m_tolerance'] = 'Tolerance in pitching moment'

    settings_types['tail_cs_index'] = 'int'
    settings_default['tail_cs_index'] = 0
    settings_description['tail_cs_index'] = 'Index of control surfaces that move to achieve trim'

    settings_types['thrust_nodes'] = 'list(int)'
    settings_default['thrust_nodes'] = [0]
    settings_description['thrust_nodes'] = 'Nodes at which thrust is applied'

    settings_types['initial_alpha'] = 'float'
    settings_default['initial_alpha'] = 0.
    settings_description['initial_alpha'] = 'Initial angle of attack'

    settings_types['initial_deflection'] = 'float'
    settings_default['initial_deflection'] = 0.
    settings_description['initial_deflection'] = 'Initial control surface deflection'

    settings_types['initial_thrust'] = 'float'
    settings_default['initial_thrust'] = 0.0
    settings_description['initial_thrust'] = 'Initial thrust setting'

    settings_types['initial_angle_eps'] = 'float'
    settings_default['initial_angle_eps'] = 0.05
    settings_description['initial_angle_eps'] = 'Initial change of control surface deflection'

    settings_types['initial_thrust_eps'] = 'float'
    settings_default['initial_thrust_eps'] = 2.
    settings_description['initial_thrust_eps'] = 'Initial thrust setting change'

    settings_types['relaxation_factor'] = 'float'
    settings_default['relaxation_factor'] = 0.2
    settings_description['relaxation_factor'] = 'Relaxation factor'

    settings_types['save_info'] = 'bool'
    settings_default['save_info'] = False
    settings_description['save_info'] = 'Save trim results to text file'

    settings_types['folder'] = 'str'
    settings_default['folder'] = './output/'
    settings_description['folder'] = 'Output location for trim results'

    settings_table = settings.SettingsTable()
    __doc__ += settings_table.generate(settings_types, settings_default, settings_description)

    def __init__(self):
        self.data = None
        self.settings = None
        self.solver = None
        self.Nodal_Forces = 0
        self.Pass = False
        self.shear = 0
        self.moment = 0
        self.torsion = 0
        self.strain = 0
        self.strain_elem = None
        self.strain_compare = None
        self.EI_Opt_Nodal = None
        self.EI_Opt = None
        self.EI_Opt_Old = None
        self.GJ_Nodal = None
        self.GJ = None
        self.iter = 0
        self.thickness_hold = None
        self.Cap_Area = None
        self.A = None

        #look up values
        self.y_spar_lookup = None
        self.H = None
        self.Ixx_up_low = None
        self.Ixx_fwd_aft = None
        self.chord = None
        self.area_up_low = None
        self.area_fwd_aft = None
        self.cap_area_original = None
        self.solar_mass = None
        self.up_low_width = None
        self.up_low_thick = None
        self.fwd_aft_height = None
        self.fwd_aft_thick = None

        #values for new FEM
        self.elem_mass = None
        self.elem_stiffness = None
        self.mass_db = None
        self.stiffness_db = None
        
#----------------------------------------------------------------------------------------------------------------------------
#----------------------------LOADING ALL OF THE DATA AND SETTINGS FOR THIS SOLVER--------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
    def initialise(self, data):
        self.data = data
        # adding new values to save for 
        self.data.strain = [] 
        self.settings = data.settings[self.solver_id]
        settings.to_custom_types(self.settings, self.settings_types, self.settings_default)

        self.solver = solver_interface.initialise_solver(self.settings['solver'])
        self.solver.initialise(self.data, self.settings['solver_settings'])

#----------------------------------------------------------------------------------------------------------------------------
#----------------------------THE ITERATIVE LOOP OF BEAM OPT------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
    def run(self):
        self.solver.run() # first static couple run
        self.local_for_forces() # converts forces out of static coupled to local frame of reference
        iter_count = 0
        while self.Pass == False:
            iter_count += 1
            self.EI_GJ_Opt() # calculates updated stiffness and mass matrices
            self.generate_fem() # updates matrices for solver
            self.solver.run() # runs static coupled with updated stiffnesses
            self.local_for_forces() # converts forces back to local frame again for the checks
            self.checks() # violations such as deflection angle, margin of safety, etc.
            print('Iterative Step: ' + str(iter_count) + ' completed\n') # printing the completed iterative step
        return self.data


#----------------------------------------------------------------------------------------------------------------------------
#----------------------------All FUNCTIONS CALLED DURING THE ITERATIVE LOOP--------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------

    # converts from aircraft frame of reference to local frame of reference
    def nodal_a_for_2_b_for(self, nodal, tstep, filter=np.array([True] * 6)):
        nodal_a = nodal.copy(order='F')
        for i_node in range(self.data.structure.num_node):
            # get master elem and i_local_node
            i_master_elem, i_local_node = self.data.structure.node_master_elem[i_node, :]
            crv = tstep.psi[i_master_elem, i_local_node, :]
            cab = algebra.crv2rotation(crv)
            temp = np.zeros((6,))
            temp[0:3] = np.dot(cab, nodal[i_node, 0:3])
            temp[3:6] = np.dot(cab, nodal[i_node, 3:6])
            for i in range(6):
                if filter[i]:
                    nodal_a[i_node, i] = temp[i]
        return nodal_a

    def local_for_forces(self, tstep=None):
        if tstep is None:
            tstep = self.data.structure.timestep_info[self.data.ts]
        applied_forces = self.data.structure.nodal_b_for_2_a_for(tstep.steady_applied_forces,
                                                                 tstep)

        applied_forces_copy = applied_forces.copy()
        gravity_forces_copy = tstep.gravity_forces.copy()
        for i_node in range(self.data.structure.num_node):
            applied_forces_copy[i_node, 3:6] += np.cross(tstep.pos[i_node, :],
                                                         applied_forces_copy[i_node, 0:3])
            gravity_forces_copy[i_node, 3:6] += np.cross(tstep.pos[i_node, :],
                                                         gravity_forces_copy[i_node, 0:3])

        # forces are in frame of reference a
        together = np.zeros((len(applied_forces_copy), 6))
        for i in range(len(applied_forces_copy)):
            together[i, ...] = applied_forces_copy[i] + gravity_forces_copy[i]

        self.Nodal_Forces = self.nodal_a_for_2_b_for(together, tstep)
        
        if self.Pass == True:
            self.strain = []
            
            nodes = self.settings['Beam_Opt_Nodes']
            elements = self.settings['Beam_Opt_Elems']
            self.shear = np.zeros(len(nodes))
            for i in range(len(nodes)-1, -1, -1):
                if i == len(nodes)-1:
                    self.shear[i] = self.Nodal_Forces[i][2]
                else:
                    self.shear[i] = self.Nodal_Forces[i][2] + self.shear[i+1]

            self.moment =  np.zeros(len(nodes))
            for i in range(len(nodes) - 1, -1, -1):
                if i == len(nodes) - 1:
                    self.moment[i] = 0
                else:
                    back = self.data.structure.timestep_info[-1].pos[i+1]
                    current = self.data.structure.timestep_info[-1].pos[i]
                    self.moment[i] = self.shear[i+1] * (back[1] - current[1]) + self.moment[i+1]

            self.strain = np.zeros(len(nodes))
            elem_stiff = self.data.structure.elem_stiffness
            stiff = self.data.structure.stiffness_db
            elem_stiff_counter = int(0)
            for i in range(len(nodes)):
                self.strain[i] = ((self.moment[i] * self.settings['H']*.5) / stiff[int(elem_stiff[elem_stiff_counter])][4][4])

            self.data.strain.append(self.strain) 
                    

    def EI_GJ_Opt(self):
        
        # values that allow for the evaluation of the entire planform as opposed to just a single wing
        nodes = self.settings['Beam_Opt_Nodes']
        elements = self.settings['Beam_Opt_Elems']
        self.shear = np.zeros(len(nodes))
        for i in range(len(nodes)-1, -1, -1):
            if i == len(nodes)-1:
                self.shear[i] = self.Nodal_Forces[i][2]
            else:
                self.shear[i] = self.Nodal_Forces[i][2] + self.shear[i+1]

        self.moment =  np.zeros(len(nodes))
        for i in range(len(nodes) - 1, -1, -1):
            if i == len(nodes) - 1:
                self.moment[i] = 0
            else:
                back = self.data.structure.timestep_info[-1].pos[i+1]
                current = self.data.structure.timestep_info[-1].pos[i]
                self.moment[i] = self.shear[i+1] * (back[1] - current[1]) + self.moment[i+1]

        self.strain = np.zeros(len(nodes))
        self.EI_Opt_Nodal = np.zeros(len(nodes))
        max_strain = self.settings['Max Strain']

        elem_stiff = self.data.structure.elem_stiffness
        stiff = self.data.structure.stiffness_db
        elem_stiff_counter = int(0)
        for i in range(len(nodes)):
            self.strain[i] = ((self.moment[i] * self.settings['H']*.5) / stiff[int(elem_stiff[elem_stiff_counter])][4][4]) 
            self.EI_Opt_Nodal[i] = (abs(self.moment[i]*self.settings['H']*.5) / max_strain) 
            if ( i % 2) == 0 and i != 0:
                elem_stiff_counter +=1
        
        self.data.strain.append(self.strain) 
        
        I_yy = []
        ei = self.EI_Opt_Nodal
        e = self.settings['E']
        for i in range(len(ei)):
            I_yy.append((ei[i] - self.settings['EI_everything'])/ e) 


        # solving for the thickness of the spar given Iy
        self.thickness = []
        b = self.settings['B'] # width of spar
        h = self.settings['H'] # height of spar
        
        # defining variable name t to solve using sympy
        from sympy.abc import t
        min_thick = self.settings['min_thickness']
        min_bh = min([self.settings['B'], self.settings['H']])
        max_thick = min_bh/2
        for i in range(len(I_yy)):
            if I_yy[i] >= (b*h**3)/12: # if I_yy is greater than solid recatngle of material
                self.thickness.append(0) # I_yy is solid
                self.EI_Opt_Nodal[i] = self.settings['EI_everything'] + e * (b*h**3)/12
            elif I_yy[i] <= (b*h**3)/12 - ((b-2*min_thick)*(h-2*min_thick)**3)/12: # If I_yy is smaller than I_yy with the minimum thickness 
                self.thickness.append(min_thick)  
                self.EI_Opt_Nodal[i] = self.settings['EI_everything'] + e * ((b*h**3)/ 12 - ((b-2*min_thick)*(h-2*min_thick)**3)/12)
            else: # solve for the thickness if it is in between maximum and minimum bounds
                Eqn = sp.Eq((b*h**3)/ 12 - ((b-2*t)*(h-2*t)**3)/12, I_yy[i])
                ans = sp.solve(Eqn)
                if ans[0] >= max_thick: # making sure that thickness does not exceed maximum 
                    self.thickness.append(0)    
                    self.EI_Opt_Nodal[i] = self.settings['EI_everything'] + e * ((b*h**3)/ 12) 
                else:
                    self.thickness.append(ans[0])
                    self.EI_Opt_Nodal[i] = self.settings['EI_everything'] + e * ((b*h**3)/ 12 - ((b-2*ans[0])*(h-2*ans[0])**3)/12)

        # invokes non-increasing thickness criteria
        for i in range(len(self.thickness)):
            try:
                if self.thickness[i] < self.thickness[i+1] and self.thickness[i] != 0: # zero indicates a solid piece do not need to change this 
                    hold = self.thickness[i+1]
                    for i in range(i + 1):
                        self.thickness[i] = hold
                        I_yy[i] = (b*h**3)/ 12 - ((b-2*self.thickness[i])*(h-2*self.thickness[i])**3)/12
                        self.EI_Opt_Nodal[i] = self.settings['EI_everything'] + e * I_yy[i] 
            except:
                pass

       # Average of EI to update elemental stiffness matrix
        self.EI_Opt = []
        for i in range(0, len(self.EI_Opt_Nodal) - 2, 2):
            self.EI_Opt.append(abs((self.EI_Opt_Nodal[i] + self.EI_Opt_Nodal[i + 1] + self.EI_Opt_Nodal[i + 2]) / 3))

        # section where I can update the width.
        # NEED to know original dimensions of spar cap. Specifically the height as this
        # will be used to backsolve the updated width which feeds into GJ
        if self.thickness_hold == None:
            self.thickness_hold = self.settings['og_thickness']

        # calculating new width
        # change this when you have the correct thickness dimension
        # for i in range(len(self.area_hold)):
        #     width_new = self.Cap_Area[i]/(4*.002)
        #     width_old = self.area_hold[i]/(4*.002)
        #     width_diff = 2*(width_new - width_old) # accounts for both sides width diff
        #     self.up_low_width
        #     self.up_low_width[i] = self.up_low_width[i] + width_diff

        # updating area hold for next iteration
        self.thickness_hold = self.thickness


        # GJ calculation
        self.GJ_Nodal = []
        self.area = []
        g = self.settings['G'] 
        
        for i in range(len(self.thickness)):
            if self.thickness[i] == 0:
               self.area.append(b*h)
               j = (b*h)/12 * (b**2+h**2)
               GJ_temp = g * j
               self.GJ_Nodal.append(self.settings['GJ_everything'] + GJ_temp) 
            else:   
                area = b * h - ((h-2*self.thickness[i])*(b-2*self.thickness[i]))
                self.area.append(area)
                GJ_temp =  (4 * (b*h) ** 2) / ((2*b + 2*h)/ (g*self.thickness[i])) # GJ from the spar only
                self.GJ_Nodal.append(self.settings['GJ_everything'] + GJ_temp) # GJ total including the spar and the initial GJ from all of the other materials
 
        # GJ elemental
        self.GJ = []
        for i in range(0, len(self.GJ_Nodal) - 2, 2):
            self.GJ.append(abs((self.GJ_Nodal[i] + self.GJ_Nodal[i + 1] + self.GJ_Nodal[i + 2]) / 3))


        stiff_temp = []
        elem_stiff_temp = []
        for i in range(len(elem_stiff)):
            stiff_matrix = stiff[int(elem_stiff[i])]
            stiff_temp.append(stiff_matrix)
            elem_stiff_temp.append(int(i))
            if i == len(elem_stiff) - 1:
                stiff_temp = np.array(stiff_temp)

        for i in range(len(elements)):
            stiff_temp[int(elements[i])][3][3] = self.GJ[i]
            stiff_temp[int(elements[i])][4][4] = self.EI_Opt[i]
            last = int(elements[i])
            
        if self.settings['Sym'] == True:
            for i in range(len(elements)):
                index = last + 1 + i
                stiff_temp[index][3][3] = self.GJ[i]
                stiff_temp[index][4][4] = self.EI_Opt[i]

        # Updating stiffness values
        self.stiffness_db = np.array(stiff_temp)
        self.elem_stiffness = np.array(elem_stiff_temp)

        # Updating mass values
        self.mass_update()

        # strain from builtin strain calculation for comparison
        strain = beamloads.BeamLoads()
        strain.initialise(self.data)
        strain.calculate_loads(self)
        self.strain_compare = self.data.structure.timestep_info[-1].postproc_cell['strain']
        self.strain_elem = []

        for i in range(0, len(self.strain)-1, 2):
            self.strain_elem.append((self.strain[i] + self.strain[i+1] + self.strain[i+2])/3)

        self.strain_elem = np.array((self.strain_elem))

        # convergance criteria
        for i in range(len(self.EI_Opt)):
            if self.EI_Opt_Old == None:
                break
            else:
                if abs(self.EI_Opt[i] - self.EI_Opt_Old[i]) <= .01 * self.EI_Opt[i]:
                    if i == len(self.EI_Opt)-1:
                        self.Pass = True
                else:
                    break
    
        if self.iter == self.settings['Max_iter']:
            print('Max Iter')
            # print(self.strain_elem)
            # print(self.strain_compare)
            print(self.EI_Opt)

            # fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(18, 5))
            # y = np.array(self.data.structure.timestep_info[-1].pos[:, 1])
            # test = np.linspace(0, 16, 17)
            # ax[0, 0].plot(y, self.Cap_Area, 'ro')
            # ax[0, 0].plot(y, self.cap_area_original, 'bo')
            # ax[0, 1].plot(y, self.Cap_Area, 'ro')
            # ax[1, 0].plot(y, self.mass_nodal, 'ro')
            # ax[1, 0].plot(y, self.mass_OG, 'bo')
            # ax[1, 1].plot(test, self.mass_element, 'ro')
            # ax[1, 1].plot(test, self.mass_elem_OG, 'bo')
            # ax[1,0].plot(y, self.EI_Opt)
            # plt.show()

            self.Pass = True
            
        self.EI_Opt_Old = self.EI_Opt
        self.iter += 1
    
    def mass_update(self):
        
        mass_nodal = []
        #self.Cap_Area = self.cap_area_original # COMMENT OUT WHEN DONE
        elements = self.settings['Beam_Opt_Elems']
        nodes = self.settings['Beam_Opt_Nodes']
        
        # getting starting points of each node for beam opt
        start = []
        for i in range(len(elements)):
            if i == 0:
                start.append(self.data.structure.elements[int(elements[i])].coordinates_ini[0, 1])  
                start.append(self.data.structure.elements[int(elements[i])].coordinates_ini[2, 1])  
                start.append(self.data.structure.elements[int(elements[i])].coordinates_ini[1, 1])  
            else:
                start.append(self.data.structure.elements[int(elements[i])].coordinates_ini[2, 1])  
                start.append(self.data.structure.elements[int(elements[i])].coordinates_ini[1, 1])  
                
        for i in range(len(self.area)):
                mass_nodal.append(self.settings['Mass_Everything'] + self.settings['density_spar']* self.area[i])

        self.mass_nodal = np.array(mass_nodal)

        mass_elem = []
        for i in range(0, len(mass_nodal)-1, 2):
            mass_elem.append((mass_nodal[i] + mass_nodal[i+1] + mass_nodal[i+2])/3)
        self.mass_element = np.array(mass_elem)

        mass = self.data.structure.mass_db
        elem_mass = self.data.structure.elem_mass
        mass_temp = []
        elem_mass_temp = []

        for i in range(len(elem_mass)):
            mass_matrix = mass[int(elem_mass[i])]
            mass_temp.append(mass_matrix)
            elem_mass_temp.append(int(i))

        for i in range(len(elements)):
            mass_temp[int(elements[i])][0][0] = mass_elem[i]
            mass_temp[int(elements[i])][1][1] = mass_elem[i]
            mass_temp[int(elements[i])][2][2] = mass_elem[i]
            last = int(elements[i])

        if self.settings['Sym'] == True:
            for i in range(len(elements)):
                index = last + 1 + i
                mass_temp[index][0][0] = mass_elem[i]
                mass_temp[index][1][1] = mass_elem[i]
                mass_temp[index][2][2] = mass_elem[i]

        # Updating Mass Matrix
        self.elem_mass = np.array(elem_mass_temp)
        self.mass_db = np.array(mass_temp)


    def generate_fem(self):
        dat = self.data.structure
        coordinates = np.zeros((self.data.structure.num_node, 3))
        count = 0
        for i in range(len(self.data.structure.elements)):
            for j in range(len(self.data.structure.elements[0].coordinates_ini)):
                if count == 0 or count == 1 or  count == 2:
                    if j == 0:
                        coordinates[count] = self.data.structure.elements[i].coordinates_ini[j, :]
                        count += 1
                    elif j ==1:
                        coordinates[count] = self.data.structure.elements[i].coordinates_ini[j+1, :]
                        count += 1
                    elif j ==2:
                        coordinates[count] = self.data.structure.elements[i].coordinates_ini[j-1, :]
                        count += 1
                else:
                    if j == 1:
                        coordinates[count] = self.data.structure.elements[i].coordinates_ini[j + 1, :]
                        count += 1
                    elif j == 2:
                        coordinates[count] = self.data.structure.elements[i].coordinates_ini[j - 1, :]
                        count += 1

        in_data = {
            'app_forces': np.zeros((dat.num_node, 6)),
            'beam_number': dat.beam_number,
            'boundary_conditions': dat.boundary_conditions,
            'connectivities': dat.connectivities,
            'coordinates': coordinates,
            'elem_mass': self.elem_mass,
            'elem_stiffness': self.elem_stiffness,
            'frame_of_reference_delta': dat.frame_of_reference_delta,
            'lumped_mass': dat.lumped_mass,
            'lumped_mass_inertia': dat.lumped_mass_inertia,
            'lumped_mass_nodes': dat.lumped_mass_nodes,
            'lumped_mass_position': dat.lumped_mass_position,
            'mass_db': self.mass_db,
            'num_elem': dat.num_elem,
            'num_node': dat.num_node,
            'num_node_elem': dat.num_node_elem,
            'stiffness_db': self.stiffness_db,
            'structural_twist': dat.structural_twist,
        }

        settings = self.data.structure.settings
        beams = beam.Beam()
        self.data.structure = beams
        beams.generate(in_data, settings)

    def checks(self):

        # checking Margin of Safety for torsion
        forces = self.Nodal_Forces
        nodes = self.settings['Beam_Opt_Nodes']
        M_S = []
        F_U = self.settings['Ultimate_Shear']

        # cross sectional area of the spar
        area = (self.settings['B']* self.settings['H'])
        
        for i in range(len(nodes)):
            # torsional shear
            shear_flow_torsion = forces[int(nodes[i]), 3] / (2*area) 
            #up_lwr_shear = shear_flow_torsion / self.up_low_thick[i]
            if self.thickness[i] == 0:
                shear_tors = 0
            else:
                shear_tors = shear_flow_torsion / self.thickness[i]

            # bending shear
            if i == 0:
                bending_shear = .5 * self.shear[i] / area
            else:
                bending_shear = self.shear[i] / area

            shear_total = shear_tors + bending_shear

            M_S.append((F_U / abs(shear_total)) - 1)

            if min(M_S) < 2:
                print('Torsion Margin of Safety ERROR')
                # sys.exit()
        pass