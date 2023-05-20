import numpy as np
from OSU_Contribution.utils.Rotation_Matrices import local_for_forces

# Original MS calculations with dynamic gust considerations due to multiple time steps
def ms_OG_DYN(gust_calc, spar_cap_dict, trim_set, spar_geom, chord_vals):
    if gust_calc.structure.div_flag.value != 1:
        nodes = trim_set['beam_opt_nodes']
        shear = []
        moment = []
        node_forces_total = []
        local_applied_forces_log_total = []
        for ii in range(gust_calc.ts + 1):
            node_forces, local_applied_forces_log = local_for_forces(gust_calc, ii) 
            
            node_forces_total.append(node_forces)
            local_applied_forces_log_total.append(local_applied_forces_log)
            
            shear_temp = np.zeros(len(nodes))
            for i in range(len(nodes)-1, -1, -1):
                if i == len(nodes)-1:
                    shear_temp[i] = node_forces[i][2]
                elif i == 0:
                    shear_temp[i] = .5 * node_forces[i][2] + shear_temp[i+1]
                else:
                    shear_temp[i] = node_forces[i][2] + shear_temp[i+1]
            shear.append(shear_temp)

            moment_temp =  np.zeros(len(nodes))
            for i in range(len(nodes) - 1, -1, -1):
                if i == len(nodes) - 1:
                    moment_temp[i] = 0
                else:
                    back = gust_calc.structure.timestep_info[ii].pos[i+1]
                    current = gust_calc.structure.timestep_info[ii].pos[i]
                    moment_temp[i] = shear_temp[i+1] * (back[1] - current[1]) + moment_temp[i+1]
            moment.append(moment_temp)

        strain = np.zeros((gust_calc.ts + 1, len(nodes))) 
        # elem_stiff = gust_calc.structure.elem_stiffness
        # stiff = gust_calc.structure.stiffness_db
        m_s_ftu = np.zeros((gust_calc.ts + 1, len(nodes))) 
        m_s_fcu = np.zeros((gust_calc.ts + 1, len(nodes))) 
        
        # cap structural calculations
        for i in range(len(strain)):
            elem_stiff_counter = int(0)
            for j in range(len(strain[0])):
                strain[i][j] = abs((moment[i][j] * spar_geom['h_spar'][j]) / (spar_geom['ei_span'][j])) 
                if strain[i][j] != 0:
                    m_s_ftu[i][j] = spar_cap_dict['F_tu']/(spar_cap_dict['E_cap']*strain[i][j]) - 1
                    m_s_fcu[i][j] = spar_cap_dict['F_cu']/(spar_cap_dict['E_cap']*strain[i][j]) - 1
                else:
                    m_s_ftu[i][j] = 1000
                    m_s_fcu[i][j] = 1000
                if ( j % 2) == 0 and j != 0:
                    elem_stiff_counter +=1
        
        # SHEAR MARGIN OF SAFETY
        
        m_s_su = np.zeros((gust_calc.ts + 1, len(nodes))) 
        
        elastic_axis = []
        
        # laminant shear stress calculations
        ea_temp = gust_calc.aero.aero_dict['elastic_axis']
        for i in range(len(trim_set['beam_opt_elems'])):
            if i == 0:
                elastic_axis.append(ea_temp[i][0])
                elastic_axis.append(ea_temp[i][2])
                elastic_axis.append(ea_temp[i][1])
            else:
                elastic_axis.append(ea_temp[i][2])
                elastic_axis.append(ea_temp[i][1])
                        
        torsion = np.zeros((gust_calc.ts + 1, len(nodes)))
        for i in range(len(moment)):
            for j in range(len(nodes)-1, -1, -1):
                if j == len(nodes)-1:
                    torsion[i][j] = local_applied_forces_log_total[i][j][2] * (chord_vals[j]*(elastic_axis[j]-.25))
                else:
                    torsion[i][j] = local_applied_forces_log_total[i][j][2] * (chord_vals[j]*(elastic_axis[j]-.25)) + torsion[i][j+1]
        
          
        for i in range(len(torsion)):
            for j in range(len(nodes)):
                # torsional shear    
                shear_flow_torsion = torsion[i][j] / (2*spar_geom['a_cross'][j])
                
                #up_lwr_shear = shear_flow_torsion / self.up_low_thick[i]
                fwd_shear_tors = shear_flow_torsion / spar_geom['t_fwdaft'][j]
                aft_shear_tors = -shear_flow_torsion / spar_geom['t_fwdaft'][j]

                # bending shear
                if j == 0:
                    bending_shear = .5 * shear[i][j] / spar_geom['a_fwdaft'][j]
                else:
                    bending_shear = shear[i][j] / spar_geom['a_fwdaft'][j]

                fwd_shear = fwd_shear_tors + bending_shear
                aft_shear = aft_shear_tors + bending_shear
            
                M_S_fwd = (spar_cap_dict['F_su'] / abs(fwd_shear)) - 1
                M_S_aft = (spar_cap_dict['F_su'] / abs(aft_shear)) - 1
                m_s_su[i][j] = min([M_S_fwd, M_S_aft])
        
        return m_s_fcu, m_s_ftu, m_s_su
            
        
    # if diverged, has inputs for required variables
    else:
        m_s_fcu = [[.05, .05], [.05, .05]]
        m_s_ftu = [[.05, .05], [.05, .05]]
        m_s_su = [[.05, .05], [.05, .05]]
        return m_s_fcu, m_s_ftu, m_s_su
    
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------    
def ms_OG_STATIC(static_calc, spar_cap_dict, trim_set, spar_geom, chord_vals):
    if static_calc.structure.div_flag.value != 1:
        node_forces, local_applied_forces_log = local_for_forces(static_calc, static_calc.ts)
        
        nodes = trim_set['beam_opt_nodes']

        shear = np.zeros(len(nodes))
        for i in range(len(nodes)-1, -1, -1):
            if i == len(nodes)-1:
                shear[i] = node_forces[i][2]
            elif i == 0:
                shear[i] = .5 * node_forces[i][2] + shear[i+1]
            else:
                shear[i] = node_forces[i][2] + shear[i+1]

        moment =  np.zeros(len(nodes))
        for i in range(len(nodes) - 1, -1, -1):
            if i == len(nodes) - 1:
                moment[i] = 0
            else:
                back = static_calc.structure.timestep_info[-1].pos[i+1]
                current = static_calc.structure.timestep_info[-1].pos[i]
                moment[i] = shear[i+1] * (back[1] - current[1]) + moment[i+1]


        strain = np.zeros(len(nodes)) 
        # elem_stiff = static_calc.structure.elem_stiffness
        # stiff = static_calc.structure.stiffness_db
        # elem_stiff_counter = int(0)
        m_s_ftu = [] 
        m_s_fcu = [] 

        m_s_tweb = []
        
        # cap structural calculations
        for i in range(len(nodes)):
            strain[i] = abs((moment[i] * spar_geom['h_spar'][i]) / (spar_geom['ei_span'][i])) # ei other taken out to find strain on spar
            if strain[i] != 0:
                m_s_ftu.append(spar_cap_dict['F_tu']/(spar_cap_dict['E_cap']*strain[i]) - 1)
                m_s_fcu.append(spar_cap_dict['F_cu']/(spar_cap_dict['E_cap']*strain[i]) - 1)
                m_s_tweb.append(spar_cap_dict['F_tu_web']/(spar_cap_dict['E_web']*strain[i]) - 1)
            else:
                m_s_ftu.append(1000)
                m_s_fcu.append(1000)
                m_s_tweb.append(1000)
            # if ( i % 2) == 0 and i != 0:
            #     elem_stiff_counter +=1
        
        # SHEAR MARGIN OF SAFETY
        
        m_s_su = []
        
        elastic_axis = []
        row_count = 0
        order_count = 0
        
        for i in range(len(nodes)):
            if i == 0:
                order = [0, 2, 1]
            else:
                order = [2, 1]    
            if len(order) == order_count:
                row_count += 1
                order_count = 0
            elastic_axis.append(static_calc.aero.aero_dict['elastic_axis'][row_count][order[order_count]])
            order_count += 1
            
        torsion = []
        for i in range(len(nodes)-1, -1, -1):
            if i == len(nodes)-1:
                torsion.append(local_applied_forces_log[i][2] * (chord_vals[i]*(elastic_axis[i]-.25)))
            else:
                torsion.append(local_applied_forces_log[i][2] * (chord_vals[i]*(elastic_axis[i]-.25)) + torsion[-1])
        
        torsion.reverse()
        for i in range(len(nodes)):
            # torsional shear    
            shear_flow_torsion = torsion[i] / (2*spar_geom['a_cross'][i])
            
            #up_lwr_shear = shear_flow_torsion / self.up_low_thick[i]
            fwd_shear_tors = shear_flow_torsion / spar_geom['t_fwdaft'][i]
            aft_shear_tors = -shear_flow_torsion / spar_geom['t_fwdaft'][i]

            # bending shear
            if i == 0:
                bending_shear = .5 * shear[i] / spar_geom['a_fwdaft'][i]
            else:
                bending_shear = shear[i] / spar_geom['a_fwdaft'][i]
            
            fwd_shear = fwd_shear_tors + bending_shear
            aft_shear = aft_shear_tors + bending_shear
           
            M_S_fwd = (spar_cap_dict['F_su'] / abs(fwd_shear)) - 1
            M_S_aft = (spar_cap_dict['F_su'] / abs(aft_shear)) - 1
            m_s_su.append(min([M_S_fwd, M_S_aft]))

        
        return m_s_fcu, m_s_ftu, m_s_su, m_s_tweb
            
        
    # if diverged, has inputs for required variables
    else:
        m_s_fcu = [.1]
        m_s_ftu = [.1]
        m_s_su = [.1]
        m_s_tweb = [.1]
        
        return m_s_fcu, m_s_ftu, m_s_su, m_s_tweb
    
    
def warren_check(static_calc, spar_cap_dict, trim_set, spar_geom, chord_vals):
    if static_calc.structure.div_flag.value != 1:
        node_forces, local_applied, local_grav, app_copy, grav_copy = local_for_forces(static_calc, static_calc.ts)
        
        nodes = trim_set['beam_opt_nodes']
        
        shear = np.zeros(len(nodes))
        for i in range(len(nodes)-1, -1, -1):
            if i == len(nodes)-1:
                shear[i] = node_forces[i][2]
            elif i == 0:
                shear[i] = .5 * node_forces[i][2] + shear[i+1]
            else:
                shear[i] = node_forces[i][2] + shear[i+1]

        moment =  np.zeros(len(nodes))
        for i in range(len(nodes) - 1, -1, -1):
            if i == len(nodes) - 1:
                moment[i] = 0
            else:
                back = static_calc.structure.timestep_info[-1].pos[i+1]
                current = static_calc.structure.timestep_info[-1].pos[i]
                moment[i] = shear[i+1] * (back[1] - current[1]) + moment[i+1]


        strain = np.zeros(len(nodes)) 
        # elem_stiff = static_calc.structure.elem_stiffness
        # stiff = static_calc.structure.stiffness_db
        # elem_stiff_counter = int(0)
        m_s_ftu = [] 
        m_s_fcu = [] 

        m_s_tweb = []
        
        # cap structural calculations
        for i in range(len(nodes)):
            strain[i] = abs((moment[i] * spar_geom['h_spar'][i]) / (spar_geom['ei_span'][i])) # ei other taken out to find strain on spar
            if strain[i] != 0:
                m_s_ftu.append(spar_cap_dict['F_tu']/(spar_cap_dict['E_cap']*strain[i]) - 1)
                m_s_fcu.append(spar_cap_dict['F_cu']/(spar_cap_dict['E_cap']*strain[i]) - 1)
                m_s_tweb.append(spar_cap_dict['F_tu_web']/(spar_cap_dict['E_web']*strain[i]) - 1)
            else:
                m_s_ftu.append(1000)
                m_s_fcu.append(1000)
                m_s_tweb.append(1000)
            # if ( i % 2) == 0 and i != 0:
            #     elem_stiff_counter +=1
        
        # SHEAR MARGIN OF SAFETY
        
        m_s_su = []
        
        elastic_axis = []
        row_count = 0
        order_count = 0
        
        for i in range(len(nodes)):
            if i == 0:
                order = [0, 2, 1]
            else:
                order = [2, 1]    
            if len(order) == order_count:
                row_count += 1
                order_count = 0
            elastic_axis.append(static_calc.aero.aero_dict['elastic_axis'][row_count][order[order_count]])
            order_count += 1
            
        torsion = []
        restor_moment_trans = []
        restor_moment_copy = []
        # m_ea
        moment_trans = []
        moment_copy = []
        # m_struct and m gondola
        m_grav_trans = []
        m_grav_copy = []
        for i in range(len(nodes)-1, -1, -1):
            torsion.append(node_forces[i][3])
            restor_moment_trans.append(local_applied[i][2] * (chord_vals[i]*(elastic_axis[i]-.25)))
            restor_moment_copy.append(app_copy[i][2] * (chord_vals[i]*(elastic_axis[i]-.25)))
            
            moment_trans.append(local_applied[i][3])
            moment_copy.append(app_copy [i][3])

            m_grav_trans.append(local_grav[i][3])
            m_grav_copy.append(grav_copy[i][3])
            # if i == len(nodes)-1:
            #     torsion.append(node_forces[i][3])
            # else:
            #     torsion.append(node_forces[i][3] + torsion[-1])
            # if i == len(nodes)-1:
            #     torsion.append(local_applied_forces_log[i][2] * (chord_vals[i]*(elastic_axis[i]-.25)))
            # else:
            #     torsion.append(local_applied_forces_log[i][2] * (chord_vals[i]*(elastic_axis[i]-.25)) + torsion[-1])
        
        torsion.reverse()
        for i in range(len(nodes)):
            # torsional shear    
            shear_flow_torsion = torsion[i] / (2*spar_geom['a_cross'][i])
            
            #up_lwr_shear = shear_flow_torsion / self.up_low_thick[i]
            fwd_shear_tors = shear_flow_torsion / spar_geom['t_fwdaft'][i]
            aft_shear_tors = -shear_flow_torsion / spar_geom['t_fwdaft'][i]

            # bending shear
            if i == 0:
                bending_shear = .5 * shear[i] / spar_geom['a_fwdaft'][i]
            else:
                bending_shear = shear[i] / spar_geom['a_fwdaft'][i]
            
            fwd_shear = fwd_shear_tors + bending_shear
            aft_shear = aft_shear_tors + bending_shear
           
            M_S_fwd = (spar_cap_dict['F_su'] / abs(fwd_shear)) - 1
            M_S_aft = (spar_cap_dict['F_su'] / abs(aft_shear)) - 1
            m_s_su.append(min([M_S_fwd, M_S_aft]))

        
        return m_s_fcu, m_s_ftu, m_s_su, m_s_tweb, shear, moment, torsion, restor_moment_trans, restor_moment_copy, moment_trans, moment_copy, m_grav_trans, m_grav_copy
            
        
    # if diverged, has inputs for required variables
    else:
        m_s_fcu = [.1]
        m_s_ftu = [.1]
        m_s_su = [.1]
        m_s_tweb = [.1]
        shear = 0
        moment = 0
        
        return m_s_fcu, m_s_ftu, m_s_su, m_s_tweb, shear, moment, torsion
        