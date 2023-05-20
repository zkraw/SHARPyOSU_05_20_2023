
def scaled_SI2_OG(beam_opt_nodes, fem, spar_cap_dict, f, chord_vals): 
    
    max_span =  fem['coordinates'][beam_opt_nodes[-1]][1] 
    
    thickness = []
    span_loc_log = []
    h_spar = []
    b_spar = []
    # finding ei along the span
    x_tap = spar_cap_dict['x_taper_start']
    h_1 = spar_cap_dict['h_spar_1']
    h_2 = spar_cap_dict['h_spar_2'] 
    b_1 = spar_cap_dict['b_spar_1'] 
    b_2 = spar_cap_dict['b_spar_2'] 
    
    for i in range(len(beam_opt_nodes)):
        span_loc_log.append(fem['coordinates'][i][1])
        thickness.append(f(fem['coordinates'][i][1]/max_span))
        if fem['coordinates'][i][1] <= x_tap:
            h_spar.append(h_1)
            b_spar.append(b_1)
        else:
            h_spar.append(((h_2 - h_1)*(fem['coordinates'][i][1] - x_tap)/(max_span - x_tap) + h_1))
            b_spar.append(((b_2 - b_1)*(fem['coordinates'][i][1] - x_tap)/(max_span - x_tap) + b_1))
        
    
    # stiffness and mass value updates
    ei_span = []
    gj_span = []
    m_span = []
    # torsion ms values
    a_cross = []
    a_fwdaft = []
    a_topbot = []
    t_fwdaft = []
    H_couple = []
    # constant values of the geometry 
    for i in range(len(thickness)):
        thick = thickness[i]
        # new width values 
        b_cap = .15*b_spar[i]
        b_lam_up_low = b_spar[i] - 2 * b_cap
       
        # EI Calculation
        H_couple.append(h_spar[i]/2)
        
        #topbot sections
        a_topbot.append(thick * b_lam_up_low)
        EI_topbot = spar_cap_dict['E_lam'] * a_topbot[-1] * H_couple[-1]**2
        
        #fwd/aft sections
        t_fwdaft.append(thick)
        h_fwdaft = h_spar[i] - 2 * thick
        a_fwdaft.append(t_fwdaft[-1]*h_fwdaft)
        EI_fwdaft = spar_cap_dict['E_lam'] * t_fwdaft[-1] * h_fwdaft ** 3 / 12
        
        #cap sections
        a_cap = thick * b_cap
        EI_cap = spar_cap_dict['E_cap'] * a_cap * H_couple[-1]**2
        # total EI by section
        ei_span.append(2 * EI_topbot + 2 * EI_fwdaft + 4 * EI_cap)
        
        # GJ calculation
        s = 2 * (h_spar[i] - 2*thick) + 2 * (b_lam_up_low)
        a_cross.append(h_spar[i] * b_spar[i])
        gj_span.append(4 * a_cross[-1] ** 2 / (s / (spar_cap_dict['G'] * thick)))
        
        # spar mass contribution calculation
        m_span.append((spar_cap_dict['density']*(4*a_cap + 2*a_topbot[-1] + 2*a_fwdaft[-1])))
 
    ei_final = []
    #taking the average of the nodal ei
    for i in range(0, len(ei_span) - 2, 2):
        ei_final.append(abs((ei_span[i] + ei_span[i + 1] + ei_span[i + 2]) / 3)) 

    gj_final = []
    #taking the average of the nodal gj
    for i in range(0, len(gj_span) - 2, 2):
        gj_final.append(abs((gj_span[i] + gj_span[i + 1] + gj_span[i + 2]) / 3))
        
    
    # MASS CONTRIBUTION
    mass_final = []
    for i in range(0, len(m_span)-2, 2):
        m_solar = 0.437 * chord_vals[i]
        m_lete = .1043 * chord_vals[i]**2 * 3.722 
        mass_final.append(((m_span[i] + m_span[i+1] + m_span[i+2])/3)*(1+spar_cap_dict['Misc_Factor']) + m_solar + m_lete)
    
    # goemtry needed in margin of safety calculations
    spar_geom = {}
    spar_geom['a_cross'] = a_cross
    spar_geom['t_fwdaft'] = t_fwdaft
    spar_geom['a_fwdaft'] = a_fwdaft
    spar_geom['ei_span'] = ei_span 

    for i in range(len(h_spar)):
        h_spar[i] = h_spar*.5
    spar_geom['h_spar'] = h_spar

    return ei_final, gj_final, mass_final, spar_geom

#-------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------

def t_tail_og(beam_opt_nodes, fem, spar_cap_dict, f, chord_vals, scale): 
    
    max_span =  fem['coordinates'][beam_opt_nodes[-1]][1] 
    
    thickness = []
    span_loc_log = []
    h_spar = []
    b_spar = []
    # finding ei along the span
    h_1 = (1/scale)*spar_cap_dict['h_spar']
    b_1 = (1/scale)*spar_cap_dict['b_spar'] 

    
    for i in range(len(beam_opt_nodes)):
        span_loc_log.append(fem['coordinates'][i][1])
        # f() is the linear fit of thickness along the span found by the control point inputs
        thickness.append(f(fem['coordinates'][i][1]/max_span))
        h_spar.append(h_1)
        b_spar.append(b_1)
        
    
    # stiffness and mass value updates
    ei_span = []
    gj_span = []
    m_span = []
    # torsion ms values
    a_cross = []
    a_fwdaft = []
    a_topbot = []
    t_fwdaft = []
    H_couple = []
    # constant values of the geometry 
    for i in range(len(thickness)):
        thick = thickness[i]
        # new width values 
        b_cap = .15*b_spar[i]
        b_lam_up_low = b_spar[i] - 2 * b_cap
       
        # EI Calculation
        H_couple.append(h_spar[i]/2 - thick/2)
        
        #topbot sections
        a_topbot.append(thick * b_lam_up_low)
        EI_topbot = spar_cap_dict['E_lam'] * a_topbot[-1] * H_couple[-1]**2
        
        #fwd/aft sections
        t_fwdaft.append(thick)
        h_fwdaft = h_spar[i] - 2 * thick
        a_fwdaft.append(t_fwdaft[-1]*h_fwdaft)
        EI_fwdaft = spar_cap_dict['E_lam'] * t_fwdaft[-1] * h_fwdaft ** 3 / 12
        
        #cap sections
        a_cap = thick * b_cap
        EI_cap = spar_cap_dict['E_cap'] * a_cap * H_couple[-1]**2
        # total EI by section
        ei_span.append(2 * EI_topbot + 2 * EI_fwdaft + 4 * EI_cap)
        
        # GJ calculation
        s = 2 * (h_spar[i] - 2*thick) + 2 * (b_lam_up_low)
        a_cross.append(h_spar[i] * b_spar[i])
        gj_span.append(4 * a_cross[-1] ** 2 / (s / (spar_cap_dict['G'] * thick)))
        
        # spar mass contribution calculation
        m_span.append((spar_cap_dict['density']*(4*a_cap + 2*a_topbot[-1] + 2*a_fwdaft[-1])))
 
    ei_final = []
    #taking the average of the nodal ei
    for i in range(0, len(ei_span) - 2, 2):
        ei_final.append(abs((ei_span[i] + ei_span[i + 1] + ei_span[i + 2]) / 3)) 

    gj_final = []
    #taking the average of the nodal gj
    for i in range(0, len(gj_span) - 2, 2):
        gj_final.append(abs((gj_span[i] + gj_span[i + 1] + gj_span[i + 2]) / 3))
        
    
    # MASS CONTRIBUTION 
    mass_final = []
    # mass from solar panels and other components of wing
    m_solar = .75 * chord_vals[0]
    m_lete = .10226 * chord_vals[0]**2 * 3.5 
    for i in range(0, len(m_span)-2, 2):
        mass_final.append((m_span[i] + m_span[i+1] + m_span[i+2])/3 + m_solar + m_lete)
    
    # goemtry needed in margin of safety calculations
    spar_geom = {}
    spar_geom['a_cross'] = a_cross
    spar_geom['t_fwdaft'] = t_fwdaft
    spar_geom['a_fwdaft'] = a_fwdaft
    spar_geom['ei_span'] = ei_span 
    
    # divide by half so that y_max from centroid is known for strain calculation
    for i in range(len(h_spar)):
        h_spar[i] = h_spar*.5
    spar_geom['h_spar'] = h_spar

    return ei_final, gj_final, mass_final, spar_geom 

#-------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------

def t_tail_three_variable(beam_opt_nodes, fem, spar_cap_dict, f, chord_vals, scale, max_height):
    max_span =  fem['coordinates'][beam_opt_nodes[-1]][1] 
    
    span_loc_log = []
    b_cap = []
    h_spar = []
    b_spar = []
    
    for i in range(len(beam_opt_nodes)):
        span_loc_log.append(fem['coordinates'][i][1])
        # f() is the linear fit of thickness along the span found by the control point inputs
        b_cap.append(f[0](fem['coordinates'][i][1]/max_span))
        b_spar.append(f[1](fem['coordinates'][i][1]/max_span))
        h_spar.append(f[2](fem['coordinates'][i][1]/max_span))
        # h_spar.append(max_height)
        
    # stiffness and mass value updates
    ei_span = []
    gj_span = []
    m_span = []
    # torsion ms values
    a_cross = []
    a_fwdaft = []
    a_topbot = []
    H_couple = []
    t_fwdaft = []
    core_area_top = []
    core_area_aft = []
    # constant values of the geometry 
    for i in range(len(beam_opt_nodes)):
        thick = spar_cap_dict['thickness']
        core_laminate_thick = spar_cap_dict['thickness'] + spar_cap_dict['core']
        # EI Calculation
        H_couple.append(h_spar[i] - core_laminate_thick/2)
        b_lam_up_low = 2*(b_spar[i] - b_cap[i])
        
        #topbot sections
        a_topbot.append(thick * b_lam_up_low)
        core_area_top.append(spar_cap_dict['core']*b_lam_up_low)
        EI_topbot = spar_cap_dict['E_web'] * a_topbot[-1] * H_couple[-1]**2
        
        #fwd/aft sections
        t_fwdaft.append(thick)
        h_fwdaft = 2*(h_spar[i] - core_laminate_thick)
        a_fwdaft.append(thick*h_fwdaft)
        core_area_aft.append(spar_cap_dict['core']*h_fwdaft)
        EI_fwdaft = (spar_cap_dict['E_web'] * thick * h_fwdaft ** 3) / 12
        
        #cap sections
        a_cap = core_laminate_thick * b_cap[i]
        EI_cap = spar_cap_dict['E_cap'] * a_cap * H_couple[-1]**2
        # total EI by section
        ei_span.append(2 * EI_topbot + 2 * EI_fwdaft + 4 * EI_cap)
        
        # GJ calculation
        s = 2 * h_fwdaft + 2 * b_lam_up_low
        a_cross.append(2*h_spar[i] * 2*b_spar[i])
        gj_span.append((4 * a_cross[-1] ** 2) / (s / (spar_cap_dict['G'] * thick)))
        
        # spar mass contribution calculation
        m_span.append((spar_cap_dict['density']*(4*a_cap + 2*a_topbot[-1] + 2*a_fwdaft[-1])) + spar_cap_dict['core_density']*(2*core_area_aft[-1] + 2*core_area_top[-1]))
 
    ei_final = []
    #taking the average of the nodal ei
    for i in range(0, len(ei_span) - 2, 2):
        ei_final.append(abs((ei_span[i] + ei_span[i + 1] + ei_span[i + 2]) / 3)) 

    gj_final = []
    #taking the average of the nodal gj
    for i in range(0, len(gj_span) - 2, 2):
        gj_final.append(abs((gj_span[i] + gj_span[i + 1] + gj_span[i + 2]) / 3))
        
    
    # MASS CONTRIBUTION 
    mass_final = []
    # mass from solar panels and other components of wing
    m_solar = .75*chord_vals[0]
    m_lete = .10226*chord_vals[0]**2 * 3.5 
    for i in range(0, len(m_span)-2, 2):
        mass_final.append((m_span[i] + m_span[i+1] + m_span[i+2])/3 + m_solar + m_lete)
    
    # goemtry needed in margin of safety calculations
    spar_geom = {}
    spar_geom['a_cross'] = a_cross
    spar_geom['t_fwdaft'] = t_fwdaft
    spar_geom['a_fwdaft'] = a_fwdaft
    spar_geom['ei_span'] = ei_span 
    spar_geom['h_spar'] = h_spar

    return ei_final, gj_final, mass_final, spar_geom, m_span

 
def warren_interp(beam_opt_nodes, fem, spar_cap_dict, f, chord_vals, scale):
    span_loc_log = []
    ei_span = []
    gj_span = []
    
    a_cross = []
    t_fwdaft = []
    a_fwdaft = []
    t_fwdaft = []
    t_fwdaft = []
    h_spar = [] 
    
    for i in range(len(beam_opt_nodes)):
        span_loc_log.append(fem['coordinates'][i][1])
        # f() is the linear fit of thickness along the span found by the control point inputs
        ei_span.append(f[0](fem['coordinates'][i][1]))
        gj_span.append(f[1](fem['coordinates'][i][1]))

        a_cross.append(f[2](fem['coordinates'][i][1]))
        t_fwdaft.append(f[3](fem['coordinates'][i][1]))
        a_fwdaft.append(f[4](fem['coordinates'][i][1]))
        h_spar.append(f[5](fem['coordinates'][i][1]))
    
    # goemtry needed in margin of safety calculations
    spar_geom = {}
    spar_geom['a_cross'] = a_cross
    spar_geom['t_fwdaft'] = t_fwdaft
    spar_geom['a_fwdaft'] = a_fwdaft
    spar_geom['ei_span'] = ei_span 
    spar_geom['h_spar'] = h_spar


    ei_final = []
    #taking the average of the nodal ei
    for i in range(0, len(ei_span) - 2, 2):
        ei_final.append(abs((ei_span[i] + ei_span[i + 1] + ei_span[i + 2]) / 3)) 

    gj_final = []
    #taking the average of the nodal gj
    for i in range(0, len(gj_span) - 2, 2):
        gj_final.append(abs((gj_span[i] + gj_span[i + 1] + gj_span[i + 2]) / 3))

    return ei_span, gj_span, spar_geom 