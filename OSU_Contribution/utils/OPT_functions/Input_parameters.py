
def scaled_SI2_inputs():

    spar_cap_dict = {}
    spar_cap_dict['E_cap'] = 1.783e11
    spar_cap_dict['E_lam'] = 1.212e10
    spar_cap_dict['G'] = 4.51e10
    spar_cap_dict['Misc_Factor'] = .59

    spar_cap_dict['x_taper_start'] = span * 0.686 / 2

    # height and width of spar before taper
    spar_cap_dict['h_spar_1'] = 0.388
    spar_cap_dict['b_spar_1'] = 0.703
    #height and widrh of spar after taper
    spar_cap_dict['h_spar_2'] = 0.357
    spar_cap_dict['b_spar_2'] = 0.477

    spar_cap_dict['density'] = 1540 

    # margin of safety calculations
    spar_cap_dict['F_tu'] = 1.283e9
    spar_cap_dict['F_cu'] = 4.61e8
    spar_cap_dict['F_su'] = 2.35e8

    return spar_cap_dict


def T_Tail_inputs():

    spar_cap_dict = {}
    spar_cap_dict['E_cap'] = 1.9e11
    spar_cap_dict['E_lam'] = 1.3e10
    spar_cap_dict['G'] = 4.8e10

    spar_cap_dict['h_spar'] = 0.129
    spar_cap_dict['b_spar'] = 0.2

    spar_cap_dict['density'] = 1540 #163.934

    spar_cap_dict['mass_other'] = .5 # ratio of the mass that comes from sources other than the spar 
    # margin of safety calculations
    spar_cap_dict['F_tu'] = 1.35e9
    spar_cap_dict['F_cu'] = 4.85e8
    spar_cap_dict['F_su'] = 2.5e8

    return spar_cap_dict