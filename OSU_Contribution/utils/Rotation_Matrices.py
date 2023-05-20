import numpy as np
import sharpy.utils.algebra as algebra

# converts from aircraft frame of reference to local frame of reference
def nodal_a_for_2_b_for(data, nodal, tstep, filter=np.array([True] * 6)):
    nodal_a = nodal.copy(order='F')
    for i_node in range(data.structure.num_node):
        # get master elem and i_local_node
        i_master_elem, i_local_node = data.structure.node_master_elem[i_node, :]
        crv = tstep.psi[i_master_elem, i_local_node, :]
        cab = algebra.crv2rotation(crv)
        temp = np.zeros((6,))
        temp[0:3] = np.dot(cab, nodal[i_node, 0:3])
        temp[3:6] = np.dot(cab, nodal[i_node, 3:6])
        for i in range(6):
            if filter[i]:
                nodal_a[i_node, i] = temp[i]
    return nodal_a

def local_for_forces(data, time):
    if time is None:
        tstep = data.structure.timestep_info[data.ts]
    else:
        tstep = data.structure.timestep_info[time]
    applied_forces = data.structure.nodal_b_for_2_a_for(tstep.steady_applied_forces,
                                                                tstep)

    applied_forces_copy = applied_forces.copy()
    global lift
    lift = applied_forces_copy
    #lift[0] = lift[0]*2 
    gravity_forces_copy = tstep.gravity_forces.copy()
    
    
    for i_node in range(data.structure.num_node):
        applied_forces_copy[i_node, 3:6] += np.cross(tstep.pos[i_node, :],
                                                        applied_forces_copy[i_node, 0:3])
        gravity_forces_copy[i_node, 3:6] += np.cross(tstep.pos[i_node, :],
                                                        gravity_forces_copy[i_node, 0:3])


    # forces are in frame of reference a
    together = np.zeros((len(applied_forces_copy), 6))
    for i in range(len(applied_forces_copy)):
        if i == 0:
            together[i, ...] = applied_forces_copy[i] + gravity_forces_copy[i]
        else:
            together[i, ...] = applied_forces_copy[i] + gravity_forces_copy[i]
            
    # # only need half at the beginning of span
    # applied_forces_copy[0] = applied_forces_copy[0] * .5
    # gravity_forces_copy[0] = gravity_forces_copy[0] * .5

    # for the logger 
    local_grav_forces_log = nodal_a_for_2_b_for(data, gravity_forces_copy, tstep)
    local_applied_forces_log = nodal_a_for_2_b_for(data, applied_forces_copy, tstep)
    
    return nodal_a_for_2_b_for(data, together, tstep), local_applied_forces_log


def local_for_forces_warren(data, time):
    if time is None:
        tstep = data.structure.timestep_info[data.ts]
    else:
        tstep = data.structure.timestep_info[time]
    applied_forces = data.structure.nodal_b_for_2_a_for(tstep.steady_applied_forces,
                                                                tstep)

    applied_forces_copy = applied_forces.copy()
    global lift
    lift = applied_forces_copy
    #lift[0] = lift[0]*2 
    gravity_forces_copy = tstep.gravity_forces.copy()
    
    
    for i_node in range(data.structure.num_node):
        applied_forces_copy[i_node, 3:6] += np.cross(tstep.pos[i_node, :],
                                                        applied_forces_copy[i_node, 0:3])
        gravity_forces_copy[i_node, 3:6] += np.cross(tstep.pos[i_node, :],
                                                        gravity_forces_copy[i_node, 0:3])


    # forces are in frame of reference a
    together = np.zeros((len(applied_forces_copy), 6))
    for i in range(len(applied_forces_copy)):
        if i == 0:
            together[i, ...] = applied_forces_copy[i] + gravity_forces_copy[i]
        else:
            together[i, ...] = applied_forces_copy[i] + gravity_forces_copy[i]
            
    # # only need half at the beginning of span
    # applied_forces_copy[0] = applied_forces_copy[0] * .5
    # gravity_forces_copy[0] = gravity_forces_copy[0] * .5

    # for the logger 
    local_grav_forces_log = nodal_a_for_2_b_for(data, gravity_forces_copy, tstep)
    local_applied_forces_log = nodal_a_for_2_b_for(data, applied_forces_copy, tstep)
    
    return nodal_a_for_2_b_for(data, together, tstep), local_applied_forces_log, local_grav_forces_log, applied_forces_copy, gravity_forces_copy