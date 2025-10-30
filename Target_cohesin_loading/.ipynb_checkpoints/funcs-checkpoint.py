from lattice_translocators import LEFTranslocator, LEFTranslocatorDynamicBoundary

import numpy as np

def make_site_array(site_types, 
                    values, 
                    at_ids=None, 
                    number_of_replica=1, 
                    **kwargs):
    
    assert site_types.max() < len(values), ('Number of values (%d) incompatible with number of site types (%d)'
                                            % (len(values), site_types.max()))
    
    prop_array = np.zeros(len(site_types), dtype=np.double)
    
    for i, value in enumerate(values):
        prop_array[site_types == i] = value
        
    if isinstance(at_ids, np.ndarray):
        mask = np.zeros(len(site_types), dtype=bool)
        mask[at_ids] = True
        
        prop_array[~mask] = 0
        
    return np.tile(prop_array, number_of_replica)



def make_CTCF_arrays(site_types,
                     CTCF_left_positions,
                     CTCF_right_positions,
                     CTCF_facestall,
                     CTCF_backstall,
                     **kwargs):
    
    stall_left_array = make_site_array(site_types, CTCF_facestall, at_ids=CTCF_left_positions, **kwargs)
    stall_right_array = make_site_array(site_types, CTCF_facestall, at_ids=CTCF_right_positions, **kwargs)
    
    stall_left_array += make_site_array(site_types, CTCF_backstall, at_ids=CTCF_right_positions, **kwargs)
    stall_right_array += make_site_array(site_types, CTCF_backstall, at_ids=CTCF_left_positions, **kwargs)
    
    return [stall_left_array, stall_right_array]


def make_CTCF_dynamic_arrays(site_types,
                             CTCF_lifetime,
                             CTCF_offtime,
                             sites_per_monomer,
                             velocity_multiplier,
                             **kwargs):
    
    CTCF_lifetime_array = make_site_array(site_types, CTCF_lifetime, **kwargs)
    CTCF_offtime_array = make_site_array(site_types, CTCF_offtime, **kwargs)
    
    CTCF_death_array = 1./ CTCF_lifetime_array / (velocity_multiplier * sites_per_monomer)
    CTCF_birth_array = 1./ CTCF_offtime_array / (velocity_multiplier * sites_per_monomer)

    return [CTCF_death_array, CTCF_birth_array]


def make_LEF_arrays(site_types,
                    LEF_lifetime,
                    LEF_stalled_lifetime,
                    LEF_birth,
                    LEF_pause,
                    sites_per_monomer,
                    velocity_multiplier,
                    **kwargs):
    
    lifetime_array = make_site_array(site_types, LEF_lifetime, **kwargs)
    stalled_lifetime_array = make_site_array(site_types, LEF_stalled_lifetime, **kwargs)
    
    birth_array = make_site_array(site_types, LEF_birth, **kwargs)
    pause_array = make_site_array(site_types, LEF_pause, **kwargs)
    
    death_array = 1./ lifetime_array / (velocity_multiplier * sites_per_monomer)
    stalled_death_array = 1./ stalled_lifetime_array / (velocity_multiplier * sites_per_monomer)

    return [death_array, stalled_death_array, birth_array, pause_array]

def make_translocator(extrusion_engine, 
                      site_types,
                      CTCF_left_positions,
                      CTCF_right_positions,
                      **kwargs):

    LEF_separation = kwargs['LEF_separation']    
    velocity_multiplier = kwargs['velocity_multiplier'] 
    
    sites_per_monomer = kwargs['sites_per_monomer'] 
    
    number_of_replica = kwargs['number_of_replica'] 
    monomers_per_replica = kwargs['monomers_per_replica'] 

    number_of_monomers = number_of_replica * monomers_per_replica
    number_of_LEFs = number_of_monomers // LEF_separation
    
    sites_per_replica = monomers_per_replica*sites_per_monomer

    assert len(site_types) == sites_per_replica, ("Site type array (%d) doesn't match replica lattice size (%d)"
                                                  % (len(site_types), sites_per_replica))

    # Create arrays
    LEF_arrays = make_LEF_arrays(site_types, **kwargs)
    
    CTCF_arrays = make_CTCF_arrays(site_types, CTCF_left_positions, CTCF_right_positions, **kwargs)
    CTCF_dynamic_arrays = make_CTCF_dynamic_arrays(site_types, **kwargs)
    # pick what to forward
    allowed_engine_kwargs = ['LEF_ipause', 'initalize_at_equilibrium_occupancy']
    optional_args_for_engine = {key: kwargs[key] for key in allowed_engine_kwargs if key in kwargs}
    LEFTran = extrusion_engine(number_of_LEFs, *LEF_arrays, *CTCF_arrays, *CTCF_dynamic_arrays, **optional_args_for_engine)

    if not isinstance(LEFTran, LEFTranslocatorDynamicBoundary):
        LEFTran.stallProbLeft = 1 - (1 - LEFTran.stallProbLeft) ** (1. / velocity_multiplier)
        LEFTran.stallProbRight = 1 - (1 - LEFTran.stallProbRight) ** (1. / velocity_multiplier)

    return LEFTran

def paramdict_to_filename(paramdictx,paramdict_keys):
    
    filename='file'
    for i in range(len(paramdictx)):
        filename += ('_'+paramdict_keys[list(paramdictx)[i][:]]+'_'+str(paramdictx[list(paramdictx)[i]]))
    chars = ['[',']']
    filename_new = ''.join(i for i in filename if not i in chars)    
    return filename_new

