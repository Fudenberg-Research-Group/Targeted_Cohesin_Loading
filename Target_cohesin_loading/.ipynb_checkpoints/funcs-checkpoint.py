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
                    LEF_ipause,
                    sites_per_monomer,
                    velocity_multiplier,
                    **kwargs):
    
    lifetime_array = make_site_array(site_types, LEF_lifetime, **kwargs)
    stalled_lifetime_array = make_site_array(site_types, LEF_stalled_lifetime, **kwargs)
    
    birth_array = make_site_array(site_types, LEF_birth, **kwargs)
    pause_array = make_site_array(site_types, LEF_pause, **kwargs)
    ipause_array = make_site_array(site_types, LEF_ipause, **kwargs)
    
    death_array = 1./ lifetime_array / (velocity_multiplier * sites_per_monomer)
    stalled_death_array = 1./ stalled_lifetime_array / (velocity_multiplier * sites_per_monomer)

    return [death_array, stalled_death_array, birth_array, pause_array, ipause_array]

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

    LEFTran = extrusion_engine(number_of_LEFs, *LEF_arrays, *CTCF_arrays, *CTCF_dynamic_arrays)

    if not isinstance(LEFTran, LEFTranslocatorDynamicBoundary):
        LEFTran.stallProbLeft = 1 - (1 - LEFTran.stallProbLeft) ** (1. / velocity_multiplier)
        LEFTran.stallProbRight = 1 - (1 - LEFTran.stallProbRight) ** (1. / velocity_multiplier)

    return LEFTran

def paramdict_to_filename(paramdict, paramdict_keys=None):
    if paramdict_keys is None:
        paramdict_keys = {
            'CTCF_facestall':'face',
            'CTCF_backstall':'back',
            'CTCF_lifetime':'clife',
            'CTCF_offtime':'cof',
            'CTCF_number':'cnum',
            'LEF_lifetime':'life',
            'LEF_stalled_lifetime':'slife',
            'LEF_birth':'birth',
            'targetsnum':'targetsnum',
            'deltactcf':'deltactcf',
            'LEF_ipause':'ipause',
            'LEF_pause':'pause',
            'LEF_separation':'sep',
            'sites_per_monomer':'site',
            'monomers_per_replica':'monomer',
            'number_of_replica':'replica',
            'steps':'steps',
            'velocity_multiplier':'vel'
        }

    filename = 'file'
    for key, value in paramdict.items():
        short_key = paramdict_keys.get(key, key)  # fallback if not found
        if isinstance(value, list):
            value_str = '_'.join(str(v) for v in value)
        else:
            value_str = str(value)
        filename += f'_{short_key}_{value_str}'

    # Remove characters that might interfere with file naming
    filename = filename.replace('[', '').replace(']', '').replace(' ', '')
    return filename

def adjust_LEF_density(paramdict, base_loading=0.0001):
    """
    Computes adjusted LEF_separation based on new loading on target.
    Returns an integer value for the updated LEF_separation.
    """
    birth = paramdict['LEF_birth'][1] if isinstance(paramdict['LEF_birth'], list) else paramdict['LEF_birth']
    monomers_per_replica = int(paramdict['monomers_per_replica'])
    sites_per_monomer = int(paramdict['sites_per_monomer'])
    number_of_replica = int(paramdict['number_of_replica'])

    total_sites = monomers_per_replica * sites_per_monomer * number_of_replica
    density_multiplier = 1 + (((birth - base_loading) / base_loading) / total_sites)
    
    return int(paramdict['LEF_separation'] / density_multiplier)


