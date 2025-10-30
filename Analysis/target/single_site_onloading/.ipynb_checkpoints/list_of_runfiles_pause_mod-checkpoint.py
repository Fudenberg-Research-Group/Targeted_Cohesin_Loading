import glob
import ast 
import datetime

LEF_pause= 0#0.9

pause_multiplier = 1/(1-LEF_pause)

lifetimes_b = [50, 66]
lifetimes = [pause_multiplier*life for life in lifetimes_b]

velocities = [1]

face_stalls = [1.0]
back_stalls = [0]
#ctcf_life_range = [  0.15, 0.3, 1.5, 5.0, 15, 30, 50, 150, 300, 500, 1500, 5000, 15000]
#ctcf_oftime_range = [0.15, 0.3, 1.5, 5.0, 15, 30, 50, 150, 300, 500, 1500, 5000, 15000]
CTCF_lifetimes_b =  [ 17.0]

#CTCF_lifetimes_b =  [3.0]
CTCF_lifetimes = [pause_multiplier*clife for clife in CTCF_lifetimes_b]
#CTCF_offtimes_b = [  1.7, 17.0, 170]
CTCF_offtimes_b = [17.0]
#CTCF_offtimes_b = [0.3]
CTCF_offtimes = [pause_multiplier*cof for cof in CTCF_offtimes_b]

stall_dists = [100]
#LEF_births = [0.0001,  0.00025, 0.0005,  0.00075,  0.001, 0.0025,  0.005,  0.0075,  0.01, 0.025, 0.05,  0.075, 0.1]
LEF_births = [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5]
#delta_list = [ 4, 8, 10, 20, 40,  80, 100,  200, 400,  800, 1000, 1500, 3000, 5000] 
#delta_list = [4, 8, 10, 20, 40, 100, 200, 400, 1000, 1500, 3000, 5000]
delta_list = [100, 1000, 3000]
#LEF_separations = [100, 500, 3000]
LEF_separations = [74]
sites_per_monomer = 10
replication_number = 1
monomer_per_replica = 1000

steps_b = 200
steps = steps_b/pause_multiplier

already_processed = []
for fname  in glob.glob('/home1/start/polychrom/projects/targeted_loading_cohesin/*_STALL*'):
    already_processed.append(fname.split('/')[-1])

with open(str(datetime.date.today())+'_runfile_all.txt','w') as f:
    for lifetime in lifetimes:
        for velocity_multiplier in velocities:
            for lef_birth in LEF_births:
                for delta_ctcf in delta_list:
                    for face_stall in face_stalls:
                        for back_stall in back_stalls:
                            for LEF_separation in LEF_separations:
                                for CTCF_lifetime in CTCF_lifetimes:
                                    for CTCF_offtime in CTCF_offtimes:
                                        paramset = (
                                            'folder_face_'+str(face_stall)+
                                            '_back_'+str(back_stall)+
                                            '_Clife_'+str(CTCF_lifetime) +
                                            '_Cof_'+str(CTCF_offtime)+
                                            '_life_'+str(lifetime)+
                                            '_slife_'+str(lifetime)+
                                            '_birth_'+str(lef_birth)+
                                            '_deltactcf_'+str(delta_ctcf)+
                                            '_pause_'+str(LEF_pause/velocity_multiplier)+
                                            '_sep_'+str(LEF_separation)+
                                            '_site_'+str(sites_per_monomer)+
                                            '_monomer_'+str(monomer_per_replica)+
                                            '_replica_'+str(replication_number)+
                                            '_steps_'+str(steps)+
                                            '_vel_'+str(velocity_multiplier)
                                            )
                                        if paramset not in already_processed:
                                            f.write( paramset +'\n')
                                        else:
                                            print('already done')
                                
            

paramdict_keys={
                'CTCF_facestall':'face',
                'CTCF_backstall':'back',
                'CTCF_lifetime':'Clife',
                'CTCF_offtime':'Cof',
                'LEF_lifetime':'life',
                'LEF_stalled_lifetime':'slife',
                'LEF_birth':'birth',
                'LEF_pause':'pause',
                'LEF_separation':'sep',
                'sites_per_monomer':'site',
                'monomers_per_replica':'monomer',
                'number_of_replica':'replica',
                'steps':'steps',
                'velocity_multiplier':'vel'
                }
