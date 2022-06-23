#!/usr/bin/python3
import numpy as np
import inspect as isp

#-------------------------------------------------------------------------------
#For reading in the template template CMS
def writecms(tmpl_cms,pos_dict,pbc_value,new_cms_path):

    #define the coordinate start line
    all_cord_line = '1 26 39.112492 -18.447816 25.110272 5 S "P" 38 -0.94390 -0.94390 "SER " " N  " "    " 7 0 2 2E2EFF  "N1"  0 1 7 <> 0 0 0 1 33.61 1 1 1 2 1 5 1 5 <> <> <> <>'
    solute_start = '1 26 39.112492 -18.447816 25.110272 5 S "P" 38 -0.94390 -0.94390 "SER " " N  " "    " 7 0 2 2E2EFF  "N1"  0 1 7 <> 0 0 0 1 33.61 1 1 1 2 1 5 1 5 <> <> <> <>'
    sodium_start = '1 66 6.309163 89.974742 -13.766686 1 4 1.00000 1.00000 "NA  " "NA  " 11 1 1E1EE1  1  24799 1 2 0'
    solvent_start = '1 16 -8.268000 3.827300 5.370200 1 70 -0.82000 -0.82000 "SPC " " O  " 8 FF2F2F  25154 1 3 2'
    stop_line = ":::"

    with open(tmpl_cms,'r') as fh1, open(new_cms_path,'w') as fh2:
        decimal_places = 6
        edit_file = False
        pass_val = False

        for line in fh1.readlines():

            if all_cord_line in line and not pass_val:
                edit_file = True
                xyz_array_i = 0 #restart counter at every point
                input_cord_length = pos_dict['all_num']
                xyz_values = pos_dict['all_xyz']
                pass_val = True

            elif solute_start in line:
                edit_file = True
                xyz_array_i = 0
                input_cord_length = pos_dict['solute_num']
                xyz_values = pos_dict['solute_xyz']

            elif sodium_start in line:
                edit_file = True
                xyz_array_i = 0
                input_cord_length = pos_dict['sodium_num']
                xyz_values = pos_dict['sodium_xyz']

            elif solvent_start in line:
                edit_file = True
                xyz_array_i = 0
                input_cord_length = pos_dict['solvent_num']
                xyz_values = pos_dict['solvent_xyz']

            if edit_file and stop_line in line:
                edit_file = False

                if xyz_array_i != input_cord_length: #basically is the xyz_cords i fed in the same length as the template file? if not, exit and throw error
                    print("::: ARRAY LENGTH ERROR (did not use all coordinates) :::"); exit()
                    #print("error in array length")
                    #pass #pass is for testing so dont throw error

            if edit_file:
                #print(line)
                ind = line.index('"')
                #newline = list(map(float, line[:ind].split())) #actually dont need because i can just replace the strings"
                b_line2 = line[:ind] #beginning of the line not split
                b_line = line[:ind].split() #beginning of the line till first quote
                #print(b_line)
                e_line = line[ind:] #end of the line
                if xyz_array_i > input_cord_length:
                    print(xyz_array_i, input_cord_length)
                    print("::: ARRAY LENGTH ERROR (not enough coordinates) :::"); exit()
                    #xyz_array_i += 1
                    #print("error in array length2")
                    pass
                else:
                    line = "{:."+str(decimal_places)+"f}"
                    b_line[2]= line.format(xyz_values[xyz_array_i,0]) #need to cast as a string
                    b_line[3]= line.format(xyz_values[xyz_array_i,1]) #need to cast as a string
                    b_line[4]= line.format(xyz_values[xyz_array_i,2]) #need to cast as a string
                    xyz_array_i += 1

                rejoined_b = "  "+" ".join(b_line)+" " #need to join two spaces for spacing, b_line is unchanged if it cant find new coordinates
                rejoined_full = rejoined_b+e_line

                fh2.write(rejoined_full)

            else:
                if "replace_this_with_pbc" in line:
                    line = "  "+str(pbc_value)+'\n'
                fh2.write(line)

    fh1.close()
    fh2.close()

#-------------------------------------------------------------------------------
#for writing the new msj file with correct cfg name
def writemsj(new_msj_path, cfg_name):
    msj_guts = isp.cleandoc("""# Desmond standard NPT relaxation protocol
    # All times are in the unit of ps.
    # Energy is in the unit of kcal/mol.
    task {
       task = "desmond:auto"
       set_family = {
          desmond = {
             checkpt.write_last_step = no
          }
       }
    }

    simulate {
       title       = "Brownian Dynamics NVT, T = 10 K, small timesteps, and restraints on solute heavy atoms, 100ps"
       annealing   = off
       time        = 100
       timestep    = [0.001 0.001 0.003 ]
       temperature = 10.0
       ensemble = {
          class = "NVT"
          method = "Brownie"
          brownie = {
             delta_max = 0.1
          }
       }
       restrain = {
          atom = "solute_heavy_atom"
          force_constant = 50.0
       }
    }

    simulate {
       effect_if   = [["==" "-gpu" "@*.*.jlaunch_opt[-1]"] 'ensemble.method = Langevin']
       title       = "NVT, T = 10 K, small timesteps, and restraints on solute heavy atoms, 12ps"
       annealing   = off
       time        = 12
       timestep    = [0.001 0.001 0.003]
       temperature = 10.0
       restrain    = { atom = solute_heavy_atom force_constant = 50.0 }
       ensemble    = {
          class  = NVT
          method = Berendsen
          thermostat.tau = 0.1
       }

       randomize_velocity.interval = 1.0
       eneseq.interval             = 0.3
       trajectory.center           = []
    }

    simulate {
       title       = "NPT, T = 10 K, and restraints on solute heavy atoms, 12ps"
       effect_if   = [["==" "-gpu" "@*.*.jlaunch_opt[-1]"] 'ensemble.method = Langevin']
       annealing   = off
       time        = 12
       temperature = 10.0
       restrain    = retain
       ensemble    = {
          class  = NPT
          method = Berendsen
          thermostat.tau = 0.1
          barostat  .tau = 50.0
       }

       randomize_velocity.interval = 1.0
       eneseq.interval             = 0.3
       trajectory.center           = []
    }

    solvate_pocket {
       should_skip = true
       ligand_file = ?
    }

    simulate {
       title       = "NPT and restraints on solute heavy atoms, 12ps"
       effect_if   = [["@*.*.annealing"] 'annealing = off temperature = "@*.*.temperature[0][0]"'
                      ["==" "-gpu" "@*.*.jlaunch_opt[-1]"] 'ensemble.method = Langevin']
       time        = 12
       restrain    = retain
       ensemble    = {
          class  = NPT
          method = Berendsen
          thermostat.tau = 0.1
          barostat  .tau = 50.0
       }

       randomize_velocity.interval = 1.0
       eneseq.interval             = 0.3
       trajectory.center           = []
    }

    simulate {
       title       = "NPT and no restraints, 24ps"
       effect_if   = [["@*.*.annealing"] 'annealing = off temperature = "@*.*.temperature[0][0]"'
                      ["==" "-gpu" "@*.*.jlaunch_opt[-1]"] 'ensemble.method = Langevin']
       time        = 24
       ensemble    = {
          class  = NPT
          method = Berendsen
          thermostat.tau = 0.1
          barostat  .tau = 2.0
       }

       eneseq.interval   = 0.3
       trajectory.center = solute
    }

    simulate {
       cfg_file = "%s"
       jobname  = "$MASTERJOBNAME"
       dir      = "."
       compress = ""
    }"""%(cfg_name))

    with open(new_msj_path,'w') as fh:
        fh.write(msj_guts)
        fh.close()

#-------------------------------------------------------------------------------
#for writing a new cfg file with correct info
def writecfg(new_cfg_path, seed):
    cfg_guts = isp.cleandoc("""annealing = false
    backend = {
    }
    bigger_rclone = false
    checkpt = {
       first = 0.0
       interval = 240.06
       name = "$JOBNAME.cpt"
       write_last_step = true
    }
    cpu = 1
    cutoff_radius = 9.0
    elapsed_time = 0.0
    energy_group = false
    eneseq = {
       first = 0.0
       interval = 1.2
       name = "$JOBNAME$[_replica$REPLICA$].ene"
    }
    ensemble = {
       barostat = {
          tau = 2.0
       }
       class = NPT
       method = MTK
       thermostat = {
          tau = 1.0
       }
    }
    glue = solute
    maeff_output = {
       first = 0.0
       interval = 120.0
       name = "$JOBNAME$[_replica$REPLICA$]-out.cms"
       periodicfix = true
       trjdir = "$JOBNAME$[_replica$REPLICA$]_trj"
    }
    meta = false
    meta_file = ?
    pressure = [1.01325 isotropic ]
    randomize_velocity = {
       first = 0.0
       interval = inf
       seed = %4i
       temperature = "@*.temperature"
    }
    restrain = none
    simbox = {
       first = 0.0
       interval = 1.2
       name = "$JOBNAME$[_replica$REPLICA$]_simbox.dat"
    }
    surface_tension = 0.0
    taper = false
    temperature = [
       [300.0 0 ]
    ]
    time = 10000.0
    timestep = [0.002 0.002 0.006 ]
    trajectory = {
       center = []
       first = 0.0
       format = dtr
       frames_per_file = 250
       interval = 10.0
       name = "$JOBNAME$[_replica$REPLICA$]_trj"
       periodicfix = true
       write_velocity = false
    }""" %(seed))

    with open(new_cfg_path,'w') as fh:
        fh.write(cfg_guts)
        fh.write('\n\n') #need so that the last line is blank to match the schrodinger format
        fh.close()

#-------------------------------------------------------------------------------
#write msj with no relaxation stages
def writemsj_norelax(new_msj_path, cfg_name):
    msj_guts = isp.cleandoc("""# Desmond standard NPT relaxation protocol
    # All times are in the unit of ps.
    # Energy is in the unit of kcal/mol.
    task {
       task = "desmond:auto"
       set_family = {
          desmond = {
             checkpt.write_last_step = no
          }
       }
    }

    simulate {
       cfg_file = "%s"
       jobname  = "$MASTERJOBNAME"
       dir      = "."
       compress = ""
    }"""%(cfg_name))

    with open(new_msj_path,'w') as fh:
        fh.write(msj_guts)
        fh.close()
