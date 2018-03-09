"""
Script to pre-process input gramicidin structure to generate file ready for simulation.
This prepares only the gramicidin structure, the memebrane or water molecules are not placed.
"""

import os
import sys
from shutil import copy2
from pymcce.automated_mcce import MCCEParams


def main():
    curr_dir = os.getcwd()
    init_struct_dir = os.path.abspath("../data/gramicidin/init_structures")
    input_dir = os.path.abspath("../data/gramicidin/prep_structures")
    mcce_dir = os.path.abspath("../code/mcce3.5")
    ipece_dir = os.path.abspath("../code/ipece")

    os.chdir(input_dir)
    init_structs = sorted([f for f in os.listdir(init_struct_dir) if not f.startswith("noh") and f.endswith("pdb")])
    print(init_structs)
    for index, st in enumerate(init_structs):
        formatted_index = "%02d" % index
        print("Processing structure: ", formatted_index)
        st_infile = os.path.join(init_struct_dir, st)
        st_writefile = os.path.join(init_struct_dir, "noh_" + st)
        # remove hydrogens
        reduce_command = "reduce -T %s > %s" % (st_infile, st_writefile)
        os.system(reduce_command)

        # replace termini
        terminal_rep = {'CY  VAL' : "CA  FOR",
                        'OY  VAL' : "O   FOR",
                        'NT  TRP' : "N   ENA",
                        'C1  TRP' : "CA  ENA",
                        'C2  TRP' : "CB  ENA",
                        'OG  TRP' : "OG1 ENA"
                        }
        input_struct_file = "%s_%s" % (formatted_index, "input_struct.pdb")
        with open(st_writefile, "r") as oldfile, open(input_struct_file, 'w') as newfile:
            for line in oldfile:
                atm_res_name = line[13:20]
                for terminal_atom in terminal_rep.keys():
                    if terminal_atom == atm_res_name:
                        new_name = terminal_rep[terminal_atom]
                        #print(line[:13] + atm_res_name + line[20:])
                        #print(line[:13] + new_name + line[20:])
                        line = line[:13] + new_name + line[20:]

                #print(line)
                newfile.write(line)

        # generate a protein step1 file
        m = MCCEParams(mcce_dir)
        m.edit_parameters(INPDB=input_struct_file)
        m.edit_parameters(DO_PREMCCE="t", DO_ROTAMERS="f", TERMINALS="f")
        m.write_runprm(input_dir + "/")
        run_command = m.mcce_directory + "/mcce > prep_struct_run.log"
        os.system(run_command)
        copy2("step1_out.pdb", formatted_index + "_step1_out.pdb")


        # insert membrane
        ipece_parm = os.path.join(ipece_dir, input_dir, "ipece.prm")
        input_pdb = formatted_index + "_step1_out.pdb"
        output_pdb = formatted_index + "_mem_step1_out.pdb"
        # assemble ipece command and run
        run_command = ipece_dir + "/ipece " + ipece_parm + " " + input_pdb + " " + output_pdb
        #print(run_command)
        os.system(run_command)

        # run step 2
        m = MCCEParams(os.path.abspath(mcce_dir))
        m.edit_parameters(DO_PREMCCE="t", DO_ROTAMERS="t", INPDB=os.path.join(input_dir, output_pdb))
        m.write_runprm(input_dir + "/")
        run_command = m.mcce_directory + "/mcce > prep_struct_run.log"
        os.system(run_command)
        copy2("step2_out.pdb", formatted_index + "_step2_out.pdb")

        # fix membrane atom names in step2
        pdb = os.path.join(input_dir, formatted_index + "_step2_out.pdb")
        fixed_pdb = os.path.join(input_dir, formatted_index + "_final_step2_out.pdb")
        with open(pdb, "r") as f1:
            with open(fixed_pdb, "w") as f2:
                lines = f1.readlines()
                for l in lines:
                    if "EM_ _" in l:
                        l = l.replace("x", "C")
                        l = l.replace("EM_ _", "MEM M")
                        f2.write(l)
                    else:
                        f2.write(l)
        break
    os.chdir(curr_dir)


if __name__ == "__main__":
    sys.exit(main())
