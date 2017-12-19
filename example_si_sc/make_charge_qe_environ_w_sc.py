import os
import shlex
import shutil
import subprocess
import numpy as np


def main():
	tot_charge = ['0.001','0.003','0.005','0.007','0.009','0.013','0.018','0.023']

	epsilon = 11.7
	prefix = 'SiO'

	for q in tot_charge:
		os.mkdir("q_%s"%(q))
		os.chdir("q_%s"%(q))
		make_pw_in(q)
		make_environ_in(q,epsilon)
		shutil.copy("/storage/home/qjc5019/Group/test_files/2H/q_0.001/SiSi2H.q0.001.pbs","./%s.q%s.pbs"%(prefix,q))
		subprocess.call(["qsub","%s.q%s.pbs"%(prefix,q)])
		os.chdir("./../")


	return



def make_pw_in(q):
	f = open('scf.in','w')
	q_num = float(q)
	f.write("""&CONTROL
  calculation='scf',
  prefix='Si-O',
  pseudo_dir='/storage/home/qjc5019/Pseudopotentials/Pseudos/GBRV_USPP_Library/GBRV_USPP_PBE_UPF_format',
  verbosity='high',
  restart_mode='from_scratch',
  wf_collect = .TRUE.,
  nstep=200
  outdir= '/storage/home/qjc5019/scratch/tmp_Si-O_alt_1-side_relax',
/

&SYSTEM
  ibrav=0,
  celldm(1)=7.2565749368d0,
  nat=15,
  ntyp=2,
  ecutwfc=50.0d0,
  ecutrho=200.0d0,
  occupations='smearing',
  smearing='mv',
  degauss=0.03000d0,
  assume_isolated = 'slabz'
  tot_charge = 0.0,
/
&ELECTRONS
  diagonalization='david',
  conv_thr=1d-10,
  mixing_mode='local-TF',
  mixing_beta=0.500d0,
  electron_maxstep = 250
/
&IONS
  !ion_dynamics = 'damp'
/


ATOMIC_SPECIES
  Si 28.085500d0 si_pbe_v1.uspp.F.UPF
  O  15.999000d0 o_pbe_v1.2.uspp.F.UPF 

ATOMIC_POSITIONS {alat}
Si       0.000000000   0.000000000   1.802073905    0   0   0
Si       0.500000000   1.060724541   1.802078557    0   0   0
Si       0.000000000   0.353553391   2.302078557    0   0   0
Si       0.499909714   0.707106781   2.302060048    0   0   0
Si       0.000000000   0.000000000   2.802078557    0   0   0
Si       0.500000000   1.059951789   2.802078557    0   0   0
Si      -0.000852798   0.353553391   3.302078557    0   0   0
Si       0.500000000   0.707106781   3.302090970    0   0   0
Si       0.000000000   0.000000000   3.802078557    0   0   0
Si       0.500000000   1.060660172   3.802078557    0   0   0
Si       0.011264507   0.346111979   4.314003460
Si       0.492983512   0.724485161   4.333979006
Si       0.127614179  -0.162293223   4.668709662
Si       0.672109676   0.916491069   4.928620199
O        0.428409191   1.299131531   4.980203547

K_POINTS {automatic}
  5 5 1 1 1 1

CELL_PARAMETERS {alat}
  1.000000000000d0  0.000000000000d0  0.000000000000d0
  0.000000000000d0  1.414213562374d0  0.000000000000d0
  0.000000000000d0  0.000000000000d0  7.624988536903d0""")
	f.close()





	return

def make_environ_in(q,epsilon):
	f = open("environ.in",'w')
	q_num = float(q)
	f.write("""&Environ
   environ_type = 'water',
   verbose = 5,
   eps_mode = 'full',
   tolrhopol = 1d-13,
   !environ_restart = .true.,
   semiconductor_region = .TRUE.,
/
&semiconductor
   dopant_type = 'n',
   dopant_concentration = 1d18,
   slab_direction = 'z',
   sc_dielectric = %s,
   electrode_charge = %s,
   !slab_distance = 10
   sc_cutoff = 33.07,
   !sc_cutoff_spread = 1.0,
   flatband_calc = 'from_scratch',
   !flatband_fermi = -5.72,
   !chg_thr =1.0D-6
/"""%(str(epsilon),str(q)))
	f.close()



	return









if __name__=='__main__':
    main()



