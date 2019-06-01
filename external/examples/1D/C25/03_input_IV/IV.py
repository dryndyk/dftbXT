#> Calculates the current-voltage curve.
#> Calculates the T(E,V) and LDOS(E,V).
#> The initial DFTB+XT input file is asumed to be dftb_in0.hsd!

import subprocess
import math
import os

#--------------------------------------------------------------------------------------------------#
# Input parameters (to be replaced by your input values)
#--------------------------------------------------------------------------------------------------#

#> The .gen file for configuration.
genFileName = 'C25.gen' #'<filename>'

#> The Fermi level [eV].
EF =  -5.6407

#> The first voltage [Volt].
V1 = 0.

#> The last voltage [Volt].
V2 = 1.

#> The voltage step [Volt].
Vstep = 0.1

NumMPI = 100
dftbXT = 'dftbXT'

#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#

#out_current = open('current.dat','w')
#out_current.close()
#out_dIdV = open('dIdV.dat','w')
#out_dIdV.close()

Nv = int((V2-V1)/Vstep)+1

dIdV = 0.
Iprev = 0.

for i in range(1,Nv+1):
    
  V = V1+(i-1)*Vstep 
  
  os.system("rm GS/*")
  
  subprocess.call(['echo',''])
  subprocess.call(['echo','Voltage =', '%f' % V])
  subprocess.call(['echo',''])

  #----------------------------------------------------------------------------#
  # Change the dftb_in.hsd file.
  #----------------------------------------------------------------------------#
  
  inFile = open('dftb_in0.hsd','r')
  outFile = open('dftb_in.hsd','w')
 
  k = 0
  for line in inFile:
    words = line.split()
    if len(words) > 0:
      if words[0] == '<<<':
        outFile.write('  <<< '+genFileName)
        outFile.write('\n')  
      elif words[0] == 'Potential':
        k = k+1  
        if k == 1:
          outFile.write('    Potential [eV] = 0.0')
        if k == 2:
          outFile.write('    Potential [eV] = '+str(-V))  
        outFile.write('\n')
      else:
        outFile.write(line)
              
  inFile.close()    
  outFile.close()

  #----------------------------------------------------------------------------#
  
  #subprocess.call(['cp','dftb_in.hsd','dftb_in_'+str(V)+'.hsd'])
  
  #----------------------------------------------------------------------------#
  # dftb+xt calculation
  #----------------------------------------------------------------------------#  
    
  out1 = open('out1','w')   
  subprocess.call(['echo','dftbXT is started'])
  #subprocess.call(['mpirun','-np','8','dftbXT'])
  subprocess.call(['mpirun','-np',str(NumMPI),dftbXT],stdout=out1)
  out1.close()

  #----------------------------------------------------------------------------#
  # Current and dI/dV from output
  #----------------------------------------------------------------------------#
 
  out1 = open('out1','r')
  
  for line in out1:
    words = line.split()
    if len(words) > 0:
      if words[0] == 'contacts:':
        I = float(words[4])
                  
  out1.close()
  
  dIdV = (I-Iprev)/Vstep

  out_current = open('current.dat','a')
  out_current.write('%19.10e' % V)
  out_current.write('%19.10e' % I)
  out_current.write('\n')
  out_current.close()
  out_dIdV = open('dIdV.dat','a')
  out_dIdV.write('%19.10e' % V)
  out_dIdV.write('%19.10e' % dIdV)
  out_dIdV.write('\n')
  out_dIdV.close()
  
  Iprev = I
  
  #----------------------------------------------------------------------------#
  # from transmission
  #----------------------------------------------------------------------------#

  out_T = open('T.dat','a')
  tun = open('transmission.dat','r')   
  for line in tun:
    words = line.split()
    if len(words) > 0.:
      if float(words[0]) > EF+V-0.001:
        out_T.write('%19.10e' % (V))
        out_T.write('   ')
        out_T.write(str(float(words[1])))
        out_T.write('\n')
        break
  out_T.close()
  tun.close()
     
  #----------------------------------------------------------------------------#   

  subprocess.call(['cp','transmission.dat','transmission_'+str(V)+'.dat']) 
  subprocess.call(['cp','localDOS.dat','localDOS_'+str(V)+'.dat'])
  subprocess.call(['cp','charges.bin','charges_'+str(V)+'.bin'])

  os.system("cat out1 >> out")
  subprocess.call(['rm','out1'])
  subprocess.call(['rm','transmission.dat'])
  subprocess.call(['rm','localDOS.dat'])
  subprocess.call(['rm','dftb_pin.hsd'])
  subprocess.call(['rm','detailed.out'])
  #subprocess.call(['rm','band.out'])
  

subprocess.call(['rm','dftb_in.hsd'])
subprocess.call(['echo',''])
subprocess.call(['echo','Bye!'])
subprocess.call(['echo',''])

