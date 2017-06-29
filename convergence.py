import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pycbc.io import InferenceFile
import sys

print "Initialising..."

## Select file
folder="convtest1"
#parameter=sys.argv[2]
folder="final/"+folder+"/"

## Combine inputs to form variables
data_name="output.hdf"
num_walkers=5000
no_steps=50

## Function to extract posterior for a given parameter
## at a given iteration
def getParameter(parameter,iteration):
   ## Prepare to read in parameters
   datafile=folder+data_name
   fp = InferenceFile("%s" % datafile, "r")
   
   ## Take last iteration of each walker
   parameter_values=np.array([])
   for aa in range(num_walkers):
      samples = fp.read_samples("%s" % parameter, walkers=aa)
      temp=getattr(samples,parameter)
      parameter_values=np.append(parameter_values,temp[iteration])
   return parameter_values

def chi_prec(it_no):
   ## chi_p is given by (1/B1m1^2)*max(B1*S1perp,B2*S2perp)
   ## with B1=2+3/2q, B2=2+3q/2
   ## so we need m1, q, s1_a, s1_polar, s2_a, s2_polar
   ## NB chi_p should always be 0 < chi_p < 1
   chi_p=np.zeros(num_walkers)
   print "   Extracting intrinsic parameters..."
   ## Generate arrays for each parameter
   s1_a=getParameter("spin1_a",it_no)
   s1_polar=getParameter("spin1_polar",it_no)
   s2_a=getParameter("spin2_a",it_no)
   s2_polar=getParameter("spin2_polar",it_no)
   mass1=getParameter("mass1",it_no)
   mass2=getParameter("mass2",it_no)
   print "   Calculating derived parameters..."
   for aa in range(len(mass1)):  ## Standard chi_p function]
      if mass1[aa]>mass2[aa]:
         ratio=mass2[aa]/mass1[aa]
         B1=2+((3*ratio)/2)
         B2=2+(3/(ratio*2))
         spin1_plane=s1_a[aa]*np.sin(s1_polar[aa])
         spin2_plane=s2_a[aa]*np.sin(s2_polar[aa])
         arg1=B1*spin1_plane*mass1[aa]*mass1[aa]
         arg2=B2*spin2_plane*mass2[aa]*mass2[aa]
         chi_p[aa]=(max(arg1,arg2))/(mass1[aa]*mass1[aa]*B1)
      else:
         ratio=mass2[aa]/mass1[aa] # Modify function for inverted mass ratio
         B1=2+((3*ratio)/2)
         B2=2+(3/(ratio*2))
         spin1_plane=s1_a[aa]*np.sin(s1_polar[aa]) # Spin1 is smaller mass this time!
         spin2_plane=s2_a[aa]*np.sin(s2_polar[aa]) # Spin2 is larger mass this time!
         arg1=B1*spin1_plane*mass1[aa]*mass1[aa]   # Swap the B coefficients now as B1 should be on the larger mass
         arg2=B2*spin2_plane*mass2[aa]*mass2[aa]
         chi_p[aa]=(max(arg1,arg2))/(mass2[aa]*mass2[aa]*B2)
   return chi_p

def chi_effect(it_no):
   ## chi_eff is given by (S1/m1+S2/m2) dot L/M where M is total mass
   ## So for this we need m1, m2, s1_a, s2_a, s1_polar, s2_polar... fak me
   chi_eff=np.zeros(num_walkers)
   print "   Extracting intrinsic parameters..."
   ## Generate arrays for each paramter
   s1_a=getParameter("spin1_a",it_no)
   s1_polar=getParameter("spin1_polar",it_no)
   s2_a=getParameter("spin2_a",it_no)
   s2_polar=getParameter("spin2_polar",it_no)
   m1=getParameter("mass1",it_no)
   m2=getParameter("mass2",it_no)
   M=m1+m2

   ## Find spins along z-axis
   s1_z=m1*m1*s1_a*np.cos(s1_polar)
   s2_z=m2*m2*s2_a*np.cos(s2_polar)

   print "   Calculating derived parameters..."
   chi_eff=(s1_z/m1+s2_z/m2)/M
   return chi_eff

## Find relative change in mean for this parameter as a function of iteration
def findMeans(parameter):
   print "Performing convergence test for %s" % parameter
   if parameter=="chi_p":
      int_step=np.linspace(0,19999,no_steps)
      for bb in range(no_steps):
         int_step[bb]=int(int_step[bb])
      means=np.zeros(no_steps)
      for aa in range(no_steps): ## Loop over iterations from start to finish
         print "%.2f %% complete" % (100*aa/no_steps)
         lower=chi_prec(aa) ## Earlier iteration
         upper=chi_prec(aa+1) ## Later iteration
         lower=np.mean(lower) ## Mean of earlier it
         upper=np.mean(upper) ## Mean of later it
         means[aa]=abs(lower-upper) ## Relative change in mean
   elif parameter=="chi_eff":
      int_step=np.linspace(0,19999,no_steps)
      for bb in range(no_steps):
         int_step[bb]=int(int_step[bb])
      means=np.zeros(no_steps)
      for aa in range(no_steps): ## Loop over iterations from start to finish
         print "%.2f %% complete" % (100*aa/no_steps)
         lower=chi_effect(aa) ## Earlier iteration
         upper=chi_effect(aa+1) ## Later iteration
         lower=np.mean(lower) ## Mean of earlier it
         upper=np.mean(upper) ## Mean of later it
         means[aa]=abs(lower-upper) ## Relative change in mean
   ## Cannot handle any other derived parameters at the moment
   else:
      int_step=np.linspace(0,19999,no_steps)
      for bb in range(no_steps):
         int_step[bb]=int(int_step[bb])
      means=np.zeros(no_steps)
      for aa in range(no_steps): ## Loop over iterations from start to finish
         print "%.2f %% complete" % (100*aa/no_steps)
         lower=getParameter(parameter,aa) ## Earlier iteration
         upper=getParameter(parameter,aa+1) ## Later iteration
         lower=np.mean(lower) ## Mean of earlier it
         upper=np.mean(upper) ## Mean of later it
         means[aa]=abs(lower-upper) ## Relative change in mean
   ## Plot results
   savename=folder+parameter
   plt.figure()
   plt.plot(int_step,means,'bx')
   plt.title("Convergence of %s" % parameter)
   plt.xlabel("Iteration number")
   plt.ylabel("Absolute change in mean")
   plt.grid()
   plt.show("hold")
   plt.savefig("%s_conv.png" % savename)

findMeans("chi_eff")
findMeans("chi_p")
