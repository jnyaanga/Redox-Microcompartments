import numpy as np
'''Parameters of system - Inputs'''
Rc = 1.0e-5 #radius of MCP in cm
#Km values
Kmdhab = 0.9e3
Kmdhat_H = 0.06e3 #MM constant in microM
Kmdhat_P = 7.44e3 #microM
Kmdhat_N = 0.099e3 #microM
Kmdhat_D = 0.23e3 #microM
Kmicde_I = 0.6e3
Kmicde_A = 7.44e3
Kmicde_N = 0.23e3
Kmicde_D = 0.099e3 #microM
#inhibition constants
Kidhat_H = 0.06e6
Kidhat_P = 7.44e6
Kidhat_N = 0.099e6
Kidhat_D = 0.23e6
Kiicde_I = 0.6e6
Kiicde_A = 7.44e6
Kiicde_N = 0.23e6
Kiicde_D = 0.099e6
#kcat values
kcat_dhab = 100.0
kcatf_dhat = 50.0 #reactions/s
kcatr_dhat = 10 #reactions/s
kcatf_icde = 50.0
kcatr_icde = 10.0
#number of enzymes active sites
Ndhab = 2000.0 
Ndhat = 2000.0
Nicde = 2000.0

DG = 1.0e-5 #diffusivity in cm^2/s - calculate using einstein stokes eq
DH = 1.0e-5 #cm^2/s
DP = 1.0e-5 #cm^2/s
DI = 1.0e-5 #cm^2/s
DA = 1.0e-5 #cm^2/s
DN = 0
DD = 0

kcG = 1.0e-4 #permeability of MCP cm/s
kcH = kcG
kcP = kcG
kcI = kcG
kcA = kcG
kcN = 0
kcD = 0

Gcyt = 30e3 #check literature for physiological conc.
Hcyt = 0 #microM external concentration
Pcyt = 0 #microM
Icyt = 30e3 #microM
Acyt = 0 #microM
Nmcp = 15e3
Dmcp = 15e3

tstop = 3000 #s final time point for integration 
numgrid = 100.0 #number of grid points the radius is split into
MCPmil = 1.0e12 #num MCP per ml


'''Parameter combinations - do not touch'''
import numpy as np
VMCP = np.divide(4.,3.)*np.pi*np.power(Rc,3) #volume of MCP
Vratio = 20000#np.divide(1-VMCP*MCPmil,VMCP*MCPmil) #ratio of MCP volume to free bath volume per MCP 
Na = 6.022e23 #avogadro's number - constant
MCPMolar = MCPmil*1000/Na

Vdhab = np.divide(kcat_dhab*Ndhab*1.0e9,(np.divide(4,3.0)*np.pi*(Rc**3)*Na))
Vfdhat = np.divide(kcatf_dhat*Ndhat*1.0e9,(np.divide(4,3.0)*np.pi*(Rc**3)*Na))
Vrdhat = np.divide(kcatr_dhat*Ndhat*1.0e9,(np.divide(4,3.0)*np.pi*(Rc**3)*Na))
Vficde = np.divide(kcatf_icde*Nicde*1.0e9,(np.divide(4,3.0)*np.pi*(Rc**3)*Na))
Vricde = np.divide(kcatr_icde*Nicde*1.0e9,(np.divide(4,3.0)*np.pi*(Rc**3)*Na))
#nondim variables
gcyt = Gcyt/Kmdhab
hcyt = Hcyt/Kmdhat_H
pcyt = Pcyt/Kmdhat_P
icyt = Icyt/Kmicde_I
acyt = Acyt/Kmicde_A
nmcp = Nmcp/Kmdhat_N
dmcp = Dmcp/Kmdhat_D
#nondimen. parameters
alpha1 = Kmdhat_N/Kidhat_N
alpha2 = Kmdhat_D/Kidhat_D
alpha3 = Kmdhat_P/Kidhat_P
alpha4 = Kmdhat_H/Kidhat_H
alpha5 = Vfdhat/Vdhab
alpha6 = Vrdhat/Vdhab
beta1 = Kmdhat_D/Kiicde_D
beta2 = Kmdhat_N/Kiicde_N
beta3 = Kmicde_D/Kiicde_D
beta4 = Kmicde_N/Kiicde_N
beta5 = Kmicde_A/Kiicde_A
beta6 = Kmicde_I/Kiicde_I
beta7 = Vficde/Vdhab
beta8 = Vricde/Vdhab
epsilon1 = Kmdhab/Kmdhat_H
epsilon2 = Kmdhab/Kmdhat_P
epsilon3 = Kmdhab/Kmdhat_N
epsilon4 = Kmdhab/Kmdhat_D
epsilon5 = Kmdhab/Kmicde_A
epsilon6 = Kmdhab/Kmicde_I
tau = Kmdhab/Vdhab
#nondimen. kc
chiG = (Kmdhab*kcG)/(Vdhab*Rc)
chiH = (Kmdhab*kcH)/(Vdhab*Rc)
chiP = (Kmdhab*kcP)/(Vdhab*Rc)
chiI = (Kmdhab*kcI)/(Vdhab*Rc)
chiA = (Kmdhab*kcA)/(Vdhab*Rc)
chiN = (Kmdhab*kcN)/(Vdhab*Rc)
chiD = (Kmdhab*kcD)/(Vdhab*Rc)
#nondim. Diffusion
xiG = (Kmdhab*DG)/(Vdhab*(Rc**2))
xiH = (Kmdhab*DH)/(Vdhab*(Rc**2))
xiP = (Kmdhab*DP)/(Vdhab*(Rc**2))
xiI = (Kmdhab*DI)/(Vdhab*(Rc**2))
xiA = (Kmdhab*DA)/(Vdhab*(Rc**2))
xiN = (Kmdhab*DN)/(Vdhab*(Rc**2))
xiD = (Kmdhab*DD)/(Vdhab*(Rc**2))

dm = 1/float(numgrid)

tfinal = int(np.divide(np.multiply(tstop, epsilon1),tau))

nM = 3*np.int(np.reciprocal(dm))+3
initials = np.zeros(nM)

p = np.array([alpha1, alpha2, alpha3, alpha4, beta1, beta2, beta3, beta4, beta5, beta6, alpha5, alpha6, beta7, beta8, epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, epsilon6, chiG, chiH, chiP, chiI, chiA, chiN, chiD, xiG, xiH, xiP, xiI, xiA, xiN, xiD, gcyt, hcyt, pcyt, icyt, acyt, nmcp, dmcp, Vratio, dm])    