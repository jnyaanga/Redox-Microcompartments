import numpy as np
import Params

def derivAve(allvar, time):
    nsteps = int(np.floor(np.divide(1,Params.dm)))
    d = np.zeros(14)
    
    variables = allvar[:-7]
    ext = allvar[-7:]
    allder = np.zeros(21)
    
    RdhaB = np.divide(variables[0],1+variables[0])
   
    RdhaT = np.divide((Params.alpha1*Params.alpha5*variables[5]*variables[1]) - (Params.alpha2*Params.alpha6*variables[2]*variables[6]),1+Params.alpha1*(variables[5]+variables[1]+(variables[5]*variables[1])) + Params.alpha2*(variables[2]+variables[6]+(variables[2]*variables[6])) + Params.alpha1*Params.alpha2*((variables[5]*variables[2])+(variables[1]*variables[6])) + Params.alpha1*Params.alpha3*(variables[1]*variables[5]*variables[2]) + Params.alpha2*Params.alpha4*(variables[1]*variables[2]*variables[6]))
    
    RicdE = np.divide((Params.beta1*Params.beta7*variables[6]*variables[3]) - (Params.beta2*Params.beta8*variables[4]*variables[5]),1+Params.beta1*(variables[6]+(variables[6]*variables[3])) + Params.beta3*variables[3] + Params.beta2*(variables[5]+(variables[4]*variables[5])) + Params.beta4*variables[4] + Params.beta1*Params.beta4*(variables[6]*variables[4]) + Params.beta2*Params.beta3*(variables[3]*variables[5]) + Params.beta1*Params.beta5*(variables[6]*variables[3]*variables[4]) + Params.beta2*Params.beta6*(variables[3]*variables[4]*variables[5]))
    
#inside -
#G
    d[0] = np.multiply(9*((nsteps-0.5)**(4/3.0))*(variables[-7] - variables[0])*Params.xiG,Params.dm**(1/3.0)) - RdhaB
#H
    d[1] = np.multiply(9*((nsteps-0.5)**(4/3.0))*(variables[-6] - variables[1])*Params.xiH,Params.dm**(1/3.0)) + np.multiply(Params.epsilon1,(RdhaB - RdhaT))             
#P
    d[2] = np.multiply(9*((nsteps-0.5)**(4/3.0))*(variables[-5] - variables[2])*Params.xiP,Params.dm**(1/3.0)) + Params.epsilon2*RdhaT
#I    
    d[3] = np.multiply(9*((nsteps-0.5)**(4/3.0))*(variables[-4] - variables[3])*Params.xiI,Params.dm**(1/3.0)) - Params.epsilon6*RicdE
#A    
    d[4] = np.multiply(9*((nsteps-0.5)**(4/3.0))*(variables[-3] - variables[4])*Params.xiA,Params.dm**(1/3.0)) + Params.epsilon5*RicdE
#N    
    d[5] = np.multiply(9*((nsteps-0.5)**(4/3.0))*(variables[-2] - variables[5])*Params.xiN,Params.dm**(1/3.0)) + np.multiply(Params.epsilon3,(RicdE - RdhaT))
#D    
    d[6] = np.multiply(9*((nsteps-0.5)**(4/3.0))*(0 - variables[6])*Params.xiD,Params.dm**(1/3.0)) + np.multiply(Params.epsilon4,(RdhaT - RicdE))
                                  
    
#boundary - 
#G
    d[-7] = np.divide(3*((nsteps+0.5)**(2/3.0))*(ext[0] - variables[-7])*Params.chiG,Params.dm**(1/3.0)) - np.divide(9*((nsteps-0.5)**(4/3.0))*(variables[-7] - variables[0])*Params.xiG,Params.dm**(2/3.0))
#H
    d[-6] = np.divide(3*((nsteps+0.5)**(2/3.0))*(ext[1] - variables[-6])*Params.chiH,Params.dm**(1/3.0)) - np.divide(9*((nsteps-0.5)**(4/3.0))*(variables[-6] - variables[1])*Params.xiH,Params.dm**(2/3.0))
#P    
    d[-5] = np.divide(3*((nsteps+0.5)**(2/3.0))*(ext[2] - variables[-5])*Params.chiP,Params.dm**(1/3.0)) - np.divide(9*((nsteps-0.5)**(4/3.0))*(variables[-5] - variables[2])*Params.xiP,Params.dm**(2/3.0))
#I    
    d[-4] = np.divide(3*((nsteps+0.5)**(2/3.0))*(ext[3] - variables[-4])*Params.chiI,Params.dm**(1/3.0)) - np.divide(9*((nsteps-0.5)**(4/3.0))*(variables[-4] - variables[3])*Params.xiI,Params.dm**(2/3.0))
#A    
    d[-3] = np.divide(3*((nsteps+0.5)**(2/3.0))*(ext[4] - variables[-3])*Params.chiA,Params.dm**(1/3.0)) - np.divide(9*((nsteps-0.5)**(4/3.0))*(variables[-3] - variables[4])*Params.xiA,Params.dm**(2/3.0))
#N
    d[-2] = np.divide(3*((nsteps+0.5)**(2/3.0))*(ext[5] - variables[-2])*Params.chiN,Params.dm**(1/3.0)) - np.divide(9*((nsteps-0.5)**(4/3.0))*(variables[-2] - variables[5])*Params.xiN,Params.dm**(2/3.0))
#D    
    d[-1] = np.divide(3*((nsteps+0.5)**(2/3.0))*(ext[6] - variables[-1])*Params.chiD,Params.dm**(1/3.0)) - np.divide(9*((nsteps-0.5)**(4/3.0))*(variables[-1] - variables[6])*Params.xiD,Params.dm**(2/3.0))
    
#external - 
    dext = np.zeros(7)

#G
    dext[0] = (-np.divide(3*((nsteps+0.5)**(2/3.0))*(ext[0] - variables[-7])*Params.chiG,Params.dm**(1/3.0)))*np.divide(Params.dm,Params.Vratio)
#H
    dext[1] = (-np.divide(3*((nsteps+0.5)**(2/3.0))*(ext[1] - variables[-6])*Params.chiH,Params.dm**(1/3.0)))*np.divide(Params.dm,Params.Vratio)
#P
    dext[2] = (-np.divide(3*((nsteps+0.5)**(2/3.0))*(ext[2] - variables[-5])*Params.chiP,Params.dm**(1/3.0)))*np.divide(Params.dm,Params.Vratio)
#I
    dext[3] = (-np.divide(3*((nsteps+0.5)**(2/3.0))*(ext[3] - variables[-4])*Params.chiI,Params.dm**(1/3.0)))*np.divide(Params.dm,Params.Vratio)
#A
    dext[4] = (-np.divide(3*((nsteps+0.5)**(2/3.0))*(ext[4] - variables[-3])*Params.chiA,Params.dm**(1/3.0)))*np.divide(Params.dm,Params.Vratio)
#N
    dext[5] = (-np.divide(3*((nsteps+0.5)**(2/3.0))*(ext[5] - variables[-2])*Params.chiN,Params.dm**(1/3.0)))*np.divide(Params.dm,Params.Vratio)
#D
    dext[6] = (-np.divide(3*((nsteps+0.5)**(2/3.0))*(ext[6] - variables[-1])*Params.chiD,Params.dm**(1/3.0)))*np.divide(Params.dm,Params.Vratio)
    
    allder[:-7] = d
    allder[-7:] = dext
    return allder  