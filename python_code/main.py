
# coding: utf-8

# In[2]:


import math as m
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from sympy.abc import o,n

# INPUTS OF THE SCRIPT:
Uinf = 1
profile = "NACA 2408"
while type(profile) == str:
    try:
        a = int(profile)
        break
    except:
        profile = profile[1:]    

if profile[0]==" ":
    profile = profile[1:]

# PROFILE PARAMETERS
# We supose an unitary chord, with no loss of generality.
f = int(profile[0])/100
p = int(profile[1])/10
t = int(profile[2:])/100


# FUNTIONS' DEFINITION
def DVM(M,alpha,f,p,t,Uinf,xh=1,eta=0):
    N = M+1
    alpha = alpha * m.pi / 180
    eta = eta * m.pi / 180
    # PROFILE'S DESCRIPTION POINTS
    xilist = []
    zilist = []
    for i in range(0, N):
        #xi = i*(1/(N-1))    # Uncomment this line for trying a uniform distribution
        xi = 0.5 * (1-m.cos(m.pi*i/(N-1)))
        xilist.append(xi)
        if xi < p:
            zi = (f*(2*p*xi - (xi**2)))/(p**2)    
        else:
            zi = (f*(1-2*p+2*p*xi-(xi**2))) / ((1-p)**2)
        zilist.append(zi)
    
    # FLAP POINTS COORDINATE'S CORRECTION
    distance = 10
    for xi in xilist:
        if abs(xi-xh) < distance:
            distance = abs(xi-xh)
            xhfinal = xi
            zhfinal = zilist[xilist.index(xhfinal)]
    for i in range(xilist.index(xhfinal)+1,len(xilist)):
        xilist[i] = xhfinal + (xilist[i] - xhfinal)*m.cos(eta) + (zilist[i] - zhfinal)*m.sin(eta)
        zilist[i] = zhfinal - (xilist[i] - xhfinal)*m.sin(eta) + (zilist[i] - zhfinal)*m.cos(eta)
        
    # CALCULATION    
    cilist = []
    nmatr = np.zeros((M,2))
    tmatr = np.zeros((M,2))
    rvmatr = np.zeros((M,2))
    rcpmatr = np.zeros((M,2))

    for j in range(0,N-1):
        ci = m.sqrt((xilist[j+1] - xilist[j])**2 + (zilist[j+1] - zilist[j])**2)
        cilist.append(ci)
        nix = -(zilist[j+1] - zilist[j])/ci
        niy = (xilist[j+1] - xilist[j])/ci
        nmatr[j][0] = nix
        nmatr[j][1] = niy
        tix = (xilist[j+1] - xilist[j])/ci
        tiy = (zilist[j+1] - zilist[j])/ci
        tmatr[j][0] = tix
        tmatr[j][1] = tiy
        xvi = 0.25*ci*tmatr[j][0] + xilist[j]
        zvi = 0.25*ci*tmatr[j][1] + zilist[j]
        rvmatr[j][0] = xvi
        rvmatr[j][1] = zvi
        xcpi = 0.75*ci*tmatr[j][0] + xilist[j]
        zcpi = 0.75*ci*tmatr[j][1] + zilist[j]
        rcpmatr[j][0] = xcpi
        rcpmatr[j][1] = zcpi

    A = np.zeros((M,M))   
    RHSmatr = np.zeros((M,1))
    gammamatr = np.zeros((M,1))

    for i in range(0,M):
        for j in range(0,M):
            rsq = (rcpmatr[i][0]-rvmatr[j][0])**2 + (rcpmatr[i][1]-rvmatr[j][1])**2
            u_ij = (rcpmatr[i][1]-rvmatr[j][1])/(2*m.pi*rsq)
            w_ij = -(rcpmatr[i][0]-rvmatr[j][0])/(2*m.pi*rsq)
            A[i][j] = (u_ij*nmatr[i][0])+(w_ij*nmatr[i][1])
        RHSmatr[i][0] = -Uinf*(m.cos(alpha)*nmatr[i][0] + m.sin(alpha)*nmatr[i][1])

    Ainv = np.linalg.inv(A)   
    gammamatr = np.dot(Ainv, RHSmatr)

    Cl = 0
    for l in gammamatr:
        Cl = Cl + float(l)
    Cl = Cl * 2 /Uinf 

    Cm = 0
    for j in range(0,len(xilist)-1):
        gamma = gammamatr[j][0]
        xj = xilist[j]
        Cm = Cm + gamma * xj * m.cos(alpha)
        
    Cm = Cm*(-2)/Uinf


    Cp = [] 
    for j in range(0,len(xilist)-1):
        gamma = gammamatr[j][0]
        cj = cilist[j]
        Cp.append(2 * gamma / (Uinf * cj))
    return(Cl, Cm, Cp, xilist, zilist)
#                                                END OF FUNCTION: DVM 

def find_alpha0(lim_inf,lim_sup,xh,eta):
    cl = 100
    result_sup = DVM(256, lim_sup, f, p, t, Uinf,xh,eta)[0]
    result_inf = DVM(256, lim_inf, f, p, t, Uinf,xh,eta)[0]
    if (result_sup * result_inf) > 0:
        print('Error: Bad limits. Try again with other values')
    while abs(cl) > 1e-4:
        alpha = 0.5*( lim_sup + lim_inf)
        result = DVM(256,alpha,f,p,t,Uinf,xh,eta)
        result_sup = DVM(256, lim_sup, f, p, t, Uinf,xh,eta)
        result_inf = DVM(256, lim_inf, f, p, t, Uinf,xh,eta)
        cl = result[0]
        cl_sup = result_sup[0]
        cl_inf = result_inf[0]
        if cl*cl_sup > 0:
            lim_sup = alpha
        else:
            lim_inf = alpha

    return(alpha, cl)
# END OF FUNCTION: find_alpha0

def find_cl_alpha(xh, eta):
    cl_lin_sup = DVM(256,10,f,p,t,Uinf,xh,eta)[0]
    cl_lin_inf = DVM(256,-10,f,p,t,Uinf,xh,eta)[0]
    cla = ((cl_lin_sup - cl_lin_inf)/20)  # in º^-1
    return(cla)
# END OF FUNCTION: find_cl_alpha

def find_cm_0(alpha,xh,eta):
    return(DVM(256, alpha, f, p, t, Uinf, xh, eta)[1])
# END OF FUNCTION find_cm_0


# In[2]:


# ---------------------------------------- CALCULATIONS AND RESULTS -----------------------------------
#      1. Verification assessment 
# 1.1. Analitical solution:
fun = (2*p - 1 + sp.cos(o))*sp.cos(n*o)
O = m.acos(1-2*p)
A0 = 4*m.pi/180 -((f/p**2)*sp.integrate(fun, (o,0,O)).evalf(subs={n:0}) + (f/(1-p)**2)*sp.integrate(fun, (o,O,m.pi)).evalf(subs={n:0}))/m.pi
A1 = (2/m.pi)*((f/p**2)*(sp.integrate(fun, (o,0,O)).evalf(subs={n:1})) + (f/(1-p)**2)*(sp.integrate(fun, (o,O, m.pi)).evalf(subs={n:1})))
A2 = (2/m.pi)*((f/p**2)*(sp.integrate(fun, (o,0,O)).evalf(subs={n:2})) + (f/(1-p)**2)*(sp.integrate(fun, (o,O, m.pi)).evalf(subs={n:2})))
CL = m.pi * (2*A0 +A1)
CM = -0.5*m.pi*(A0+A1-0.5*A2)

# 1.2. Numerical solution:
Mlist = [4,8,16,32,64,128,256]
resultCl = []
resultCm = []
errorCl = []
errorCm = []

for M in Mlist:
    result = DVM(M,4,f,p,t,Uinf,1,0) 
    resultCl.append(result[0])
    resultCm.append(result[1])
    errorCl.append(abs(CL-result[0]) * 100 / CL)
    errorCm.append(abs((CM-result[1])) / CM * 100)    

fig1, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.set_ylabel('Cl')
ax1.set_xlabel('N')
ax2.set_ylabel('Error (%)')
plt.xlabel("N")
ax1.plot(Mlist,resultCl, color='Blue',label='Cl')
ax2.plot(Mlist,errorCl, color='Red', label='Error')
ax1.legend(loc = 'best', bbox_to_anchor=(0.105, 0.105,0.855,0.6), frameon = False)
ax2.legend(loc = 'best', bbox_to_anchor=(0.15, 0.1,0.855,0.5), frameon = False)
plt.title('Cl convergence')   
#plt.savefig("1.1.jpg", bbox_inches='tight')

fig2, ax3 = plt.subplots()
ax4 = ax3.twinx()
ax3.set_ylabel('$ Cm_{le} $')
ax3.set_xlabel('N')
ax4.set_ylabel('Error (%)')
plt.xlabel("N")
ax3.plot(Mlist,resultCm, color='Blue',label='$ Cm_{le} $')
ax4.plot(Mlist,errorCm, color='Red', label='Error')
ax3.legend(loc = 'best', bbox_to_anchor=(0.15, 0.105,0.855,0.6), frameon = False)
ax4.legend(loc = 'best', bbox_to_anchor=(0.15, 0.1,0.855,0.5), frameon = False)
plt.title('$ Cm_{le} $ convergence')   
#plt.savefig("1.2.jpg", bbox_inches='tight')


# In[4]:


#      2. Validation assessment 
# 2.1. Experimental results

alpha_0_exp = -2.24
cla_exp = 6.01 *m.pi/180
cm_0_exp = -0.051

# 2.2. Calculated results
alpha_0 = find_alpha0(-10,10,1,0)[0]
cla = find_cl_alpha(1,0)
cm_0 = find_cm_0(alpha_0,1,0)

error_alpha_0 = abs((alpha_0_exp - alpha_0) * 100 / alpha_0_exp)
error_cla = abs((cla_exp - cla) * 100 / cla_exp)
errorcm_0 = abs((cm_0_exp - cm_0) * 100 / cm_0_exp)


# In[14]:


print('Calculated:   alpha_0 = ', alpha_0,'º,       cl_alpha =',cla,'º^-1,     Cm_0 =',cm_0)
print('Exprimental:  alpha_0 = ', alpha_0_exp,'º,                cl_alpha =',cla_exp,'º^-1,     Cm_0 =',cm_0_exp )
print('Errors:    E(alpha_0) =  ', error_alpha_0, '%    E(cl_alpha) = ', error_cla, '%       E(Cm_0) = ',errorcm_0,'%')


# In[16]:


# 2.3. Flap efficiency
# 2.3.1 Efficiency, function of deflection of the flap
xhlist = list([0.85, 0.8, 0.75, 0.7])
eta_increment = 0.5
etalist = list(np.arange(0,10 + eta_increment,eta_increment))

alpha_0_matr = np.zeros((len(xhlist),len(etalist)))

for xh in xhlist:
    print(xh)
    for eta in etalist:     
        alpha_0_matr[xhlist.index(xh)][etalist.index(eta)] = find_alpha0(-40,10,xh,eta)[0]


# In[17]:


print(alpha_0_matr)


# In[18]:


print(len(etalist))


# In[19]:


result_085 = []
result_080 = []
result_075 = []
result_070 = []


for i in range(0,len(etalist)-1):
    result_085.append(abs(alpha_0_matr[0][i+1] - alpha_0_matr[0][i]) / (eta_increment))
    result_080.append(abs(alpha_0_matr[1][i+1] - alpha_0_matr[1][i]) / (eta_increment))
    result_075.append(abs(alpha_0_matr[2][i+1] - alpha_0_matr[2][i]) / (eta_increment))
    result_070.append(abs(alpha_0_matr[3][i+1] - alpha_0_matr[3][i]) / (eta_increment))


# In[20]:



etalist.pop()
plt.title('Flap efficiency vs flap deflection angle')
plt.xlabel('$\\eta$ (º)', fontsize = 10)
plt.ylabel('$\\frac{\\Delta \\alpha_{l,0}}{\\Delta \\eta}$', fontsize = 15)
plt.plot(etalist,result_085,color='red',label = 'E = 0.15')
plt.plot(etalist,result_080, color ='orange', label = 'E = 0.20')
plt.plot(etalist,result_075, color ='green', label = 'E = 0.25')
plt.plot(etalist,result_070, color = 'blue', label = 'E = 0.30')
plt.legend()
plt.savefig("2.3.1.jpg", bbox_inches='tight')


# In[21]:


# 2.3.2 Flap efficiency, function of flap-chord ratio.
E_list = [0.15, 0.20, 0.25, 0.30]
eff_list = []
thetahlist = []

for E in E_list:
    xhinge = 1 - E
    alpha_0 = find_alpha0(-10,10,xhinge,0)[0]
    alpha_0_plus = find_alpha0(-10,10,xhinge,10)[0]
    eff_list.append(abs(alpha_0_plus - alpha_0)/10)
    thetahlist.append(m.acos(1-2*xhinge))
    
theoretical_efficiency = []
for th in thetahlist:
     theoretical_efficiency. append(-(1-th/m.pi + m.sin(th)/m.pi))
print(theoretical_efficiency)



# In[ ]:


factor = 0
for i in range(0,len(eff_list)):
    f = eff_list[i]/theoretical_efficiency[i]
    factor = factor + f
    print(f)
print(0.25*factor)


# In[28]:


factor = 0.6
eff_list_2 = []
for i in range(0,len(eff_list)):
    eff_list_2.append( eff_list[i] * factor)
    
plt.title('Flap efficiency vs flap-chord ratio')
plt.xlabel('E', fontsize = 10)
plt.ylabel('$\\frac{\\Delta \\alpha_{l,0}}{\\Delta \\eta}$', fontsize = 15)
plt.plot(E_list, eff_list, '*-', color = 'blue',label = 'Calculated ')
plt.plot(E_list, eff_list_2, '*-', color = 'red',label = 'Corrected factor = 0.6 ')
plt.xlim(0.10, 0.35)
plt.ylim(0, 0.7)
plt.legend()
plt.savefig("2.3.2.jpg", bbox_inches='tight')


# In[28]:


# thetahlist = []
# for xh in xhlist:
#     thetahlist.append(m.acos(1-2*xh))
# theoretical_efficiency = []

# for th in thetahlist:
#     theoretical_efficiency. append(-(1-th/m.pi + m.sin(th)/m.pi))


# correction_085 = []
# correction_080 = []
# correction_075 = []
# correction_070 = []

# for i in range(0,len(result_085)):
#     correction_085.append(result_085[i]/theoretical_efficiency[0])
#     correction_080.append(result_080[i]/theoretical_efficiency[1])
#     correction_075.append(result_075[i]/theoretical_efficiency[2])
#     correction_070.append(result_070[i]/theoretical_efficiency[3])


# In[29]:


# plt.plot(etalist, correction_085, color ='red')
# plt.plot(etalist, correction_080, color ='yellow')
# plt.plot(etalist, correction_075, color ='green')
# plt.plot(etalist, correction_070, color ='blue')


# In[3]:


# 3.1. Discussion
f_list = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06]
p_list = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6]
alpha_0_matr_fp = np.zeros((len(f_list),len(p_list)))
cm_0_matr_fp = np.zeros((len(f_list),len(p_list)))

for f in f_list:
    print(f)
    for p in p_list:
        alpha_0_fp = find_alpha0(-15,0,1,0)[0]
        alpha_0_matr_fp[f_list.index(f)][p_list.index(p)] = alpha_0_fp
        cm_0_matr_fp[f_list.index(f)][p_list.index(p)] = find_cm_0(0,1,0) 


# In[4]:


alpha_0_f1 = []
alpha_0_f2 = []
alpha_0_f3 = []
alpha_0_f4 = []
alpha_0_f5 = []
alpha_0_f6 = []
alpha_0_f7 = []
for i in range(0,len(p_list)):
    alpha_0_f1.append(alpha_0_matr_fp[0][i])
    alpha_0_f2.append(alpha_0_matr_fp[1][i])
    alpha_0_f3.append(alpha_0_matr_fp[2][i])
    alpha_0_f4.append(alpha_0_matr_fp[3][i])
    alpha_0_f5.append(alpha_0_matr_fp[4][i])
    alpha_0_f6.append(alpha_0_matr_fp[5][i])
    alpha_0_f7.append(alpha_0_matr_fp[6][i])
    
    
plt.title('$\\alpha_{l,0}$ dependence with p')
plt.xlabel('p', fontsize = 10)
plt.ylabel('$\\alpha_{l,0}(º)$', fontsize = 15)
plt.plot(p_list, alpha_0_f1 ,label = 'f = 0.00')
plt.plot(p_list, alpha_0_f2 ,label = 'f = 0.01')
plt.plot(p_list, alpha_0_f3 ,label = 'f = 0.02')
plt.plot(p_list, alpha_0_f4 ,label = 'f = 0.03')
plt.plot(p_list, alpha_0_f5 ,label = 'f = 0.04')
plt.plot(p_list, alpha_0_f6 ,label = 'f = 0.05')
plt.plot(p_list, alpha_0_f7 ,label = 'f = 0.06')
plt.legend()
plt.savefig("3.1.jpg", bbox_inches='tight')


# In[5]:


cm_0_f1 = []
cm_0_f2 = []
cm_0_f3 = []
cm_0_f4 = []
cm_0_f5 = []
cm_0_f6 = []
cm_0_f7 = []
for i in range(0,len(p_list)):
    cm_0_f1.append(cm_0_matr_fp[0][i])
    cm_0_f2.append(cm_0_matr_fp[1][i])
    cm_0_f3.append(cm_0_matr_fp[2][i])
    cm_0_f4.append(cm_0_matr_fp[3][i])
    cm_0_f5.append(cm_0_matr_fp[4][i])
    cm_0_f6.append(cm_0_matr_fp[5][i])
    cm_0_f7.append(cm_0_matr_fp[6][i])
    
plt.title('$cm_{0}$ dependence with p')
plt.xlabel('p', fontsize = 10)
plt.ylabel('$cm_{0}$', fontsize = 15)
plt.plot(p_list, cm_0_f1 ,label = 'f = 0.00')
plt.plot(p_list, cm_0_f2 ,label = 'f = 0.01')
plt.plot(p_list, cm_0_f3 ,label = 'f = 0.02')
plt.plot(p_list, cm_0_f4 ,label = 'f = 0.03')
plt.plot(p_list, cm_0_f5 ,label = 'f = 0.04')
plt.plot(p_list, cm_0_f6 ,label = 'f = 0.05')
plt.plot(p_list, cm_0_f7 ,label = 'f = 0.06')
plt.legend()
plt.savefig("3.2.jpg", bbox_inches='tight')


# In[ ]:


# NOTES:
# - Amb flap, canvia alpha_0, però no cl_alpha


# In[1]:


print(cm_0_matr_fp)

