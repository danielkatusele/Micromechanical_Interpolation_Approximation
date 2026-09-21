import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'Times'
plt.rcParams["font.size"] =16.0
plt.rcParams["xtick.top"] = True 
plt.rcParams["ytick.right"] = True
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams["legend.fancybox"] = False
plt.rcParams["legend.framealpha"] = 1
plt.margins(0.0)


####################################################
    #############  SHEAR LOAD   ###############
####################################################


# Replace 'filename.xlsx' with the path to your Excel file
df_Numerical = pd.read_excel("Numerical_Ellipse_Shear45s_el_0.05_0.1.xlsx", engine="openpyxl")
df_Numerical2 = pd.read_excel("Numerical_Ellipse_Shear45s_el_0.1_0.1.xlsx", engine="openpyxl")
df_Numerical3 = pd.read_excel("Numerical_Ellipse_Shear45s_el_0.15_0.1.xlsx", engine="openpyxl")
#df_Numerical = pd.read_excel("Shear_LM_vol69.8.xlsx", engine="openpyxl")
#df_Numerical2 = pd.read_excel("Shear_LM_vol628.3.xlsx", engine="openpyxl")
#df_Numerical3 = pd.read_excel("Shear_LM_vol981.7.xlsx", engine="openpyxl")
df_Numerical4 = pd.read_excel("Shear_LM_vol3102.8.xlsx", engine="openpyxl")



# Display the first few rows
ShearLoadN = df_Numerical["Shear load"].to_numpy()
ShearLoadN2 = df_Numerical2["Shear load"].to_numpy()
ShearLoadN3 = df_Numerical3["Shear load"].to_numpy()
ShearLoadN4 = df_Numerical4["Shear load"].to_numpy()
NumericalS = df_Numerical["Numerical"].to_numpy()
NumericalS2 = df_Numerical2["Numerical"].to_numpy()
NumericalS3 = df_Numerical3["Numerical"].to_numpy()
NumericalS4 = df_Numerical4["Numerical"].to_numpy()

print(f"Max load/rot: S1= {ShearLoadN[-1]:.2f}/{NumericalS[-1]:.2f}")
print(f"Max load/rot: S2= {ShearLoadN2[-1]:.2f}/{NumericalS2[-1]:.2f}")
print(f"Max load/rot: S3= {ShearLoadN3[-1]:.2f}/{NumericalS3[-1]:.2f}")


def θapproxShear(a, b, λv, α):
    return np.arctan2(λv* (a**4 + 6* a**2 *b**2 + b**4 + (-a**4 + b**4)* np.cos(2*α)) ,
                    2* (a**4 + 6* a**2 *b**2 + b**4) + (a**4 - b**4)* λv* np.sin(2*α)) 


def θapproxtension(a, b, λv, α):
    return np.arctan2((a**4 - b**4)*(-(1/λv) + λv)* np.sin(2*α),\
        (a**4 + 6* a**2 *b**2 + b**4) *(1/λv + λv) + (a**4 -\
            b**4)* (-(1/λv) + λv)* np.cos(2*α))


def θpolardecomp(a, b, λv):
    return np.arctan2(-(np.sqrt(2 + λv**2 - λv* np.sqrt(4 + λv**2)) - np.sqrt(2 \
        + λv* (λv + np.sqrt(4 + λv**2))))/np.sqrt(2), ((λv + np.sqrt(4 + λv**2))*\
            np.sqrt(2 + λv**2 - λv* np.sqrt(4 + λv**2)) + (-λv + np.sqrt(4 + λv**2))*\
                np.sqrt(2 + λv* (λv + np.sqrt(4 + λv**2))))/(2* np.sqrt(2) *np.sqrt(4 + λv**2)))
    
    
def θpolardecomp(a, b, λv):
    return np.arctan2(-(np.sqrt(2 + λv**2 - λv* np.sqrt(4 + λv**2)) - np.sqrt(2 \
        + λv* (λv + np.sqrt(4 + λv**2))))/np.sqrt(2), ((λv + np.sqrt(4 + λv**2))*\
            np.sqrt(2 + λv**2 - λv* np.sqrt(4 + λv**2)) + (-λv + np.sqrt(4 + λv**2))*\
                np.sqrt(2 + λv* (λv + np.sqrt(4 + λv**2))))/(2* np.sqrt(2) *np.sqrt(4 + λv**2)))


#a_m = [0.1, 0.2, 0.4, 1.]
a_m = [0.1, 0.2, 0.4, 1.]
angles = [np.pi/4]
angles_l = ["0", "15^{\circ}", "30^{\circ}", "60^{\circ}", "90^{\circ}"]
types = ["-", "--", "-.", (0, (5, 10)), (0, (3,5,1,5,1,5))]
for i in range(len(angles)):
    λ = np.linspace(0, 3.5, 100)
    a = a_m[1]; b=0.1; agl= angles[i]

    M1 = 1.0; M2 = 0.
    θ = []
    θapp = []
    θPD = []
    θAff = []
    Diffθapp = []
    for τ in ShearLoadN:
        θapp.append(np.arctan2(M2 *np.cos(agl - θapproxShear(a, b, τ, agl)) + M1 *np.sin(agl - θapproxShear(a, b, τ, agl)) , 
                            M1* np.cos(agl - θapproxShear(a, b, τ, agl)) - M2 *np.sin(agl - θapproxShear(a, b, τ, agl))) - agl)
        
        θPD.append(np.arctan2(M2 *np.cos(agl - θpolardecomp(a, b, τ)) + M1 *np.sin(agl - θpolardecomp(a, b, τ)), 
                            M1* np.cos(agl - θpolardecomp(a, b, τ)) - M2 *np.sin(agl - θpolardecomp(a, b, τ))) - agl)
        
        θAff.append(np.arctan2(M2 *np.cos(agl) + M1*np.sin(agl), 
                            (M1 + M2 *τ) *np.cos(agl) + (-M2 + M1 *τ) *np.sin(agl)) - agl)
                                                                                        

    ar = int(a / b)
    #plt.plot(λ, θapp, linestyle=types[i], label=r"$\eta = {}$".format(ar), linewidth=3.)
    #plt.plot(λ, θapp, linestyle=types[i], label=r"$\alpha = {}$".format(angles_l[i]), linewidth=3.)


####################################################
    #############  TENSION LOAD   ###############

####################################################

## Replace 'filename.xlsx' with the path to your Excel file
# Replace 'filename.xlsx' with the path to your Excel file
df_NumericalT = pd.read_excel("Tension_Pen_vol_f_0.7_another.xlsx", engine="openpyxl")
df_NumericalT2 = pd.read_excel("Tension_Pen_vol_f_6.28_another_attempt.xlsx", engine="openpyxl")
df_NumericalT3 = pd.read_excel("Tension_Pen_vol_f_9.8_another_attempt.xlsx", engine="openpyxl")
df_NumericalT4 = pd.read_excel("Tension_Pen_vol_f_31_another.xlsx", engine="openpyxl")

# Display the first few rows
TensileLoadN = df_NumericalT["Tensile load"].to_numpy()
TensileLoadN2 = df_NumericalT2["Tensile load"].to_numpy()
TensileLoadN3 = df_NumericalT3["Tensile load"].to_numpy()
TensileLoadN4 = df_NumericalT4["Tensile load"].to_numpy()
#
## Display the first few rows
NumericalT = df_NumericalT["Numerical"].to_numpy()
NumericalT2 = df_NumericalT2["Numerical"].to_numpy()
NumericalT3 = df_NumericalT3["Numerical"].to_numpy()
NumericalT4 = df_NumericalT4["Numerical"].to_numpy()

#
λT = np.linspace(1., 2., 100)
a = 0.2; b=0.1; agl=np.pi/4
M1 = 1.0; M2 = 0.
#
θappT = []
θAffT = []
for τ in TensileLoadN:
    θappT.append(np.arctan2(M2 *np.cos(agl - θapproxtension(a, b, τ, agl)) + M1 *np.sin(agl - θapproxtension(a, b, τ, agl)) , 
                           M1* np.cos(agl - θapproxtension(a, b, τ, agl)) - M2 *np.sin(agl - θapproxtension(a, b, τ, agl))) - agl)
        
    θAffT.append(np.arctan2((M2* np.cos(agl) + M1* np.sin(agl))/τ, 
                                       τ* (M1* np.cos(agl) - M2* np.sin(agl)))- agl)
                              
θapp2 = []                             
for τ in TensileLoadN2:
    θapp2.append(np.arctan2(M2 *np.cos(agl - θapproxShear(a, b, τ, agl)) + M1 *np.sin(agl - θapproxShear(a, b, τ, agl)) , 
                        M1* np.cos(agl - θapproxShear(a, b, τ, agl)) - M2 *np.sin(agl - θapproxShear(a, b, τ, agl))) - agl)
   
θapp3 = []                             
for τ in TensileLoadN3:
    θapp3.append(np.arctan2(M2 *np.cos(agl - θapproxShear(a, b, τ, agl)) + M1 *np.sin(agl - θapproxShear(a, b, τ, agl)) , 
                        M1* np.cos(agl - θapproxShear(a, b, τ, agl)) - M2 *np.sin(agl - θapproxShear(a, b, τ, agl))) - agl)
    
θapp4 = []                             
for τ in TensileLoadN4:
    θapp4.append(np.arctan2(M2 *np.cos(agl - θapproxShear(a, b, τ, agl)) + M1 *np.sin(agl - θapproxShear(a, b, τ, agl)) , 
                        M1* np.cos(agl - θapproxShear(a, b, τ, agl)) - M2 *np.sin(agl - θapproxShear(a, b, τ, agl))) - agl)

θapp1 = np.arctan2(M2 *np.cos(agl - θapproxtension(a, b, TensileLoadN[-1], agl)) + M1 *np.sin(agl - θapproxtension(a, b, TensileLoadN[-1], agl)) , 
                   M1* np.cos(agl - θapproxtension(a, b, TensileLoadN[-1], agl)) - M2 *np.sin(agl - θapproxtension(a, b, TensileLoadN[-1], agl))) - agl

θapp2 = np.arctan2(M2 *np.cos(agl - θapproxtension(a, b, TensileLoadN2[-1], agl)) + M1 *np.sin(agl - θapproxtension(a, b, TensileLoadN2[-1], agl)) , 
                   M1* np.cos(agl - θapproxtension(a, b, TensileLoadN2[-1], agl)) - M2 *np.sin(agl - θapproxtension(a, b, TensileLoadN2[-1], agl))) - agl

θapp3 = np.arctan2(M2 *np.cos(agl - θapproxtension(a, b, TensileLoadN3[-1], agl)) + M1 *np.sin(agl - θapproxtension(a, b, TensileLoadN3[-1], agl)) , 
                   M1* np.cos(agl - θapproxtension(a, b, TensileLoadN3[-1], agl)) - M2 *np.sin(agl - θapproxtension(a, b, TensileLoadN3[-1], agl))) - agl

θapp4 = np.arctan2(M2 *np.cos(agl - θapproxtension(a, b, TensileLoadN4[-1], agl)) + M1 *np.sin(agl - θapproxtension(a, b, TensileLoadN4[-1], agl)) , 
                   M1* np.cos(agl - θapproxtension(a, b, TensileLoadN4[-1], agl)) - M2 *np.sin(agl - θapproxtension(a, b, TensileLoadN4[-1], agl))) - agl


Err1 = np.arctan2(M2 *np.cos(agl - θapproxShear(a, b, τ, agl)) + M1 *np.sin(agl - θapproxShear(a, b, τ, agl)) , 
                        M1* np.cos(agl - θapproxShear(a, b, τ, agl)) - M2 *np.sin(agl - θapproxShear(a, b, τ, agl))) - agl
                                                                                        
                                                                                        
print("N Shape: ", NumericalT2.shape, "N TensileLoad: ", TensileLoadN.shape, "N θapp: ", np.array(θappT).shape)
print("N Shape: ", NumericalT4.shape, "N TensileLoad: ", TensileLoadN4.shape, "N θapp: ", np.array(θapp2).shape)
rel_l2_ = np.linalg.norm(θappT  - NumericalT) / np.linalg.norm(θappT)
rel_l2_2 = np.linalg.norm(θapp2  - NumericalT2) / np.linalg.norm(θapp2)
rel_l2_3 = np.linalg.norm(θapp3  - NumericalT3) / np.linalg.norm(θapp3)
rel_l2_4 = np.linalg.norm(θapp4  - NumericalT4) / np.linalg.norm(θapp4)
print("Relative L2 error for shear load: ", rel_l2_*100, rel_l2_2*100, rel_l2_3*100, rel_l2_4*100)

print(f"Max load/rot: S1= {ShearLoadN[-1]:.2f}/{NumericalS[-1]:.2f}")
print(f"Max load/rot: S2= {ShearLoadN2[-1]:.2f}/{NumericalS2[-1]:.2f}")
print(f"Max load/rot: S3= {ShearLoadN3[-1]:.2f}/{NumericalS3[-1]:.2f}")
err1 = np.abs(NumericalS2[-1] - NumericalS[-1])/np.abs(NumericalS[-1])
err2 = np.abs(NumericalS3[-1]- NumericalS[-1])/np.abs(NumericalS[-1])
print(f"Relative L2 error element size (mid/fine): {err1:0.3f}")
print(f"Relative L2 error element size (max/fine): {err2:0.3f}")



plt.plot(ShearLoadN, -NumericalS3 ,"--", label=r"$3542$ Elements", linewidth=3, zorder=10, color="orange")
plt.plot(ShearLoadN2, -NumericalS2 ,"-.", label=r"$13508$ Elements", linewidth=3, zorder=10, color="blue")
plt.plot(ShearLoadN3, -NumericalS ,linestyle=(0, (3, 1, 1, 1, 1, 1)), label=r"$146970$ Elements", linewidth=3.5, zorder=9, color="green")
#plt.plot(ShearLoadN4, -NumericalS4 ,"--", label=r"$\phi = 31 \%$", linewidth=3., zorder=8, color="purple")
#plt.plot(ShearLoadN, -NumericalS ,"--", label=r"$\phi = 0.8 \%$", linewidth=3, zorder=10, color="orange")
#plt.plot(ShearLoadN2, -NumericalS2 ,"-.", label=r"$\phi = 6.3 \%$", linewidth=3, zorder=10, color="blue")
#plt.plot(ShearLoadN3, -NumericalS3 ,linestyle=(0, (3, 1, 1, 1, 1, 1)), label=r"$\phi = 10 \%$", linewidth=3.5, zorder=9, color="green")
#plt.plot(ShearLoadN4, -NumericalS4 ,"--", label=r"$\phi = 31 \%$", linewidth=3., zorder=8, color="purple")
plt.plot(ShearLoadN, θapp, label="MIA", linewidth=3., color ="red")
#plt.plot(ShearLoadN, θPD, label="PD", linewidth=3., color ="magenta")
#plt.plot(ShearLoadN, θAff,label="MM", linewidth=3., color ="black")
#plt.plot(TensileLoadN, NumericalT ,"--", label=r"$\phi =0.7 \%$", linewidth=3., zorder=7, color="green")
#plt.plot(TensileLoadN2, NumericalT2 ,"-.", label=r"$\phi =6.3 \%$", linewidth=3., zorder=7, color="orange")
#plt.plot(TensileLoadN3, NumericalT3 ,linestyle=(0, (3, 1, 1, 1, 1, 1)), label=r"$\phi =17.5 \%$", linewidth=3., zorder=7, color="blue")
#plt.plot(TensileLoadN4, NumericalT4 ,"--", label=r"$\phi =31 \%$", linewidth=3., zorder=7, color="purple")
#plt.plot(TensileLoadN, θappT, label="MIA", linewidth=3., color ="red", zorder=7)
#plt.plot(TensileLoadN, θAffT, label="MM", linestyle=(0, (5, 1)),linewidth=3., color ="black", zorder=7)

plt.xlabel(r"Shear load $\tau$")
plt.xticks([0, 1, 2.,3., 3.5], [r"$0$", r"$1$", r"$2$", r"$3$", r"$3.5$"])
plt.yticks([0, -np.pi/8, -np.pi/4, -3*np.pi/8, -np.pi/2], 
           ["0", r"$-22.5^{\circ}$", r"$-45^{\circ}$", r"$-67.5^{\circ}$", r"$-90^{\circ}$"])
#plt.yticks([0, -np.pi/16, -np.pi/8], 
#           ["0", r"$-11.25^{\circ}$", r"$-22.5^{\circ}$"])
##plt.plot(λT, θappT, label="MIA", linewidth=3., color ="blue")
#plt.plot(TensionLoadN, NumericalT,"--" , label="Numerical", linewidth=3., color ="red")
#plt.plot(λT, θAffT,"--", label="Material", linewidth=3., color ="green")
#plt.xlabel(r"Tensile load $\lambda$")
##
#plt.xticks([1, 1.5, 2.], [r"$1$", r"$1.5$", r"$2$"])
#plt.yticks([0, -np.pi/16,-np.pi/8, -3*np.pi/16], 
#           ["0", r"$-11.25^{\circ}$", r"$-22.5^{\circ}$", r"$-33.75^{\circ}$"])
plt.ylabel(r"Inclusion rotation $ \theta$")
#plt.xticks([0, 1, 2.,3, 4, 5, 6], [r"$0$", r"$1$", r"$2$", r"$3$", r"$4$", r"$5$", r"$6$"])
#plt.yticks([0, -np.pi/8, -np.pi/4, -3*np.pi/8, -np.pi/2], 
#           ["0", r"$-22.5^{\circ}$", r"$-45^{\circ}$", r"$-67.5^{\circ}$", r"$-90^{\circ}$"])
plt.text(2, -0.2, r"$\eta = 2$", fontsize=18)
plt.text(2, -0.3, r"$\alpha = 45^{\circ}$", fontsize=18)
plt.grid()

plt.legend()
plt.show()