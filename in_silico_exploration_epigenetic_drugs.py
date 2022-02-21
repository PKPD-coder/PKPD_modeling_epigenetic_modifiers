# -*- coding: utf-8 -*-
"""
Created on Mon Jan  3 16:52:56 2022

@author: vandevya
"""
"""

NOTE TO NEW USERS: the file is set-up so that simply running the full code (by pressing F5) 
will generate all figures from the article sequentially. This may not be ideal if you want 
to study the effects of custom dosing regimens or different parameter values.
The first three code blocks should in principle not be changed by the user.
To run a code block, press Ctrl + Enter

"""

#%% import all required packages
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.interpolate import interp1d
from sklearn.utils.extmath import cartesian
import pandas as pd
import math
import matplotlib.patches as patches
import seaborn as sns
from PIL import Image
from matplotlib import rc
#%% Defining the system of coupled ODEs

hmax =1 #max allowed stepsize in ode integration: should in principle not be change except in case of stability issues (in which case the problem is probably caused by poor choice of parameters or initials)
endtime = 1700 #simulation stop time
steps =int(endtime/hmax + 1)
t = np.linspace(0,endtime,steps) #time array


#Defining system of ODEs
def EpigeneticModulator_repeat(SV,t,param,Dose, inputs):
    Depot = SV[0] #umol
    Central = SV[1] #umol/L
    Peripheral= SV[2]
    LSD1b = SV[3] 

    GRP = SV[4]
    P = SV[5] 
    Q = SV[6] 
 
    Prol = SV[7] #cells/L
    Transit1 = SV[8] #cells/L
    Transit2 = SV[9] #cells/L
    Transit3 = SV[10] #cells/L
    Circ = SV[11] #cells/L
    P_control = SV[12]
    
    Exposure = SV[13]
    CircCtrl = SV[14]
    Circ_AUC = SV[15]
    CircCtrl_AUC =SV[16]
    
    [LSD1_0,kP,k50P, kmaxPQ, nPQ, kmaxQP,nQP,bs, k50QP,b,kdeg,Vm,Ki, kinact,n,k50,k50_BM, Km,ka,Vc, kel, Qc,kmax,kprol,g,ktr] = param
    TE = LSD1b/LSD1_0*100
    LSD1u = LSD1_0 - LSD1b
    vP =kP*(1-P/(k50P+P))*P
    perGRP = 1 - GRP

    BM = perGRP
    BMPQ = BM

    vPQ = kmaxPQ*(((BMPQ**nPQ)/((k50PQ**nPQ) + (BMPQ**nPQ))))*P
    BMQP = 1 - BM

    vQP = kmaxQP*(((BMQP**nQP)/((k50QP**nQP)+(BMQP**nQP)))+bs)*Q

    m = b -kdeg
    kmax = -m*GRP + b


    dydt = [(Dose*inputs(t)-ka*Depot),#
            ka*Depot - CL/Vc*Central - Qc/Vc*Central + Qc/Vp*Peripheral,
            Qc/Vc*Central - Qc/Vp*Peripheral,
            
            -LSD1b*(Vm/(Km+LSD1b))+ ((kinact*Central/Vc)/(Ki+Central/Vc))*LSD1u,
            kmax*(k50**n/(k50**n+TE**n))-kdeg*GRP,
            vP - vPQ + vQP,
            vPQ - vQP,

            kprol*Prol*(1-Emax*(Central/Vc**n)/(k50_BM**n+Central/Vc**n))*(Circ0/Circ)**g - ktr*Prol,
            ktr*(Prol-Transit1),
            ktr*(Transit1 - Transit2),
            ktr*(Transit2 - Transit3),
            ktr*(Transit3-Circ),
    
            kP*(1-P_control/(k50P+P_control))*P_control,
            Central,
            0,
            Circ,
            CircCtrl]
            
    return dydt

#%% Define function to generate custom dosing regimens

def regimens(total_days, on_dosing, off_dosing,Dose_mgkg, param):
    hmax=1
    total_days = total_days*24
    days_with_dose = on_dosing
    on_days = days_with_dose


    off_days=off_dosing*24
    with_dosing_days = on_days
    no_dosing_days = off_days

    on = []
    on.append((np.arange(0+j*(off_days+days_with_dose*24),days_with_dose*24+j*(off_days+days_with_dose*24),24) for j in range(math.ceil(total_days/((with_dosing_days*24)+no_dosing_days)))))#math.ceil(total_days/(off_days+on_days)))))
    on_days = [element for lis in on for element in lis]
    on_days= np.concatenate(on_days)
    dosing_days = [i for i in on_days]# if i < 100*24]#not in np.concatenate((np.arange(5,27), np.arange(33,50),np.arange(56,80)))] #list of all dosing days, changed units to hours (same unit as time t)

    t = np.linspace(0,total_days,int(total_days/hmax)+1)
    dosing_events = np.zeros(len(t))
    number_doses = len(dosing_days)
    for i in range(len(t)):
        for day in dosing_days:
            if int(t[i]) == int(day):
                dosing_events[i] = 1           
    inputs = interp1d(t,dosing_events,bounds_error=False,kind='previous',fill_value=0)
    Dosing = BW*((Dose_mgkg*10**-3)/MW)*10**6/number_doses
    simulation = odeint(EpigeneticModulator_repeat, SV0, t, args =(param,Dosing,inputs,),full_output=0,hmax=hmax,rtol=1e-10)


    Depot = simulation[:,0] #umol
    Central = simulation[:,1] #umol/L
    Peripheral= simulation[:,2]
    LSD1b = simulation[:,3] 

    GRP = simulation[:,4] #cells/L
    P = simulation[:,5] 
    Q = simulation[:,6] 

    Prol = simulation[:,7] #cells/L
    Transit1 = simulation[:,8] #cells/L
    Transit2 = simulation[:,9] #cells/L
    Transit3 = simulation[:,10] #cells/L
    Circ = simulation[:,11] #cells/L
    P_control = simulation[:,12]

    TE = LSD1b/LSD1_0*100
    
    Tumor_vol = P+Q
    Tumor_vol_ctrl = P_control

    return Depot, Central, Peripheral, LSD1b, GRP, P, Q, Prol, Transit1,Transit2, Transit3, Circ, P_control,TE, Tumor_vol, Tumor_vol_ctrl,t # P_AUC, P_controlAUC, Exposure,



#%% Define parameter values and initial conditions
"""
------------------------------------------------------------
PARAMETERS
------------------------------------------------------------
"""
    
BW = 0.025 #kg

ktr= 4/120
kprol = ktr
Circ0 = 300*10**9 #platelets/L
g = 0.145
ka =0.88205 #1/h


Vc = 0.57749 #L
kel=0.177 #1/h
CL=kel*Vc
MW=200 #Dalton
Qc =CL*2
Vp = 102*BW
Emax = 0.4

kmax =1 #unitless


LSD1_0 = 3.0833 #nM
Ki =0.5 #uM
kinact = 0.87052
Vm = 0.3515 #nM/h
Km = 0.9028 #nM
k50 = 50
n = 2
kdeg = 0.025567
b = 0.083048
kP = 0.0037
k50P = 100000
kmaxPQ = 3
kmaxQP = 4.4908
bs= 0.033328 #Bias towards Q-t-P transition
k50PQ = 0.8836
k50QP = 0.99897
nPQ = 3
nQP = 37
k50_BM  = 0.5


Dose_mgkg = 800
Dose = BW*((Dose_mgkg*10**-3)/MW)*10**6 #total dose in mg/kg converted to umol
#Dose = Dose/number_doses #total dose distributed over number of total doses

#putting all parameters in list param (required as input for functions 'regimens()' and 'EpigeneticModulator_repeat()')
param = [LSD1_0,kP,k50P, kmaxPQ, nPQ, kmaxQP,nQP,bs, k50QP,b,kdeg,Vm,Ki, kinact,n,k50,k50_BM, Km,ka,Vc, kel, Qc,kmax,kprol,g,ktr]

"""
----------------------------------------------------------
Initial Conditions
----------------------------------------------------------
"""
Depot_0 = 0
Central_0 = 0
Peripheral_0 = 0
LSD1b_0 =0 #initial concentration of bound LSD1
GRP_0 = 1 #initial GRP relative activity
P_0 = 70 #initial proliferating tumor volume 
Q_0 = 0 #initial quiescent tumor volume 
Circ0 = 300*10**9 #initial circulating platet count

SV0 = [Depot_0,Central_0,Peripheral_0,LSD1b_0, GRP_0,P_0,Q_0, Circ0,Circ0,Circ0,Circ0,Circ0,P_0,0,Circ0,0,0] #to ensure equilibrium, either simulate for +1000h before adding drug or put all cell conc equal to Circ0 including transit compartments



#%% --FIGURE 1C-- 
    
"""Fill in"""
total_days=30
off_dosing=0
on_dosing =1
Dose_mgkg =50
"""------"""

Depot, Central, Peripheral, LSD1b, GRP, P, Q, Prol, Transit1,Transit2, Transit3, Circ, P_control, TE, Tumor_vol, Tumor_vol_ctrl,t = regimens(total_days, on_dosing, off_dosing,Dose_mgkg, param) #P_AUC, P_controlAUC, Exposure,

figure1b = plt.figure(figsize=(12,2))
fig1BPK = plt.subplot2grid((1,4),(0,0))
fig1BPK.plot(t/24,Central/Vc,lw=2.5)
fig1BPK.set_ylabel('PK')
fig1BPK.set_ylim(0,max(Central/Vc)*1.3)

fig1Bplatelets = plt.subplot2grid((1,4),(0,3))
fig1Bplatelets.plot(t/24,Circ,lw=2.5)
fig1Bplatelets.set_ylabel('Circulating Platelets')
fig1Bplatelets.set_ylim(1e11,4e11)

fig1Bbiomarker = plt.subplot2grid((1,4),(0,1))
fig1Bbiomarker.plot(t/24, GRP,lw=2.5)
fig1Bbiomarker.set_ylabel('Biomarker Response')
fig1Bbiomarker.set_ylim(top=max(GRP)*1.1)

fig1Btumor = plt.subplot2grid((1,4),(0,2))
fig1Btumor.plot(t/24, P+Q,lw=2.5)
fig1Btumor.plot(t/24, P_control,lw=2.5, color='lightgrey',ls=':')
fig1Btumor.set_ylabel('Tumor Volume')
fig1Btumor.set_ylim(top=max(P_control)*1.1)

subs = [fig1BPK,fig1Bplatelets,fig1Bbiomarker,fig1Btumor]
for j in subs:
    j.yaxis.set_visible(True)
    plt.setp(j.get_yticklabels(),visible=False)
    j.set_yticks([])
    j.xaxis.set_ticks([])
    j.set_xlabel('Time')
    j.spines['right'].set_visible(False)
    j.spines['top'].set_visible(False)

    
#plt.savefig('figures\Figure 1B.jpeg',dpi=1000) #--> uncomment this line if you want to save the figure
#%% --FIGURE 2--


"""Fill in"""
total_days=30 #total treatment time in days
Reg1_off_dosing=6 #dose regimen 1: days without dose per on/off cycle
Reg2_off_dosing=20 #dose regimen 2: days without dose per on/off cycle
Reg1_on_dosing=1 #dose regimen 1: days with dose per on/off cycle
Reg2_on_dosing=1 #dose regimen 2: days with dose per on/off cycle
Low_Dose_mgkg = 100 #low dose in mg/kg
High_Dose_mgkg = 1000 #high dose in mg/kg
"""------"""


rc('text', usetex=False)
figure2=plt.figure(figsize=(15,15))

"""figure A will not be shown since this picture was generated with Biorender.com and put inside the original figure 2 as shown in the publication
-------------------------------------------------------------
A = plt.subplot2grid((12,15),(0,0), colspan=15,rowspan=7)
#fill in picture file path --> A.imshow(Image.open(r'C:\....'))
A.spines['right'].set_visible(False)
A.spines['left'].set_visible(False)
A.spines['top'].set_visible(False)
A.spines['bottom'].set_visible(False)
A.yaxis.set_visible(False)
A.xaxis.set_visible(False)
A.axis('tight')
A.axis('off')
#"""

C1 = plt.subplot2grid((12,15),(7,1), colspan=7,rowspan=4)
C2 = plt.subplot2grid((12,15),(7,8), colspan=7,rowspan=4)

Depot, Central, Peripheral, LSD1b, GRP, P, Q, Prol, Transit1,Transit2, Transit3, Circ, P_control,TE, Tumor_vol, Tumor_vol_ctrl,t = regimens(total_days, Reg2_on_dosing, Reg2_off_dosing,Low_Dose_mgkg, param) #P_AUC, P_controlAUC, Exposure,
ac=C1.plot(t/24,P+Q,color='dodgerblue',ls='--',lw=3, label = 'Dosing regimen 1: low dose')
C2.plot(t/24, Circ,color='dodgerblue',ls='--',lw=3)

Depot, Central, Peripheral, LSD1b, GRP, P, Q, Prol, Transit1,Transit2, Transit3, Circ, P_control, TE, Tumor_vol, Tumor_vol_ctrl,t = regimens(total_days, Reg2_on_dosing, Reg2_off_dosing,High_Dose_mgkg, param) #P_AUC, P_controlAUC, Exposure,
ad=C1.plot(t/24, P+Q,color='mediumblue',ls='--',lw=3, label = 'Dosing regimen 1: high dose')
C2.plot(t/24, Circ,color='mediumblue',ls='--',lw=3)

Depot, Central, Peripheral, LSD1b, GRP, P, Q, Prol, Transit1,Transit2, Transit3, Circ, P_control, TE, Tumor_vol, Tumor_vol_ctrl,t = regimens(total_days, Reg1_on_dosing, Reg1_off_dosing,Low_Dose_mgkg, param) #P_AUC, P_controlAUC, Exposure,
aa=C1.plot(t/24, P+Q, color='rosybrown',ls=':',lw=3, label = 'Dosing regimen 2: low dose')
C2.plot(t/24, Circ, color='rosybrown',ls=':',lw=3)

Depot, Central, Peripheral, LSD1b, GRP, P, Q, Prol, Transit1,Transit2, Transit3, Circ, P_control, TE, Tumor_vol, Tumor_vol_ctrl,t = regimens(total_days, Reg1_on_dosing, Reg1_off_dosing,High_Dose_mgkg, param) #P_AUC, P_controlAUC, Exposure,
ab=C1.plot(t/24, P+Q,color='firebrick',ls=':',lw=3, label = 'Dosing regimen 2: high dose')
C2.plot(t/24, Circ,color='firebrick',ls=':',lw=3)



C1.set_xlabel('Time [days]', fontsize=15)
C1.set_ylabel('Tumor Volume [$mm^3$]', fontsize=16)
C1.tick_params(axis='both',labelsize=15)


C2.set_xlabel('Time [days]', fontsize=15)
C2.set_ylabel('Platelet count [cells/L]', fontsize=16)
C2.tick_params(axis='both',labelsize=15)

figure2.text(0.02,0.97, 'A',size = 20, weight = 'bold')
figure2.text(0.02,0.42,'B',size = 20, weight = 'bold')
figure2.text(0.52,0.42, 'C',size = 20, weight = 'bold')


plt.tight_layout()

#the legend is treated as a separate 'plot' for ease of layout
all_plots = ac+ad+aa+ab
all_legends = [lab.get_label() for lab in all_plots]
C3 = plt.subplot2grid((12,15),(11,1),colspan =13)
C3.legend(all_plots,all_legends,ncol=2, bbox_to_anchor=(0.92,1), fontsize =20)
C3.axis('off')

plt.tight_layout()
#figure2.savefig('figures\Figure2.jpeg',dpi=1000) #--> uncomment this line if you want to save the figure

#%% --FIGURE 3--

"""Fill in"""
total_days=100 #total treatment time in days
Reg1_off_dosing=0 #dose regimen 1: days without dose per on/off cycle
Reg2_off_dosing=13 #dose regimen 2: days without dose per on/off cycle
Reg3_off_dosing=27 #dose regimen 3: days without dose per on/off cycle
Reg1_on_dosing=1 #dose regimen 1: days with dose per on/off cycle
Reg2_on_dosing=1 #dose regimen 2: days with dose per on/off cycle
Reg3_on_dosing=1 #dose regimen 3: days with dose per on/off cycle
Dose_mgkg = 1000 #dose in mg/kg
Dose_control = 0 #no dose for control group
"""------"""


figure3 = plt.figure(figsize=(10,7.5))
A = plt.subplot2grid((2,1),(0,0)) #plot for TGI
B = plt.subplot2grid((2,1),(1,0)) #plot for platelet counts

#Regimen 1
Depot, Central, Peripheral, LSD1b, GRP, P, Q, Prol, Transit1,Transit2, Transit3, Circ, P_control, TE, Tumor_vol, Tumor_vol_ctrl,t = regimens(total_days, Reg1_on_dosing, Reg1_off_dosing,Dose_mgkg, param) #P_AUC, P_controlAUC, Exposure,
A.plot(t/24,P+Q,label= 'QD',lw=3)
B.plot(t/24,Circ,label= 'QD' ,lw=3)

#Regimen 2
Depot, Central, Peripheral, LSD1b, GRP, P, Q, Prol, Transit1,Transit2, Transit3, Circ, P_control, TE, Tumor_vol, Tumor_vol_ctrl,t = regimens(total_days, Reg2_on_dosing, Reg2_off_dosing,Dose_mgkg, param) #P_AUC, P_controlAUC, Exposure,
A.plot(t/24,P+Q, label = 'Q2W',lw=3)
B.plot(t/24,Circ, label = 'Q2W',lw=3)

#Regimen 3
Depot, Central, Peripheral, LSD1b, GRP, P, Q, Prol, Transit1,Transit2, Transit3, Circ, P_control, TE, Tumor_vol, Tumor_vol_ctrl,t = regimens(total_days, Reg3_on_dosing, Reg3_off_dosing,Dose_mgkg, param)#P_AUC, P_controlAUC, Exposure,
A.plot(t/24,P+Q, label = 'Q4W',lw=3)
B.plot(t/24,Circ, label = 'Q4W',lw=3)

#Control
Depot, Central, Peripheral, LSD1b, GRP, P, Q, Prol, Transit1,Transit2, Transit3, Circ, P_control, TE, Tumor_vol, Tumor_vol_ctrl,t = regimens(total_days, Reg1_on_dosing, Reg1_off_dosing,Dose_control, param)#P_AUC, P_controlAUC, Exposure,
A.plot(t/24,P_control, label='Control',lw=3,ls=':')
B.plot(t/24,Circ, label='Control',lw=3,ls=':')

#figure styling
A.set_ylabel('tumor volume [mm$^3$]',fontsize=12,labelpad=None)
A.tick_params(axis='both',labelsize=12)
A.legend()
A.text(-10,140000,'A',fontsize=15, weight='bold')
B.set_ylabel('platelets [cells/L]',fontsize=12,labelpad=30)
B.set_xlabel('time [days]',fontsize=12)
B.tick_params(axis='both',labelsize=12)

#add colored rectangles depicting thrombocytopenia grade per CTCAE v5.0 hematologic criteria
rectG1 = patches.Rectangle((-300, 0.75e11), 3000, 0.75e11, linewidth=1, edgecolor='none', facecolor='darkseagreen', alpha=99)
rectG2 = patches.Rectangle((-300, 0.5e11), 3000, 0.25e11, linewidth=1, edgecolor='none', facecolor='moccasin', alpha=99)
rectG3 = patches.Rectangle((-300, 0.25e11), 3000, 0.25e11, linewidth=1, edgecolor='none', facecolor='coral', alpha=99)
rectG4 = patches.Rectangle((-300, 0e11), 3000, 0.25e11, linewidth=1, edgecolor='none', facecolor='pink', alpha=99)

gradess = [rectG1,rectG2,rectG3,rectG4]
for a in gradess:
    B.add_patch(a)
B.text(-10,3.5e11,'B',fontsize=15, weight='bold')

#figure3.savefig(r'figures\figure3.jpeg',dpi=1000) #--> uncomment this line if you want to save the figure

#%% Guide on generating custom dosing schemes
'''
Dosing regimens should be interpreted as follows. 1xON means 1 day with administered dose. When this is followed by 0xOFF, this means there are no days without dose, so 1xON, 0xOFF refers to once daily dosing. 
Custom dosing schemes can be tested when this notion is followed.
As of current, the code does not allow to simulate multiple doses per day or more irregular dosing schemes.

Examples:
1) QD (1xON, 0xOFF) -> once daily dosing
2) QW (1xON, 6xOFF) -> once weekly dosing
3) weekdays (5xON, 2xOFF) -> dosing on weekdays, no dosing over weekend
4) 3 days (3xON, 4xOFF) -> dosing for three days, followed by four days without dose
5) Q2W (1xON, 13xOFF) -> once dosing every second week
6) Q3W (1xON, 20xOFF) -> once dosing every third week
7) Q4W (1xON, 27xOFF) -> once dosing every fourth week
    
The following section takes multiple dosing regimens to be simulated. The following lines:
days_with_dose = [1,1,5]
days_without_dose = [0,6,2]
will allow the code to simulate the regimens [(1xON, 0xOFF),(1xON,6xOFF),(2xON, 5xOFF)]
or in other words: [daily,weekly,weekdays]
the line:
Dose_range = [10,100,1000] will sequentially simulate the aformentioned regimens at the specified doses
'''

#%% --Simulation of multiple doses and regimens (needed for figure 4 and 5)-----


"""fill in"""
##################################################################

total_days = 100 #total days on treatment (i.e., ON days + OFF days)
endtime = total_days #total simulated days: can be longer than treatment time (e.g., to study recovery of platelets after permanent treatment cessation)
hmax =1 #max allowed stepsize in ODE integration: should in principle not be change except in case of stability issues (in which case the problem is probably caused by poor choice of parameters or initials)
steps =int(endtime/hmax + 1) #defining stepsize of time array
t = np.linspace(0,endtime,steps)

#define sequential ON-dosing and OFF-dosing days:
days_with_dose = [1,1,5,3,1,1,1]
days_without_dose = [0,6,2,4,13,20,27]
#define list of total doses to be simulated
Dose_range =[10,30,100,300,1000,3000,10000]

####################################################################
""""""

#initialisation of lists that will be appended with the simulation results from each run
TGI = []
TGIctrl= []
Platelets = []
Total_exposure = []
Platelets_ctrl=[]
PlateletsAUC=[]
Platelets_ctrlAUC=[]
Dosing_limit=[]


ons = np.tile(days_with_dose,len(Dose_range)) #list of on_days extended: to give it same lenght as cartesian product 'combinations'
count = 0
count_list =[]
combinations = cartesian((Dose_range,days_without_dose))
for combos,i in zip(combinations,ons):
    
    Dose_mgkg = combos[0]
    Dose = BW*((Dose_mgkg*10**-3)/MW)*10**6
    
    off_days = combos[1]*24
    off_dosing_days = off_days
    on_dosing_days=i*24
    number_doses = on_dosing_days*math.ceil(total_days*24/(on_dosing_days+off_dosing_days))
    on = []
    on.append((np.arange(0+j*(off_days+i*24),i*24+j*(off_days+i*24),24) for j in range(math.ceil(total_days*24/(on_dosing_days+off_dosing_days)))))
    on_days = [element for lis in on for element in lis]
    on_days= np.concatenate(on_days)
    dosing_days = [i for i in on_days if i < 100*24]
    
    t = np.linspace(0,total_days*24,total_days*24+1)
    dosing_events = np.zeros(len(t))
    number_doses = len(dosing_days)
    for i in range(len(t)):
        for day in dosing_days:
            if t[i] == day:
                dosing_events[i] = 1           
    inputs = interp1d(t,dosing_events,bounds_error=False,kind='previous',fill_value=0)
    
    #the following snippet of code is needed to check later whether platelets were sufficiently recovered when the new dose was given (required for figure 5)
    force_break = []
    force_break.append((np.arange(0+j*(off_days+i*24),i*24+j*(off_days+i*24),24) for j in range(math.ceil(total_days/(on_dosing_days+off_dosing_days)))))
    force_break_day = [element for lis in force_break for element in lis]
    force_break_day= np.concatenate(force_break_day)
    breaking_days = [i for i in force_break_day if i < 100*24]

    break_checks = np.zeros(len(t))
    number_break_checks = len(breaking_days)
    for i in range(len(t)):
        for day in breaking_days:
            if t[i] == day:
                break_checks[i] = 1           
    first_on_day_check = interp1d(t,break_checks,bounds_error=False,kind='previous',fill_value=0)
    
    
    
    Dose = Dose/number_doses
    param = [LSD1_0,kP,k50P, kmaxPQ, nPQ, kmaxQP,nQP, bs,k50QP,b,kdeg,Vm,Ki, kinact,n,k50,k50_BM, Km,ka,Vc, kel, Qc,kmax,kprol,g,ktr]
    #call ODE model and store results in 'simulation_interval
    simulation_interval = odeint(EpigeneticModulator_repeat, SV0,t, args=(param,Dose,inputs,), hmax=hmax)
    
    #appending the relevant results in the premade lists
    TGI.append(simulation_interval[:,5]+simulation_interval[:,6])
    TGIctrl.append(simulation_interval[:,12])
    Platelets.append(simulation_interval[:,11])

    count_list.append(number_doses)#record number of doses administered in each separate run
    Total_exposure.append(simulation_interval[:,13][-1])# to check that each regimen results in the same exposure by end of treatment
    Platelets_ctrl.append(simulation_interval[:,14][-1])
    PlateletsAUC.append(simulation_interval[:,15][-1])
    Platelets_ctrlAUC.append(simulation_interval[:,16][-1])
    Dosing_limit.append(any((list(map(lambda s: s !=0 and s < 1.0e11,inputs(t)*simulation_interval[:,11]))))) #the platelet count limit for allowing redosing may depend on expert opinion. In this simulation a limit of 100e9 cells/L was choses as presented by Warren W. Piette, Candace M. Broussard-Steinberg, in Comprehensive Dermatologic Drug Therapy (Fourth Edition), 2021) 

#%% Calculate usable metrics for figure 4 and 5: Percentage TGI, AUC ratio of platelets, and boolean for sufficient platelet recovery before next dose  

TGI_percentage=[]
for i in range(len(TGI)):
    TGI_percentage.append(1-(TGI[i]/TGIctrl[i])[-1])
    
#Thrombocytopenia grades according to CTCAE v5.0 hematologic criteria
Thrombocytopenia_G1 = 150*10**9
Thrombocytopenia_G2 = 75*10**9
Thrombocytopenia_G3 = 50*10**9
Thrombocytopenia_G4 = 25*10**9

Weight = []
for i in range(len(Platelets)):
    Safety_weight = 0
    if Thrombocytopenia_G2 <= min(Platelets[i]) < Thrombocytopenia_G1 :
        Safety_weight = 1
    elif Thrombocytopenia_G3 <= min(Platelets[i]) < Thrombocytopenia_G2:
        Safety_weight = 2
    elif Thrombocytopenia_G4 <= min(Platelets[i]) < Thrombocytopenia_G3:
        Safety_weight = 3
    elif min(Platelets[i]) < Thrombocytopenia_G4:
        Safety_weight = 4
    Weight.append(Safety_weight)
    
    
Platelets_AUCratio = np.array(PlateletsAUC)/np.array(Platelets_ctrlAUC)
Sufficient_recovery = np.array([not x for x in Dosing_limit])#returns 0 (FALSE) when next dosing attempt happens when platelets are still low

#%% Dataframe with simulation outputs

v=list(combinations[:,0])
for i in range(len(v)):
    v[i] = str(v[i]) + " mg/kg"
 
#give more recognizable names to regimens    
reg = list(combinations[:,1])
for index, element in enumerate(reg):
    if element == 0:
        reg[index] = 'QD' 
    elif element == 4:
        reg[index] = '3-ON/4-OFF'
    elif element == 2:
        reg[index] = 'weekdays [5xON/2xOFF]'
    elif element == 6:
        reg[index] = 'QW [1xON/6xOFF]'
    elif element == 13:
        reg[index] = 'Q2W [1xON/13xOFF]'
    elif element == 20:
        reg[index] = 'Q3W [1xON/20xOFF]'
    elif reg[index] == 27:
        reg[index] = 'Q4W [1xON/27xOFF]'

#invert list of booleans so that those regimens with premature dosing are designated True
Premature_dosing = [not element for element in Sufficient_recovery]
    
ClinEff = pd.DataFrame(list(zip(v,reg,100*np.array(TGI_percentage),100*np.array(Platelets_AUCratio), Weight,Premature_dosing)),
                       columns=['Total dose','Dosing regimen','Total TGI [%]','Platelets from baseline [%]', 'Grade Thrombocytopenia','Premature_dosing'])




#%% --FIGURE 4--
figure4=plt.figure(figsize=(12,8))

#step 1: scatter plot of the simulated doses and regimens with respect to the level of TGI achieved as well as the exposure to platelet count suppression
a=sns.scatterplot(data=ClinEff, y="Platelets from baseline [%]", x='Total TGI [%]',hue='Total dose', style='Dosing regimen',s=400, 
            palette= 'rainbow')

#step 2: connect identical regimens with lines
for regimen in ClinEff['Dosing regimen'].unique():
    sns.lineplot( x=ClinEff['Total TGI [%]'].loc[ClinEff['Dosing regimen'] == regimen], y = ClinEff['Platelets from baseline [%]'].loc[ClinEff['Dosing regimen'] == regimen],ls='-', color='grey',lw=2)

#styling
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=12,markerscale=2.5)
plt.xlabel('Efficacy\n(Tumor Growth Inhibition [%])',fontsize=15)
plt.xticks(fontsize=14)
plt.ylabel('Toxicity\n(Platelets from baseline [%])',fontsize=15, multialignment='center')
plt.yticks(fontsize=14)
plt.tight_layout()
#figure4.savefig(r'figures\figure4.jpeg',dpi=1000,bbox_inches='tight') #--> uncomment this line if you want to save the figure

#%% --FIGURE 5--
figure5=plt.figure()
figure5, axs = plt.subplots(ncols=2,figsize=(12,6))

#A: swarm plot of whether there is risk for premature dosing in function of platelet levels
a=sns.swarmplot(data=ClinEff, y="Platelets from baseline [%]",x="Premature_dosing",s=9, palette='flare',hue='Dosing regimen',ax=axs[0])

a.set_ylabel("Platelets from baseline [%]",fontsize=15)
a.set_xlabel("Premature Dosing",fontsize=15)
a.tick_params(axis='x', labelsize=14)
a.legend([],[], frameon=False)
a.text(-0.65,95,'A',fontsize=15,weight='bold')

#B: swarm plot of grades of thrombocytopenia in function of platelet levels
b=sns.swarmplot(data=ClinEff, y="Platelets from baseline [%]",x="Grade Thrombocytopenia",s=9, palette='flare',hue='Dosing regimen', ax=axs[1])

b.set_ylabel("Platelets from baseline [%]",fontsize=15)
b.set_xlabel("Thrombocytopenia Grades",fontsize=15)
b.tick_params(axis='x', labelsize=14)

b.text(-0.85,95,'B',fontsize=15,weight='bold')
b.legend(bbox_to_anchor=(0.88, -0.15), loc=0, borderaxespad=1.5,ncol=3,fontsize=12)
figure5.subplots_adjust(bottom=0.3)
#figure5.savefig('figures\Figure5.jpeg',dpi=1000) #--> uncomment this line if you want to save the figure

#%%  --SUPPLEMENTARY FIGURE S1--

"""Fill in"""
total_days=100
off_dosing=13
on_dosing =1
Dose_mgkg =1000
Dose_control = 0
"""------"""


figureS1 = plt.figure(figsize=(9,6))
A1 = plt.subplot2grid((2,2),(0,0),colspan=2)
A2 = plt.subplot2grid((2,2),(1,0),colspan=2)


Depot, Central, Peripheral, LSD1b, GRP, P, Q, Prol, Transit1,Transit2, Transit3, Circ, P_control,TE, Tumor_vol, Tumor_vol_ctrl,t = regimens(total_days, on_dosing, off_dosing,Dose_mgkg, param)
A1.plot(t/24,Circ,color='darkgreen',label='Q2W dosing')
A1.fill_between(t/24,Circ, color='darkgreen', alpha=0.5)
A2.plot(t/24,P+Q,color='darkblue', label= 'Q2W dosing')
A2.fill_between(t/24,P+Q, color='darkblue', alpha=0.99)

Depot, Central, Peripheral, LSD1b, GRP, P, Q, Prol, Transit1,Transit2, Transit3, Circ, P_control, TE, Tumor_vol, Tumor_vol_ctrl,t = regimens(total_days, on_dosing, off_dosing,Dose_control, param)
A1.plot(t/24,Circ,color='pink',label='control')
A1.fill_between(t/24,Circ, color='pink', alpha=0.5)
A2.plot(t/24,P_control,color='orange',label='Control')
A2.fill_between(t/24,P_control, color= 'orange', alpha=0.5)

#Figure styling
A1.set_ylabel('Platelet count [cells/L]', fontsize=12, labelpad=45)
A1.tick_params(axis='both',labelsize=12)
A1.legend(title='Circulating platelets', loc= 'upper center')

A2.set_ylabel('Tumor volume [$mm^3$]', fontsize=12)
A2.set_xlabel('Time [days]', fontsize=12)
A2.tick_params(axis='both',labelsize=12)
A2.legend(title='Tumor volume', loc= 'upper center')


#figureS1.savefig(r'figures\supp_figure_S1.jpeg',dpi=1000) #--> uncomment this line if you want to save the figure
