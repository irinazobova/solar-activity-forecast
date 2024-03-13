import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import math
import scipy
from sympy import fwht, ifwht
from sympy.combinatorics import GrayCode
from scipy import stats


### importing the data: montly mean SSN from input.csv and dates of minima and maxima from timeframe-wiki.txt
data = pd.read_csv("/input.csv", sep=';', header = 0)
# plt.plot(data['yearfrac'],data['ssn']) # visualising SSN time series
f = open('/timeframe-wiki.txt')
s = f.read().splitlines()
f.close()
s = list(map(lambda x:x.split(), s))


### preprocessing the data: splitting the SSN time series into 24 cycles
sc_minima_indices = []
for i in range(len(s)):
  sc_minima_indices.append(i+1)
  sc_minima_indices.append(data[(data['year'] == int(s[i][1])) & (data['month'] == int(s[i][2]))].index[0])
list2 = [] # to split the cycle number and the start point of the cycle in the observation file input.csv
list2.append(sc_minima_indices[::2])
list2.append(sc_minima_indices[1::2])
cycles_data = dict() # creating dictionary for all cycles
cycles_binwidth = list() # rounded down amount of months in one phase for each cycle
cycles_width = list()
cycles_length = list()
for i in range(24):
  cycles_data['num_%s' % (i+1)] = data[(list2[1][i]):(list2[1][i+1])]
  cycles_binwidth.append(int(cycles_data['num_%s' % (i+1)].shape[0] / 16))
  cycles_width.append(cycles_data['num_%s' % (i+1)].shape[0] / 16)
  cycles_length.append(cycles_data['num_%s' % (i+1)].shape[0])
#for i in range(24): # visualising SSN time series after dividing it into cycles
#  plt.plot(cycles_data['num_%s' % (i+1)]['yearfrac'], cycles_data['num_%s' % (i+1)]['ssn'])
#plt.xlim(1750, 2030), plt.xlabel('Year')
#plt.ylim(0, 450), plt.ylabel('SSN')

  
### preprocessing the data: splitting cycles into 16 phases
weights = [] # list that stores the weight of all months within each solar cycle
phase_numbers = [] # list of phase numbers for all months of all cycles
for i in range(24):
    temp = [] # weights for i-th cycle
    phnum = [] # phase numbers for i-th cycle
    for m in range(cycles_length[i]):
      temp.append(0)
      phnum.append(0)
    #phase 1
    temp[:cycles_binwidth[i]] = [1]*cycles_binwidth[i]
    phnum[:cycles_binwidth[i]] = [1]*cycles_binwidth[i]
    temp[cycles_binwidth[i]] = cycles_width[i] - cycles_binwidth[i]
    phnum[cycles_binwidth[i]] = 1.5
    lst = temp[cycles_binwidth[i]] # save the weight of the last month in the phase
    #phases 2-15
    for q in range(14):
      nxt = 1 - lst # calculate the weight of the same month in the next phase
      nmnth = int(cycles_width[i]-nxt) # calculate the amount of months with weight of 1 in each phase
      for p in range(nmnth):
        temp[temp.index(0)] = 1
        phnum[phnum.index(0)] = q+2
      temp[temp.index(0)] = cycles_width[i]-nxt-nmnth
      phnum[phnum.index(0)] = q+2.5
      lst = cycles_width[i]-nxt-nmnth
    #phase 16
    while (temp.index(0) < cycles_length[i]-1):
      temp[temp.index(0)] = 1
    while (phnum.index(0) < cycles_length[i]-1):
      phnum[phnum.index(0)] = 16
    temp[-1] = 1
    phnum[-1] = 16
    weights.append(temp)
    phase_numbers.append(phnum)
name_ssn = "phase_ssn_"
ssn_matrix = list()
for k in range(16): #averaging SSN inside each phase
  exec(name_ssn + "%s = list()" % (k+1))
  for i in range(24):
    wd = cycles_width[i]
    ssn_wsum = 0
    for j in range(cycles_length[i]): #for all months in the i-th cycle we weigh ssn and sum it up
      if (abs(phase_numbers[i][j]-1 - k) <= 0.5): #if this month is inside the phase or on its border
        ssn_wsum += cycles_data['num_%s' % (i+1)]['ssn'].values[j]*weights[i][j]
      else:
        ssn_wsum += 0
    exec(name_ssn + "%s.append(round(ssn_wsum/wd, 2))" % (k+1))
  exec("ssn_matrix.append(" + name_ssn + "%s)" % (k+1))
#for j in range(16): # to see the multidimensional time series in the matrix form. in the row from left to right - the ssn of cycles from the 1st to the 24th
#  for i in range(24): #in columns from top to bottom - ssn of phases from 1st to 16th
#    print(ssn_matrix[j][i], end='\t') 
#  print('\n')
plt.plot(ssn_matrix) # visualising all cycles' shapes and their progression with phases
plt.xlabel('Phase number (starting with 0)')
plt.ylabel('Weighted SSN')
plt.show()
ssn_transposed = np.transpose(ssn_matrix) 
n = 24 # to compare real montly mean SSN with weighted average SSN for solar cycle 24
fig, ax1 = plt.subplots()
ax2 = ax1.twiny()
ax1.plot(range(1,17), ssn_transposed[n-1], color='red')
ax2.plot(cycles_data['num_%s' % n]['yearfrac'], cycles_data['num_%s' % n]['ssn'], color='black')
ax1.set_xlim(1,16)
ax2.set_xlim(cycles_data['num_%s' % n]['yearfrac'].values[0],cycles_data['num_%s' % n]['yearfrac'].values[-1])
ax1.set_xlabel('Phase number'), ax2.set_xlabel('Year')
ax1.set_ylabel('SSN of solar cycle %s' % n)
plt.show()


### spectral analysis: decomposition into the generalized Fourier series with Walsh basis functions
ssn_tmp = []
for i in range(24): #we apply the Walsh transform to the matrix 24*16 line by line, the transformed one is located according to Hadamard numbering
  ssn_tmp.append(fwht(ssn_transposed[i]))
numer = [0,8,12,4,6,14,10,2,3,11,15,7,5,13,9,1] #place the values inside the tmp in this order: Walsh's basic functions' numbers, arranged according to the Gray code
ssn_walsh = []
tmpp = []
for i in range(24):
  tmpp = ssn_tmp[i]
  ssn_walshh = []
  for j in numer:
    ssn_walshh.append(tmpp[j])
  ssn_walsh.append(ssn_walshh)
#for i in range(16):
#  for j in range(24):
#    print(format(np.transpose(ssn_walsh)[i][j], '.5g'), end='\t')
#  print('\n')


### trend decomposition
# qoddt - decomposition coefficients for odd cycles after trend extraction using tablecurve program. for even cycles in the 3rd line there would be: range(1,25,2)
qoddt1, qoddt2, qoddt3, qoddt4, qoddt5, qoddt6, qoddt7, qoddt8, qoddt9, qoddt10, qoddt11, qoddt12, qoddt13, qoddt14, qoddt15, qoddt16 = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
ssn_waltp = np.transpose(ssn_walsh)
for i in range(0, 24, 2):
  qoddt1.append(np.array(ssn_walsh[i][0]-1111.234-23.7181*(i+1)))
  qoddt2.append(np.array(ssn_walsh[i][1]+60.7837-169.0191*np.log(i+1)))
  qoddt3.append(np.array(ssn_walsh[i][2]-575.45+12.31*(i+1)))
  qoddt4.append(np.array(ssn_walsh[i][3]+84.41+29.83*np.log(i+1)))
  qoddt5.append(np.array(ssn_walsh[i][4]))
  qoddt6.append(np.array(ssn_walsh[i][5]+52.96+8.13*(i+1)))
  qoddt7.append(np.array(ssn_walsh[i][6]+175.33+10.38*(i+1)))
  qoddt8.append(np.array(ssn_walsh[i][7]-65.65+5.32*(i+1)))
  qoddt9.append(np.array(ssn_walsh[i][8]+20.03-2.59*(i+1)))
  qoddt10.append(np.array(ssn_walsh[i][9]-51.25+20.25*np.log(i+1)))
  qoddt11.append(np.array(ssn_walsh[i][10]))
  qoddt12.append(np.array(ssn_walsh[i][11]))
  qoddt13.append(np.array(ssn_walsh[i][12]))
  qoddt14.append(np.array(ssn_walsh[i][13]-3.16+19.84*np.log(i+1)))
  qoddt15.append(np.array(ssn_walsh[i][14]+108.56+4.10*(i+1)))
  qoddt16.append(np.array(ssn_walsh[i][15]+56.91-12.22*np.log(i+1)))
qoddt = np.array([qoddt1, qoddt2, qoddt3, qoddt4, qoddt5, qoddt6, qoddt7, qoddt8, qoddt9, qoddt10, qoddt11, qoddt12, qoddt13, qoddt14, qoddt15, qoddt16])
for j in range(16): # let's look at the R^2 offset of the series before and after the trend selection
  fig1 = stats.probplot(qoddt[j], plot=plt, rvalue=True)
  fig2 = stats.probplot(ssn_waltpn[j], plot=plt, rvalue=True)
  if round(fig1[-1][-1]**2-fig2[-1][-1]**2, 3) <= 0: # in case trend subtraction somehow makes residuals less random, we develop a special way to mitigate the effect
      qoddt[j] = ssn_waltpn[j]
      print('Q' + str(j+1) + ' before trend subtraction (best result): r^2 = ' + format(fig2[-1][-1]**2, '.4g'))
  else:
      print('Q' + str(j+1) + ' R^2 improvement after trend subtraction: +' + format(fig1[-1][-1]**2-fig2[-1][-1]**2, '.3g') + ' resulting in r^2 = ' + format(fig1[-1][-1]**2, '.4g'))


### seasonal decomposition
qodds1, qodds2, qodds3, qodds4, qodds5, qodds6, qodds7, qodds8, qodds9, qodds10, qodds11, qodds12, qodds13, qodds14, qodds15, qodds16 = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
for i in range(12):
  qodds1.append(qoddt[0][i]-338.2801*math.sin(0.683*(2*i+1)+0.8495)-1380.2966)
  qodds2.append(qoddt[1][i]-302.3497*math.sin(0.714*(2*i+1)-0.4368)+27.965)
  qodds3.append(qoddt[2][i]-198.1385*math.sin(0.757*(2*i+1)+3.0164)+1147.2838)
  qodds4.append(qoddt[3][i]-149.9748*math.sin(0.8727*(2*i+1)+1.9902)+151.7084)
  qodds5.append(qoddt[4][i]-92.5597*math.sin(0.757*(2*i+1)+2.1645)+84.7953)
  qodds6.append(qoddt[5][i]-214.1079*math.sin(0.748*(2*i+1)+2.3752)-10.6195)
  qodds7.append(qoddt[6][i]-75.3903*math.sin(0.6614*(2*i+1)-1.0671)+305.5487)
  qodds8.append(qoddt[7][i]-64.8825*math.sin(0.748*(2*i+1)-0.7352)+1.3626)
  qodds9.append(qoddt[8][i]-25.5982*math.sin(0.5984*(2*i+1)+0.7234)-8.1498)
  qodds10.append(qoddt[9][i]-48.1767*math.sin(1.0134*(2*i+1)+1.6129)-8.3122)
  qodds11.append(qoddt[10][i]-89.1441*math.sin(1.0472*(2*i+1)-1.308)-77.635)
  qodds12.append(qoddt[11][i]-61.8904*math.sin(0.8727*(2*i+1)-0.973)-46.4486)
  qodds13.append(qoddt[12][i]-27.6512*math.sin(1.0649*(2*i+1)+3.0701)+36.5165)
  qodds14.append(qoddt[13][i]-81.5248*math.sin(0.6981*(2*i+1)+2.3238)+31.8725)
  qodds15.append(qoddt[14][i]-56.0413*math.sin(0.3927*(2*i+1)-2.7903)-11.4992)
  qodds16.append(qoddt[15][i]-23.9952*math.sin(0.4189*(2*i+1)+2.4171)-4.2817)
qodds = np.array([qodds1, qodds2, qodds3, qodds4, qodds5, qodds6, qodds7, qodds8, qodds9, qodds10, qodds11, qodds12, qodds13, qodds14, qodds15, qodds16])
for j in range(16): # let's look at the R^2 offset of the series before and after the seasonal selection
  fig1 = stats.probplot(qoddt[j], plot=plt, rvalue=False)
  fig2 = stats.probplot(qodds[j], plot=plt, rvalue=False)
  if round(fig2[-1][-1]**2-fig1[-1][-1]**2, 3) <= 0:
      qodds[j] = qoddt[j]
      print('Q' + str(j+1) + ' before season subtraction (best result): r^2 = ' + format(fig1[-1][-1]**2, '.4g'))
  else:
      print('Q' + str(j+1) + ' R^2 improvement after season subtraction: +' + format(fig2[-1][-1]**2-fig1[-1][-1]**2, '.3g') + ' resulting in r^2 = ' + format(fig2[-1][-1]**2, '.4g'))


### extrapolation
# unfortunately, it's done manually with excel to make sure that the trend & season functions are valid. four cases were considered, resulting in ssn_extra list:
# 1. R^2 improved after both trend and season subtraction. we extrapolate both trend and season and sum the results (the trend-season model is additive)
# 2. R^2 didn't improve after both trend and season subtraction. worst-case scenario is to not predict this coefficient at all
# in case 2, we extrapolate both trend and season function and weigh the extrapolations by their respective r^2 (the r^2 for approximation)
# 3. R^2 improved after trend subtraction but not after season subtraction. for some of these cases r^2 of season approximation was too significant to neglect season
# therefore in case 3 we acted the same as in case 2
# 4. R^2 didn't improve after trend subtraction but improved after season subtraction. we extrapolate the season and neglect the trend extrapolation
ssn_extra = [1295.349561, 155.6813, -804.412, -289.8425581, -12.19604234, -73.5463, -414.3491047, -57.90837766, 26.09183344, 39.41312584, 57.96008169, 100.2450119, -63.84666414, 24.63560442, -161.5289, 3.248601]
ssn_walsh.append(ssn_extra)

### convolution
numer1 = [0, 15, 7, 8, 3, 12, 4, 11, 1, 14, 6, 9, 2, 13, 5, 10] # moving on from Gray's numbering to Hadamard's numbering
ssn_forecasted = []
tmppp = []
for i in range(25):
  tmppp = ssn_walsh[i]
  ssn_forecastedd = []
  for j in numer1:
    ssn_forecastedd.append(tmppp[j])
  ssn_forecasted.append(ssn_forecastedd)
ssn_res = []
for i in range(25):  #we apply the inverse Walsh_Hadamard transform to the matrix 25*16 line by line, the transformed one is located according to Hadamard
  ssn_res.append(ifwht(ssn_forecasted[i]))
for j in range(16):
  if ssn_res[24][j] < 0:
    ssn_res[24][j] = 0 # we check for inconsistencies with the physics: the number of sunspots cannot be negative

### forecasting results
# It remains to associate the height of the cycle with its length. This can be done according to statistical patterns (such as the Waldmeier effect), but we prefer another
# because the Waldmeier effect didn't work for cycle 24, and the anti-correlation between the ascending branch duration and the amplitude of cycle corresponds to r^2 = 0.77
# just 0.77 isn't enough to make decision for Waldmeier effect, and there aren't any more significant correlations for montly mean SSN in solar science
# therefore, we fit our forecast by two points of real data. point 1 is the beginning of the cycle in December 2019, obviously
# point 2 is July 2023, when there were 124 spots on average per month, which corresponds to averaged SSN of 124 in the fifth phase of our forecast
# 44 months have passed since the beginning of the cycle until July 2023 => 44*16/5 = 140.8 months or 11 years and 8.8 months
# then the peak of 184 spots is predicted for September 2025, and the end is for August 2031
data_25 = pd.read_csv("sc25input.csv", sep=';', header = 0)
fig, ax1 = plt.subplots()
ax2 = ax1.twiny()
ax1.plot(range(1,17), ssn_res[24], color='red')
ax2.plot(data_25['yearfrac'], data_25['ssn'], color='black')
#ax1.axhline(y=124, xmax = 0.26, color='black', linestyle='--')
ax1.set_xlim(1,16)
ax1.set_ylim(0,200)
ax2.set_xlim(2019.958, 2019.958+16*8.8/12)
ax2.set_xlabel('Year')
ax1.set_ylabel('SSN for solar cycle 25')
ax1.tick_params(axis='both', direction='in', labelbottom=0, bottom=0), ax2.tick_params(axis='both', direction='in')
ax2.grid(visible='True', which='major', linestyle='--')
plt.show()
