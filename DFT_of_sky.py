from math import *
from scipy import array, exp
import random
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import *
import matplotlib.cm as cm




c = 
kb =             # Boltzman's constant in SI
jansky = 
bandwidth =       #30.72 MHz
accTime =               #accTime is a correlator accumulation time
frequency =    
wavelength =  
mega = 

###################################### Generating Noise for the 496 visibilities with rmsFluxDensity  #####################################
# note : noise(noiserms, bandwidth, accTime) is the function for estimating a thermal noise
randomnoise = np.random.normal(0,noise(noiserms, bandwidth, accTime),496) 


############################## Generating random position for each antenna for a BL 1m to (less)3km distance ##############################
n = 32
coordsx = []
coordsy = []
p = 0
for p in range (0,n):
    coordsx.append(random.randint(1, 3000)) 
    coordsy.append(random.randint(1, 3000)) 

i=0
j=0;
while i < n:
	j = i + 1
	i = i+1
	while j < n-1:
		j=j+1
		#print i,j
		if sqrt((coordsx[j] - coordsx[i])**2 + (coordsy[j] - coordsy[i])**2)>=3000:
			coordsx[i] = random.randint(1, 3000)
			coordsy[i] = random.randint(1, 3000)
			i=0	
			j=n

i=i+ 1
j=j+ 1


########################################## reading the 40 sources properties    ###################################

x = open('any_astronomical_catalouge.dat')
y = x.readlines()
x.close()

flux = []
ra = []
dec = []
freq = []

for i in y:
        flux.append(float((i.split()[0])))     
	ra.append(float(i.split()[1]))             
	dec.append(float(i.split()[2]))
        freq.append(float((i.split()[3]))*mega)

ra0 = 60 * pi/180.0
dec0 = -30 * pi/180.0

ll = []
mm = []
for d in range(0,40):
	ll.append(cos(dec[d])*sin(ra[d]-ra0))
	mm.append(sin(dec[d])*cos(dec0) - cos(dec[d])*sin(dec0)*cos(ra[d]-ra0))

print flux
print ra
######################################################## defining u and v by the phase center   ########################################################

u = []
v = []
for i in  range (0,n): 
	for j in  xrange (i+1 ,n):   
		u.append((sin(ra0)*(coordsx[j] - coordsx[i])+cos(ra0)*(coordsy[j] - coordsy[i]))/wavelength) 
		v.append((-sin(dec0)*cos(ra0)*(coordsx[j] - coordsx[i])+sin(dec0)*sin(ra0)*(coordsy[j] - coordsy[i]))/wavelength)


vis = []

for k in  range (0,496):
        visibility = randomnoise[k]
	for d in range(0,40):
        	visibility += (flux[d] * exp(-2 * pi * 1j * (u[k]*ll[d] + v[k]*mm[d])))
   	vis.append(visibility)    	

print vis
print len(vis)

####################################################  DFT (Discrete Fourier Transform) of the visibilities   ##################################################

matrix = []
i = 0
for m in numpy.arange(-0.08,0.2,0.0002): #  resolution :0.0007 
        pix = []
	j = 0
	for l in numpy.arange(-0.08,0.2,0.0002):
        	intensity = 0
		for k in range(0,len(vis)):
                	intensity += vis[k] * exp(2 * pi * 1j * (u[k]*l + v[k]*m))
        	pix.append((intensity/len(vis)).real)
	matrix.append(pix)        
	
	j=j+1
i=i+1

#print matrix	
print len(pix)
#print len(matrix)
########################################################### Plotting ###########################################################

plt.imshow(matrix,origin= 'lower',extent = [-0.08,0.2,-0.08,0.2], aspect='auto')
#plt.grid(True,color='white')
plt.title('This is DFT of (Vis)')
plt.xlabel('J2000 Right Ascension')
plt.ylabel('J2000 Declination')
plt.colorbar()
plt.show()

