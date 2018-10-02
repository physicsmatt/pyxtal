import numpy as np
import matplotlib.pyplot as plt
from scipy.misc import imread
import matplotlib.cbook as cbook


#np.random.seed(0)
#x = np.random.uniform(0.0,10.0,15)
#y = np.random.uniform(0.0,10.0,15)
#
#img = imread("cat.jpg")
##plt.scatter(x,y,zorder=1)
#plt.imshow(img, zorder=0, extent=[0.5, 8.0, 1.0, 7.0])
#plt.show()

fig,ax = plt.subplots()
x = range(300)
ax.imshow(imread("cat.jpg"), extent=[0, 300, 0, 300],zorder=0)
ax.plot(x, x, '--', linewidth=5, color='firebrick',zorder=1)
plt.show()
