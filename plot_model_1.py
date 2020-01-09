import numpy as np
import matplotlib.pyplot as plt

#read from file and convert to matrix
data = np.loadtxt('extended_model_1_v3.txt', delimiter=',')
print('loaded')

# plot results
fig = plt.figure()

#plotting results of first extended model
ax = fig.add_subplot(111)
ax.imshow(data, aspect='auto', interpolation='none', cmap='Pastel1')
plt.show()