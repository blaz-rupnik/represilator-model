import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt('extended_model_1_v2.txt', delimiter=',')
print(data)
plt.imshow(data, cmap='Pastel1', interpolation='none', aspect='auto')
plt.yticks([])
plt.show()