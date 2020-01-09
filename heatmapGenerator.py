import matplotlib.pyplot as plt
import numpy as np


def drawOscilationGraph(file):
    data = np.loadtxt(file, delimiter=',')
    print(data)
    plt.imshow(data, cmap="gray", interpolation='none', aspect='auto')
    plt.yticks([])
    plt.show()

def drawAmplitudeGraph(file):
    data = np.loadtxt(file, delimiter=',')
    print(data)
    plt.imshow(data, cmap="viridis", interpolation='none', aspect='auto')
    plt.yticks([])
    plt.show()

def drawPeriodGraph(file):
    data = np.loadtxt(file, delimiter=',')
    print(data)
    plt.imshow(data, cmap="viridis", interpolation='none', aspect='auto')
    plt.yticks([])
    plt.show()

#drawOscilationGraph('extended_model_1_v3.txt')
#drawOscilationGraph('extended_model_2_v3.txt')
drawAmplitudeGraph("extended_model_1_amplitudes.txt")
drawAmplitudeGraph("extended_model_1_periods.txt")