import matplotlib.pyplot as plt
import numpy as np
import time


def drawOscilationGraph(file, name, show):
    data = np.loadtxt(file, delimiter=',')
    print(data)
    plt.imshow(data, cmap="gray", interpolation='none', aspect='auto')
    plt.yticks([])
    if show:
        plt.show()
    timestr = time.strftime("%Y%m%d-%H%M%S")
    plt.savefig('img/' + name + "_" + timestr + ".png")

def drawAmplitudeGraph(file, name, show):
    data = np.loadtxt(file, delimiter=',')
    print(data)
    plt.imshow(data, cmap="viridis", interpolation='none', aspect='auto')
    plt.yticks([])
    if show:
        plt.show()
    timestr = time.strftime("%Y%m%d-%H%M%S")
    plt.savefig('img/' + name + "_" + timestr + ".png")

def drawPeriodGraph(file, name, show):
    data = np.loadtxt(file, delimiter=',')
    print(data)
    plt.imshow(data, cmap="viridis", interpolation='none', aspect='auto')
    plt.yticks([])
    if show:
        plt.show()
    timestr = time.strftime("%Y%m%d-%H%M%S")
    plt.savefig('img/'+name+"_"+timestr+".png")

drawOscilationGraph('extended_model_1_v3.txt', "model_1_oscilations", False)
drawOscilationGraph('extended_model_2_v3.txt', "model_2_oscilations", False)
drawAmplitudeGraph("extended_model_1_amplitudes.txt", "model_1_amplitude", False)
drawPeriodGraph("extended_model_1_periods.txt", "model_1_periods", False)
