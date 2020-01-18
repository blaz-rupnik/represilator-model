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
    plt.imshow(data, cmap="YlGnBu", interpolation='none', aspect='auto')
    plt.colorbar()
    plt.yticks([])
    if show:
        plt.show()
    timestr = time.strftime("%Y%m%d-%H%M%S")
    plt.savefig('img/' + name + "_" + timestr + ".png")

def drawPeriodGraph(file, name, show):
    data = np.loadtxt(file, delimiter=',')
    print(data)
    plt.imshow(data, cmap="YlGnBu", interpolation='none', aspect='auto')
    plt.colorbar()
    plt.yticks([])
    if show:
        plt.show()
    timestr = time.strftime("%Y%m%d-%H%M%S")
    plt.savefig('img/'+name+"_"+timestr+".png")

drawOscilationGraph('extended_model_2_test.txt', "model_2_oscilations", False)
#drawOscilationGraph('extended_model_2_v3.txt', "model_2_oscilations", False)
drawAmplitudeGraph("extended_model_2_amplitudes_test.txt", "model_2_amplitude", False)
drawPeriodGraph("extended_model_2_periods_test.txt", "model_2_periods", False)
