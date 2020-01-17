import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy import signal
import time

# function that returns dZ/dt
def model(Z, t, alpha, alpha0, beta, n):
    a = Z[0]
    b = Z[1]
    c = Z[2]
    A = Z[3]
    B = Z[4]
    C = Z[5]

    # mRNAs
    dadt = alpha / (1 + C ** n) + alpha0 - a
    dbdt = alpha / (1 + A ** n) + alpha0 - b
    dcdt = alpha / (1 + B ** n) + alpha0 - c

    # proteins
    dAdt = - beta * (A - a)
    dBdt = - beta * (B - b)
    dCdt = - beta * (C - c)

    dzdt = [dadt, dbdt, dcdt, dAdt, dBdt, dCdt]

    return dzdt


def extended_model_1(Z, t, alpha, alpha0, beta, n, nA, nB, nC, kA, kB, kC):
    a = Z[0]
    b = Z[1]
    c = Z[2]
    A = Z[3]
    B = Z[4]
    C = Z[5]

    # mRNAs
    dadt = alpha / (1 + C ** n + (A / kA) ** nA) + alpha0 - a
    dbdt = alpha / (1 + A ** n + (B / kB) ** nB) + alpha0 - b
    dcdt = alpha / (1 + B ** n + (C / kC) ** nC) + alpha0 - c

    # proteins
    dAdt = - beta * (A - a)
    dBdt = - beta * (B - b)
    dCdt = - beta * (C - c)

    dzdt = [dadt, dbdt, dcdt, dAdt, dBdt, dCdt]

    return dzdt


def extended_model_2(Z, t, alpha, alpha0, beta, n, nA, nB, nC, kA, kB, kC):
    a = Z[0]
    b = Z[1]
    c = Z[2]
    A = Z[3]
    B = Z[4]
    C = Z[5]

    # mRNAs
    dadt = (alpha * ((A / kA) ** nA)) / (1 + C ** n + (A / kA) ** nA) + alpha0 - a
    dbdt = (alpha * ((B / kB) ** nB)) / (1 + A ** n + (B / kB) ** nB) + alpha0 - b
    dcdt = (alpha * ((C / kC) ** nC)) / (1 + B ** n + (C / kC) ** nC) + alpha0 - c

    # proteins
    dAdt = - beta * (A - a)
    dBdt = - beta * (B - b)
    dCdt = - beta * (C - c)

    dzdt = [dadt, dbdt, dcdt, dAdt, dBdt, dCdt]

    return dzdt


def max_diff(a):
    diff = 0
    temp = a[:]
    temp.sort()
    diff = abs(temp[0] - temp[-1])
    return diff


def osc_detect(a):
    peaks = signal.find_peaks(a)
    num_of_peaks = len(peaks[0])

    if num_of_peaks <= 3:
        return 0, 0, 0

    peak_values = a[peaks[0]]

    if (abs(peak_values[-3] - peak_values[-2])) <= 0.001:
        # If this signal is valid, we need to find the period of oscilation
        peaks_i = peaks[0]
        print("Peaks:", peaks_i)
        f = [x - peaks_i[i - 1] for i, x in enumerate(peaks_i)][2:-1]
        print("Diffs:", f)
        print("Values:", peak_values)
        diff = max_diff(f)
        print("Max diff:", diff)
        if diff > 10:
            return 0, 0, 0
        return 1, f[0], peak_values[-2]

    return 0, 0, 0