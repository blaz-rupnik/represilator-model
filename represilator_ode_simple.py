import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy import signal
import time

# function that returns dZ/dt
def model(Z,t,alpha, alpha0, beta, n):
    a = Z[0]
    b = Z[1]
    c = Z[2]
    A = Z[3]
    B = Z[4]
    C = Z[5]
    
    # mRNAs
    dadt = alpha/(1 + C**n) + alpha0 - a
    dbdt = alpha/(1 + A**n) + alpha0 - b
    dcdt = alpha/(1 + B**n) + alpha0 - c
    
    # proteins
    dAdt = - beta * (A - a)
    dBdt = - beta * (B - b)
    dCdt = - beta * (C - c)
    
    
    dzdt = [dadt, dbdt, dcdt, dAdt, dBdt, dCdt]
    
    return dzdt

def extended_model_1(Z,t,alpha, alpha0, beta, n, nA, nB, nC, kA, kB, kC):
    a = Z[0]
    b = Z[1]
    c = Z[2]
    A = Z[3]
    B = Z[4]
    C = Z[5]
    
    # mRNAs
    dadt = alpha/(1 + C**n + (A/kA)**nA) + alpha0 - a
    dbdt = alpha/(1 + A**n + (B/kB)**nB) + alpha0 - b
    dcdt = alpha/(1 + B**n + (C/kC)**nC) + alpha0 - c
    
    # proteins
    dAdt = - beta * (A - a)
    dBdt = - beta * (B - b)
    dCdt = - beta * (C - c)
    
    
    dzdt = [dadt, dbdt, dcdt, dAdt, dBdt, dCdt]
    
    return dzdt 

def extended_model_2(Z,t,alpha, alpha0, beta, n, nA, nB, nC, kA, kB, kC):
    a = Z[0]
    b = Z[1]
    c = Z[2]
    A = Z[3]
    B = Z[4]
    C = Z[5]
    
    # mRNAs
    dadt = (alpha * ((A/kA)**nA))/(1 + C**n + (A/kA)**nA) + alpha0 - a
    dbdt = (alpha * ((B/kB)**nB))/(1 + A**n + (B/kB)**nB) + alpha0 - b
    dcdt = (alpha * ((C/kC)**nC))/(1 + B**n + (C/kC)**nC) + alpha0 - c
    
    # proteins
    dAdt = - beta * (A - a)
    dBdt = - beta * (B - b)
    dCdt = - beta * (C - c)
    
    
    dzdt = [dadt, dbdt, dcdt, dAdt, dBdt, dCdt]
    
    return dzdt

# initial condition
#Z0 = [10,0,0,0,0,0]
Z0 = [ 7.02028482, 1.57288111, 58.50876295, 8.69333201, 1.40155783, 53.49915358]

# number of time points
t_end = 100
n = t_end * 10

# time points
t = np.linspace(0, t_end, n)


alpha = 216
alpha0 = 0.001 * alpha
n = 2
beta = 5


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
        #If this signal is valid, we need to find the period of oscilation
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

def simulate_extended_model_1():
    with open('extended_model_1_v3.txt', 'w') as f, \
            open('extended_model_1_amplitudes.txt', 'w') as fa, \
            open('extended_model_1_periods.txt', 'w') as fp:
        result = []
        for i in range(1, 5):
            new = []
            new_amp = []
            new_per = []
            logspace_vector = np.logspace(0,5,1000).astype(int)

            for j in range(0, len(logspace_vector)):
                params = (alpha, alpha0, beta, n, i, i, i, logspace_vector[j], logspace_vector[j], logspace_vector[j])
                Z = odeint(extended_model_1, Z0, t, args=params)
                A = Z[:, 3]
                is_osc, period, amplitude = osc_detect(A)
                new.append(is_osc)
                new_amp.append(amplitude)
                new_per.append(period)
            f.write("%s\n" % ','.join(str(x) for x in new))
            fa.write("%s\n" % ','.join(str(x) for x in new_amp))
            fp.write("%s\n" % ','.join(str(x) for x in new_per))
            result.append(new)
            print("-- starting new row --")
    f.close()
    print("file was closed")
    return result

def simulate_extended_model_2():
    with open('extended_model_2_v3.txt', 'w') as f, \
            open('extended_model_2_amplitudes.txt', 'w') as fa, \
            open('extended_model_2_periods.txt', 'w') as fp:
        result = []
        for i in range(1, 5):
            new = []
            new_amp = []
            new_per = []
            logspace_vector = np.logspace(0, 5, 1000).astype(int)

            for j in range(0, len(logspace_vector)):
                params = (alpha, alpha0, beta, n, i, i, i, logspace_vector[j], logspace_vector[j], logspace_vector[j])
                Z = odeint(extended_model_2, Z0, t, args=params)
                A = Z[:, 3]
                is_osc, period, amplitude = osc_detect(A)
                new.append(is_osc)
                new_amp.append(amplitude)
                new_per.append(period)
            f.write("%s\n" % ','.join(str(x) for x in new))
            fa.write("%s\n" % ','.join(str(x) for x in new_amp))
            fp.write("%s\n" % ','.join(str(x) for x in new_per))
            result.append(new)
            print("-- starting new row --")
    f.close()
    print("file was closed")
    return result

#params= (alpha, alpha0, beta, n)

#Z = odeint(model, Z0, t, args=params)

#params_1 = (alpha, alpha0, beta, n, 2, 2, 2, 100, 100, 100)
#Z_1 = odeint(extended_model_1, Z0, t, args=params_1)

#Z_2 = odeint(extended_model_2, Z0, t, args=params_1)

#test = simulate_extended_model_1()
#print(test)

#A = Z[:,3]
#B = Z[:,4]
#C = Z[:,5]

start = time.time()
first_test = simulate_extended_model_1()
end = time.time()
print(end - start)
#print(first_test)

# plot results
#fig = plt.figure()

#plotting results of first extended model
#ax = fig.add_subplot(111)
#ax.imshow(first_test, aspect='auto', interpolation='nearest', cmap=plt.cm.gray)
#ax.plot(t,A,'g:',label='A(t)')
#ax.plot(t,B,'b-',label='B(t)')
#ax.plot(t,C,'r--',label='C(t)')
#ax.set_ylabel('concentrations')
#ax.set_xlabel('time')
#ax.legend(loc='best')
#plt.show()

##################
##################
#################

# Generate amplitude heatmaps

