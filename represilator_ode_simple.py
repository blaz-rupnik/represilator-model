import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy import signal
import time
import libs





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

def simulate_extended_model_1():
    with open('results/extended_model_1_v3.txt', 'w') as f, \
            open('results/extended_model_1_amplitudes.txt', 'w') as fa, \
            open('results/extended_model_1_periods.txt', 'w') as fp:
        result = []
        for i in range(1, 5):
            new = []
            new_amp = []
            new_per = []
            logspace_vector = np.logspace(0,5,1000).astype(int)

            for j in range(0, len(logspace_vector)):
                params = (alpha, alpha0, beta, n, i, i, i, logspace_vector[j], logspace_vector[j], logspace_vector[j])
                Z = odeint(libs.extended_model_1, Z0, t, args=params)
                A = Z[:, 3]
                is_osc, period, amplitude = libs.osc_detect(A)
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
    with open('results/extended_model_2_v3.txt', 'w') as f, \
            open('results/extended_model_2_amplitudes.txt', 'w') as fa, \
            open('results/extended_model_2_periods.txt', 'w') as fp:
        result = []
        for i in range(1, 5):
            new = []
            new_amp = []
            new_per = []
            logspace_vector = np.logspace(0, 5, 1000).astype(int)

            for j in range(0, len(logspace_vector)):
                params = (alpha, alpha0, beta, n, i, i, i, logspace_vector[j], logspace_vector[j], logspace_vector[j])
                Z = odeint(libs.extended_model_2, Z0, t, args=params)
                A = Z[:, 3]
                is_osc, period, amplitude = libs.osc_detect(A)
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

start = time.time()
first_test = libs.simulate_extended_model_1()
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

