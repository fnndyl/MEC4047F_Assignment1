import matplotlib.pyplot as plt
import numpy as np

# Start time & end time

n_tot = 1000
start_time = 0
end_time = 10

# Options: noDamp, damp, sinF, omegaF

mode = 'sinF'

if mode == 'userDef':
    def force(n):
        return 30*np.sin((n*h))
    
    [m, c, k] = [50, 10, 200]
    [y_0, v_0] = [5, 3]

    m_0 = 6
    e = 0.225
    
    [y_0, v_0] = [1, 0]
    
elif mode == 'noDamp':
    def force(t):
        return 0
    
    [m, c, k] = [50, 0, 250]
    [y_0, v_0] = [1, 2]
    
    # Real soln
    x_s = np.arange(start_time, end_time, 0.01)
    y_s = np.cos(np.sqrt(5)*x_s) + (2/np.sqrt(5))*np.sin(np.sqrt(5)*x_s)
    
elif mode == 'damp':
    force = lambda t : 0
    
    [m, c, k] = [50, 10, 150]
    [y_0, v_0] = [1, 2]
    
    # Real soln
    x_s = np.arange(start_time, end_time, 0.01)
    y_s = (np.cos(np.sqrt(299)/10*x_s) + (21/np.sqrt(299))*np.sin(np.sqrt(299)/10*x_s))*np.exp((-1/10)*x_s)
    
elif mode == 'sinF':
    def force(n):
        return 10*np.sin((n*h))
    
    [m, c, k] = [1, 5, 6]
    [y_0, v_0] = [0, 5]
    
    # Real soln
    x_s = np.arange(start_time, end_time, 0.01)
    y_s = (-6*np.exp(-3*x_s) + 7*np.exp(-2*x_s) + np.sin(x_s) - np.cos(x_s))

# Function definitions

def cent_dif(y, n):
    y[n+1] = (((m/(h*h))*(2*y[n] - y[n-1])) + ((c/(2*h))*y[n-1]) + (-k*y[n]) + force(n))/((m/(h*h)) + (c/(2*h)))
    return y

# Critical array definitions

h = (end_time - start_time)/n_tot
t = np.linspace(start_time, end_time, n_tot)
t = np.append(t, -h)

y = np.zeros(n_tot+1)
y[0] = y_0
y[-1] = y_0 - (v_0 - h*((force(0) - c*v_0 - k*y_0)/m))*h

# Main loop

n = 0

forces = []

plt.figure(dpi=300)

while (n < n_tot):
    y = cent_dif(y, n)    
    n += 1 

# Plot output

if mode != 'userDef':
    plt.plot(x_s, y_s, color='darkmagenta')
    
#if mode == 'sinF':
#    plt.plot(t[:-1], force(t[:])[:-1])

plt.scatter(t[:-1], y[:-1], s=1, color='violet', marker='o')
plt.show()
