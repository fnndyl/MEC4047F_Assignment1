import matplotlib.pyplot as plt
import numpy as np

# Start time & end time

n_tot = 565
start_time = 0
end_time = 10

h = (end_time - start_time)/n_tot
t = np.linspace(start_time, end_time, n_tot)

# Options: noDamp, damp, sinF

mode = 'noDamp'
testing = True

if mode == 'userDef':
    
    [m, c, k] = [50, 10, 200]
    [y_0, v_0] = [5, 3]

    [m_0, e] = [6, 0.225]
    
    [y_0, v_0] = [1, 0]
    
    def force(n):
        return 30*np.sin((n*h))
    
elif mode == 'noDamp':
    def force(t):
        return 0
    
    [m, c, k] = [50, 0, 250]
    [y_0, v_0] = [1, 2]
    
    # Real soln
    y_s = np.cos(np.sqrt(5)*t) + (2/np.sqrt(5))*np.sin(np.sqrt(5)*t)
    
elif mode == 'damp':
    force = lambda t : 0
    
    [m, c, k] = [50, 10, 150]
    [y_0, v_0] = [1, 2]
    
    # Real soln
    x_s = np.arange(start_time, end_time, 0.01)
    y_s = (np.cos(np.sqrt(299)/10*t) + (21/np.sqrt(299))*np.sin(np.sqrt(299)/10*t))*np.exp((-1/10)*t)
    
elif mode == 'sinF':
    def force(n):
        return 10*np.sin((n*h))
    
    [m, c, k] = [1, 5, 6]
    [y_0, v_0] = [0, 5]
    
    # Real soln
    x_s = np.arange(start_time, end_time, 0.01)
    y_s = (-6*np.exp(-3*t) + 7*np.exp(-2*t) + np.sin(t) - np.cos(t))

# Function definitions

def cent_dif(y, n):
    y[n+1] = (((m/(h*h))*(2*y[n] - y[n-1])) + ((c/(2*h))*y[n-1]) + (-k*y[n]) + force(n))/((m/(h*h)) + (c/(2*h)))
    return y

def error_check(y, y_s):
    tolerance = 0.01
    
    y_peak, y_rebound, y_s_peak, y_s_rebound = max(y), min(y), max(y_s), min(y_s)

    peak_error = abs((y_s_peak - y_peak)/y_s_peak)
    rebound_error = abs((y_s_rebound - y_rebound)/y_s_rebound)
    
    if ((peak_error > tolerance) or (rebound_error > tolerance)):
        print("\nOut of tolerance - peak or rebound error greater than tolerance error")
        print("Tolerance = ", tolerance*100, "%", sep='')
        valid = False
    else:
        print("\nFirst peak and first rebound error within tolerance")
        print("Tolerance = ", tolerance*100, "%", sep='')
        valid = True
    print("Peak error = ", round(peak_error*100, 3), "%, rebound error = ", round(rebound_error*100, 3), "%\n", sep='')
    return valid

# Position array definitions

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
    plt.plot(t[:-1], y_s, color='darkmagenta')

plt.scatter(t[:-1], y[:-1], s=1, color='violet', marker='o')
plt.show()

# Error checking

y = y[:-1]

error_check(y, y_s)

# Time step checking
omega_n = np.sqrt(k/m)
print("Time step = ", round(h,6), ", natural frequency = ", round(omega_n,3), sep='')
print("Time step ", round((h/omega_n)*100,3), "% of natural frequency", sep='')
    