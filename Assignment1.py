import matplotlib.pyplot as plt
import numpy as np

# Start time & end time

n_tot = 2000
start_time = 0
end_time = 2.5

h = (end_time - start_time)/n_tot
t = np.linspace(start_time, end_time, n_tot)

# Function definitions

def cent_dif(y, n):
    y[n+1] = (((m/(h*h))*(2*y[n] - y[n-1])) + ((c/(2*h))*y[n-1]) + (-k*y[n]) + force(n))/((m/(h*h)) + (c/(2*h)))
    return y

def anal_soln(m, c, k, m0, e, omega, v0, y0, t):
    
    # System variables
    
    zeta = c/(2*np.sqrt(m*k))
    omega_n = np.sqrt(k/m)
    omega_d = omega_n*np.sqrt(1-zeta**2)
    
    # Derived variables
    
    f_0 = m0*e*(omega**2)
    r = omega/omega_n
    psi = np.arctan(2*zeta*r/(1-(r**2)))
    
    psi = psi + np.pi if r>1 else psi 
    
    X_mag = f_0/np.sqrt((k-m*omega**2)**2 + c**2*omega**2)
    
    # Forced response
    
    y_c = X_mag*np.cos(omega*t - psi)
    
    # More variables woop
    
    A_3 = -X_mag*np.cos(-psi)
    A_4 = (omega*X_mag*np.sin(-psi) + zeta*omega_n*A_3)/omega_d
    y_p = (A_3*np.cos(omega_d*t) + A_4*np.sin(omega_d*t))*np.exp(-zeta*omega_n*t);
    
    r = omega/omega_n
    
    return (y_p + y_c), r

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

# Function definitions

testing = False
set_by_r = True
numerical_estimate = True
auto_time_step = False

[m, c, k] = [54, 320, 26000]
[m_0, e, omega] = [5.5, 0.225, 400*np.pi/30]

[y_0, v_0] = [0, 0]

m = m + m_0

r = 0.5 if set_by_r else omega/np.sqrt(k/m)
omega = r*np.sqrt(k/m) if set_by_r else omega
    
def force(n):
    return m_0*e*(omega**2)*np.cos(omega*(n*h))
    
y_s, r = anal_soln(m=m, k=k, c=c, m0=m_0, e=e, omega=omega, v0=v_0, y0=y_0, t=t)

# Position array definitions, previous position approximation

t = np.append(t, -h)

y = np.zeros(n_tot+1)
y[0] = y_0
y[-1] = (y_0 - (v_0 - (h**2/2)*((force(0) - c*v_0 - k*y_0)/m)))

# Main loop

n = 0

while (n < n_tot):
    y = cent_dif(y, n)    
    n += 1 

# Plot output
fig, ax = plt.subplots(dpi=300)
plt.plot(t[:-1], y_s*1000, color='dodgerblue', label="Analytical solution")
if numerical_estimate:
    plt.scatter(t[:-1], y[:-1]*1000, s=1, color='darkred', marker='o', label="Numerical approximation")
plt.title("Washing machine displacement, omega = "+ str(round(omega,2)) +" rad/s, r = " + str(round(r, 3)))
plt.legend(loc="upper right")
ax.margins(y=0.3)
plt.xlabel("Time (s)")
plt.ylabel("Displacement (mm)")
plt.show()

# Error checking

y = y[:-1]

# Time step checking
if testing:
    error_check(y, y_s)
    omega_n = np.sqrt(k/m)
    print("m c k =", m, c, k)
    print("Number of points = ", n_tot,", Time step = ", round(h,6), ", natural frequency = ", round(omega_n,3), sep='')
    print("Time step ", round((h/omega_n)*100,3), "% of natural frequency", sep='')
    