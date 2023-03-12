import matplotlib.pyplot as plt
import numpy as np

# Physical constants

m = 50
c = 10
k = 250

m_0 = 6
e = 0.225

# Initial conditions

y_0 = 1
v_0 = 0


y_neg1 = (0 - c*(v_0) - k*(y_0)) # This is still terrible


y = [0.9997, y_0]
print(y)


# Function definitions

def cent_dif(y, n):
    
    next_pos = (((m/(h*h))*(2*y[n] - y[n-1])) + ((c/(2*h))*y[n-1]) + (-k*y[n]))/((m/(h*h)) + (c/(2*h)))
    y.append(next_pos)

    return y

# Time definitions

n_tot = 1000

start_time = 0
end_time = 10

h = (end_time - start_time)/n_tot

t_plot = [0, 0]
t_plot[0] = 0 - h


# Main loop

t = start_time
n = 1

while (t < end_time):
    
    t = n*h
    
    y = cent_dif(y, n)
    t_plot.append(t)
    
    n += 1 

# Print plot

fig, ax = plt.subplots()
ax.plot(t_plot, y)
