if mode == 'userDef':
    
   [m, c, k] = [54, 320, 26000]
   [m_0, e, omega] = [5.5, 0.225, 400*np.pi/30]
   
   [y_0, v_0] = [0, 0]
   
   m = m + m_0
   
   def force(n):
       return m_0*e*(omega**2)*np.cos(omega*(n*h))
    
elif mode == 'noDamp':
    def force(n):
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
    y_s = (np.cos(np.sqrt(299)/10*t) + (21/np.sqrt(299))*np.sin(np.sqrt(299)/10*t))*np.exp((-1/10)*t)
    
elif mode == 'sinF':
    def force(n):
        return 10*np.sin((n*h))
    
    [m, c, k] = [1, 5, 6]
    [y_0, v_0] = [0, 5]
    
    # Real soln
    y_s = (-6*np.exp(-3*t) + 7*np.exp(-2*t) + np.sin(t) - np.cos(t))