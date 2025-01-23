""" 
General behind the scenes functions for homework notebooks when the backend is not important for course material
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from numpy import e, pi, sin, cos
plt.rcParams["animation.html"] = "jshtml"  # needed to show animations inline 


def wave_eq_animation(length_of_string=1.5, duration_of_simulation=3,
                      junction=0.7, c_1=1.0, c_2=0.5, source_location=0,
                      source_location_2=None, source_1_sign=1, source_2_sign=1,
                      left_bound_cond="neumann", right_bound_cond="dirichlet",
                      show_junction=False):
    """
    Solve the 1D wave equation with finite differences and then plot an animation

    .. note::
        This is not my function, it was copied and modifed from:
        https://github.com/sachabinder/wave_equation_simulations/blob/main/1D_WAVE-EQ_variable-velocity.py
    
        Along with helper functions anim_1D() and gaussian()
    
        All credit goes to the original author, Sacha Binder
    
        Citation: Sacha BINDER. Étude de l’observation et de la modélisation des ondes 
        de surface en eau peu profonde. Physique-Informatique.TIPE session 2021.
        
        I have modified the codes to work as a callable function with modifiable input 
        variables so that we may use it in a Jupyter notebook for UAF GEOS604. I also changed 
        some of the variable names to be easier to read. I have not touched any of the 
        corefunctionality of the underlying FD scheme.

        
    Solve the 1D Wave equation and plot
    This file was built to solve numerically a classical PDE, 1D wave equation. 

    The numerical scheme is based on finite difference method. This program is also 
    providing several boundary conditions. More particularly the Neumann, Dirichlet 
    and Mur boundary conditions.
    
    Copyright - © SACHA BINDER - 2021
    """
    assert(left_bound_cond in ["neumann", "dirichlet", "mur"])
    assert(right_bound_cond in ["neumann", "dirichlet", "mur"])

    #Spatial mesh - i indices
    L_x = length_of_string #Range of the domain according to x [m]
    dx = 0.01 #Infinitesimal distance
    N_x = int(L_x/dx) #Points number of the spatial mesh
    X = np.linspace(0,L_x,N_x+1) #Spatial array

    #Temporal mesh with CFL < 1 - j indices
    L_t = duration_of_simulation #Duration of simulation [s]
    dt = 0.01 * dx  # Infinitesimal time with CFL condition
    N_t = int(L_t/dt) #Points number of the temporal mesh
    T = np.linspace(0,L_t,N_t+1) #Temporal array

    #Velocity array for calculation (finite elements)
    c = np.zeros(N_x+1, float)
    for i in range(0,N_x+1):
        if X[i] <= junction:
            c[i] = c_1
        else:
            c[i] = c_2

    # Calculate Constants
    C2 = (dt/dx)**2
    CFL_1 = c_1*(dt/dx)
    CFL_2 = c_2*(dt/dx)

    
    # Begin Processing loop, setup empty FD arrays
    u_jm1 = np.zeros(N_x+1,float)   #Vector array u_i^{j-1}
    u_j = np.zeros(N_x+1,float)     #Vector array u_i^j
    u_jp1 = np.zeros(N_x+1,float)   #Vector array u_i^{j+1}

    # TO DO: What does q define?
    q = np.zeros(N_x+1,float)
    q[0:N_x+1] = c[0:N_x+1]**2

    # Global solution of displacement u(x,t)
    U = np.zeros((N_x+1,N_t+1),float) 
    
    # Initial condition u(x,0), allow for 2 sources
    u_j[0:N_x+1] = source_1_sign * gaussian(X[0:N_x+1], b=source_location)
    if source_location_2:
        u_j[0:N_x+1] += source_2_sign * gaussian(X[0:N_x+1], b=source_location_2)
    U[:,0] = u_j.copy()
    
    # Initial condition u(x,t=1) with FD
    u_jp1[1:N_x] = u_j[1:N_x] + 0.5 * C2 * (0.5*(q[1:N_x] + q[2:N_x+1])*(u_j[2:N_x+1] - u_j[1:N_x]) - 
                                            0.5*(q[0:N_x-1] + q[1:N_x])*(u_j[1:N_x] - u_j[0:N_x-1]))
    
    # Impose Left Hand Side (LHS) boundary condition
    if left_bound_cond == "dirichlet":
        u_jp1[0] = 0
    elif left_bound_cond == "neumann":
        u_jp1[0] = u_j[0] + 0.5 * C2 * (0.5*(q[0] + q[0+1])*(u_j[0+1] - u_j[0]) 
                                      - 0.5 * (q[0] + q[0+1])*(u_j[0] - u_j[0+1]))
    elif left_bound_cond == "mur":
        u_jp1[0] = u_j[1] + (CFL_1-1)/(CFL_1 + 1)*( u_jp1[1] - u_j[0])
    
    # Impose Right Hand Side (RHS) boundary condition
    if right_bound_cond == "dirichlet":
        u_jp1[N_x] = 0
    elif right_bound_cond == "neumann":
        u_jp1[N_x] =  u_j[N_x] + 0.5*C2*(0.5*(q[N_x-1] + q[N_x])*(u_j[N_x-1] - u_j[N_x]) - 
                                         0.5*(q[N_x-1] + q[N_x])*(u_j[N_x] - u_j[i-1])) 
    elif right_bound_cond == "mur":
        u_jp1[N_x] = u_j[N_x-1] + (CFL_2 -1)/(CFL_2 + 1)*(u_jp1[N_x-1] - u_j[N_x])

    # Setting up for next time step
    u_jm1 = u_j.copy() 
    u_j = u_jp1.copy() 
    U[:,1] = u_j.copy()
    
    # Begin time loop starting at t=2
    for j in range(1, N_t):
        # Calculation at step j+1
        u_jp1[1:N_x] = -u_jm1[1:N_x] + 2*u_j[1:N_x] + C2*(0.5*(q[1:N_x] + q[2:N_x+1])*(u_j[2:N_x+1] - u_j[1:N_x]) - 
                                                          0.5*(q[0:N_x-1] + q[1:N_x])*(u_j[1:N_x] - u_j[0:N_x-1]))
        
        # Impose Left Hand Side (LHS) boundary condition
        if left_bound_cond == "dirichlet":
            u_jp1[0] = 0
        elif left_bound_cond == "neumann":
            u_jp1[0] = -u_jm1[0] + 2*u_j[0] + C2*(0.5*(q[0] + q[0+1])*(u_j[0+1] - u_j[0]) - 
                                                  0.5*(q[0] + q[0+1])*(u_j[0] - u_j[0+1]))       
        elif left_bound_cond == "mur":
            u_jp1[0] = u_j[1] + (CFL_1-1)/(CFL_1+1)*( u_jp1[1] - u_j[0])

        # Impose Right Hand Side (RHS) boundary condition
        if right_bound_cond == "dirichlet":
            u_jp1[N_x] = 0               
        elif right_bound_cond == "neumann":
            u_jp1[N_x] = -u_jm1[N_x] + 2*u_j[N_x] + C2*(0.5*(q[N_x-1] + q[N_x])*(u_j[N_x-1] - u_j[N_x]) - 
                                                        0.5*(q[N_x-1] + q[N_x])*(u_j[N_x] - u_j[N_x-1]))            
        elif right_bound_cond == "mur":
            u_jp1[N_x] = u_j[N_x-1] + (CFL_2 -1)/(CFL_2 + 1)*(u_jp1[N_x-1] - u_j[N_x])

        # Increment for next time step
        u_jm1[:] = u_j.copy() 
        u_j[:] = u_jp1.copy()  
        U[:,j] = u_j.copy()

    # Plot the animation 
    if not show_junction:
        junction = None
    anim1 = anim_1D(X, U, dt, 10, xlim=(0, length_of_string), 
                    ylim=(-1.0,1.5), junction=junction)

    plt.show()
    

def normal_mode_seismogram(
        string_len_m=20., velocity=3., number_modes=40, 
        source_position_m=8, rcv_position_m=16., 
        seismogram_duration_s=10, nstep=500, tau=.02, show_modes=True, 
        **kwargs):
        # Parameters from SW A8.1
        # string_len_m=1., velocity=1., number_modes=40, 
        # source_position_m=0.2, rcv_position_m=0.7, 
        # seismogram_duration_s=1.25, 
    """
    Sum normal modes to create a traveling wave and resulting seismogram.
    Transcribed from Stein and Wysession A.8.1 following SW section 2.2

    :param string_len_m: Total length of string L in unit meters
    :param velocity: Speed of the wave in meters per second (m/s)
    :param number_modes: Total number of modes to sum up from n=0
    :param source_position: Where is the input source in unit meters
    :param rc_position_m: Where is the receiver in unit meters
    :param seismogram_duration_s: Duration of the recorded waveform in unit sec.
    :param nstep: Number of time steps
    :param dt: (Optional) length of one time step in seconds
    :param tau: Term controlling the shape of the source pulse
    :type show_modes: bool
    :param show_modes: Plot the mode wiggles or just the sum
    """
    # Initialize Empty Displacement field for t=0
    time = np.linspace(0, seismogram_duration_s, nstep)
    displacement = np.zeros(nstep)
    modes = []
    
    # Each mode needs to be defined as a term in the summation of n
    for mode in range(1, number_modes + 1):
        # Spatial terms
        source_term = sin(mode * pi * source_position_m / string_len_m)
        rcv_term = sin(mode * pi * rcv_position_m / string_len_m)

        # Time Independent Terms
        omega_n = pi * velocity * mode / string_len_m  # Mode frequency 
        f_omega_n = e ** (-1 * (omega_n * tau)**2 / 4)  # Weighting Term

        # Generate the time-dependent contribution of this mode, m
        mode_disp = source_term * f_omega_n * rcv_term * cos(omega_n * time)
        modes.append(mode_disp)
        displacement += mode_disp

    if show_modes:
        plot_modes(x=time, y=displacement, modes=modes, xlabel="Time",
                   title="Seismogram", **kwargs)
    else:
        f, ax = plt.subplots(figsize=(8, 3))
        ax.plot(time, displacement, c="k")
        plt.xlabel("Time")
        plt.ylabel("Displacement")
        plt.show()


def normal_mode_snapshot(
        time=0, string_len_m=20., velocity=3., number_modes=40, 
        source_position_m=8., rcv_position_m=None, tau=.02, show_modes=True, 
        **kwargs):
    """
    Sum normal modes to create a traveling wave and resulting seismogram.
    Transcribed from Stein and Wysession A.8.1 following SW section 2.2

    :param string_len_m: Total length of string L in unit meters
    :param velocity: Speed of the wave in meters per second (m/s)
    :param number_modes: Total number of modes to sum up from n=0
    :param source_position: Where is the input source in unit meters
    :param rc_position_m: Where is the receiver in unit meters
    :param seismogram_duration_s: Duration of the recorded waveform in unit sec.
    :param nstep: Number of time steps
    :param dt: (Optional) length of one time step in seconds
    :param tau: Term controlling the shape of the source pulse
    :type show_modes: bool
    :param show_modes: Plot the mode wiggles or just the sum
    """
    nx = 1000
    # Initialize Empty Displacement field for t=0
    dx = string_len_m / nx
    distance = np.linspace(0, string_len_m, nx)
    displacement = distance * 0
    modes = []
    
    # Each mode needs to be defined as a term in the summation of n
    for mode in range(1, number_modes + 1):
        # Spatial terms
        source_term = sin(mode * pi * source_position_m / string_len_m)
        rcv_term = sin(mode * pi * distance / string_len_m)

        # Time Independent Terms
        omega_n = pi * velocity * mode / string_len_m  # Mode frequency 
        f_omega_n = e ** (-1 * (omega_n * tau)**2 / 4)  # Weighting Term

        # Generate the time-dependent contribution of this mode, m
        mode_disp = source_term * f_omega_n * rcv_term * cos(omega_n * time)
        modes.append(mode_disp)
        displacement += mode_disp

    if show_modes:
        f, axs = plot_modes(x=distance, y=displacement, modes=modes, 
                            title="Snapshot", xlabel="Distance", **kwargs)
        axs[1].plot([source_position_m], [2.5], marker="*", c="y", markersize=10)
        if rcv_position_m:
            axs[1].plot([rcv_position_m], [2.5], marker="v", c="r")
    else:
        f, ax = plt.subplots(figsize=(8, 3))
        ax.plot(distance, displacement, c="k")
        plt.xlabel("Distance")
        plt.ylabel("Displacement")
        
    plt.show()


def plot_modes(x, y, modes, xlabel, label_every=5, title="",
               figsize=(6, 8), dpi=100, show_grid=True, 
               **kwargs):
    """
    General plotting function to plot mode summation for synthetic
    seismogram or for snapshot displacements

    :type x: np.array
    :param x: time or distance array depending on function
    :type y: np.array
    :param y: displacement
    :type modes: list of arrays
    :param modes: each mode array
    :type label_every: int
    :param label_every: how often to put a text label for modes
    :type xlabel: str
    :param xlabel: what to write on the bottom axis
    """
    # Plot the seismogram created by mode summation
    f, axs = plt.subplots(figsize=figsize, dpi=dpi, nrows=2, 
                          height_ratios=[12,1], sharex=True)
    
    for m, mode in enumerate(modes):
        # Normalize modes between -0.5 and 0.5 so they fit in a space of one
        # then offset modes by one so they stack on top of each other, 
        # starting at the top with m=0 then plotting downwards
        sep = 1.25  # separation factor between adjacent modes
        scaled_mode = mode / (2 * mode.max()) + (len(modes) - sep * m)
        axs[0].plot(x, scaled_mode, c="k", lw=0.5)
        
        # Text labels, index starts at 0 but modes start at 1 so we need an offset
        label = m + 1
        # Write mode number every `label_every` modes
        if (not label % label_every) or (label in [1, len(modes)]): 
            axs[0].text(0, len(modes) - sep * m, f"{label:>2}", c="r")  

    # Plot the summed modes
    axs[1].plot(x, y, c="k")

    # We don't care about the y axis mode values, just the representations
    axs[0].get_yaxis().set_ticks([])
    if show_grid:
        axs[0].grid(ls="--", alpha=0.5)
        axs[1].grid(ls="--", alpha=0.5)

    axs[0].set_ylabel("Modes")
    axs[1].set_xlabel(xlabel)
    axs[1].set_ylabel("Sum")

    if title:
        axs[0].set_title(title)

    # Remove empty whitespace between figures
    plt.subplots_adjust(wspace=0, hspace=0)

    return f, axs
    
# =============== Wave Equation Animation Helper Functions ====================
def gaussian(x, b):
    """
    Single space variable fonction that represent the wave form at t = 0 with 
    starting location `b`
    """
    return np.exp(-(x-b)**2/0.01)
    

def wave_eq_preset(choice=0, percent=None):
    """
    Example 1 of wave eq. animation, so that student's aren't exposed to the 
    chosen variables
    """
    if choice == 0:
        wave_eq_animation()
    elif choice == 1:
        wave_eq_animation(c_2=1.5)
    elif choice == 2:
        c_1 = 1.0
        c_2 = c_1 * (1 + percent)
        wave_eq_animation(c_2=c_2)

        
def anim_1D(x, y, pas_de_temps, pas_d_images, junction=None, 
            xlim=(0, 4), ylim=(-4, 4)):
    """
    Function allowing to display an animation of the wave equation based on 
    calculation results with a given time step. This function can be used to 
    save the images sequence in the current directory.

    .. TO DO::
        Clean up parameter naming
    
    The y parameter is a list containing several functions to display : y
     = [ [f_1(x)], ... , [f_n(x)] ].
    
    (x:np.ndarray (format 1D), y:np.ndarray (format 2D), pas_de_temps:float , 
    pas_d_images:int, save:bool , myxlim:tuple , myylim:tuple) -> plot (+ .mp4)
    """
    fig = plt.figure()
    ax = plt.axes(xlim=xlim , ylim=ylim)
    line, = ax.plot([], [])
    ax.set_title("t = 0 s")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("$u$ [m]")
    ax.text(0.2, 1.2, "$v_1$", fontsize=15)
    ax.text(1.1, 1.2, "$v_2$", fontsize=15)
    ax.grid()
    if junction:
        plt.axvline(junction, ls="--", c="k")

    def init():
        line.set_data([],[])
        return line,
    
    # animation function.  This is called sequentially
    def animate(i):
        line.set_data(x, y[:,pas_d_images*i])
        ax.set_title("t = {:.2f} s".format(np.round(i*pas_d_images*pas_de_temps, 4)))
        return line,
        
    # call the animator. blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, init_func=init, 
                                   frames=y.shape[1]//pas_d_images, interval=10)

    return anim
