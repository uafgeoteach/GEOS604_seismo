# Homework 2: Waves on a String
- **Course:** UAF GEOS604 &mdash; Seismology  
- **Instructor:** Bryant Chow ([bhchow@alaska.edu](bhchow@alaska.edu))
- **Course Website:** [https://bryantchow.com/teaching/geos604](https://bryantchow.com/teaching/geos604) 
- **Last Modified:**  01/23/25

## Semester: Spring 2025
- **Total Points**: 10
- **Assigned**: Jan. 23, 2025
- **Due Date and Time**: Jan. 30, 2025 at the beginning of class


## Submitting Homework

1. If you have handwritten solutions to any problems:  
   a. include them directly in the notebook as a scan/picture in the appropriate section, OR  
   b. scan/photograph and include them separately. please also reference them in the notebook ("see handwritten notes submitted separately")
2. Save your homework as a PDF file using `File` -> `Save and Export Notebook As` -> `PDF`
3. Double check that the output PDF contains all your code, text, and any images you have included
5. Upload your PDF to the class Google Drive (see course website above) before the **Due Date and Time** listed above.

**See the class syllabus section on "Homework Policy" for late policy and information on academic integrity.**

---------------

## Motivation

This homework assignment is meant to give you some intuition of waves, their key properties, and two classic solutions to the wave equation, the harmonic wave solution (traveling waves), and the standing wave/normal mode solution.  Some things I want this homework to build your intuition and familiarity on are:
1. The 1D wave equation
2. Common solutions to the wave equation (travelling and standing waves) and their superpositions
3. Reflection and transmission behaviors, dictated by the 1D wave equation

## Background

Remember from the lecture that we derived the 1D wave equation looking at the forces acting on a string, and asking how a displaced string would respond to these forces in space and time:

$$v^2 \frac{\partial^2 u(x, t)}{\partial x^2} = \frac{\partial^2 u(x, t)}{\partial t^2}$$

- The velocity $v$ that we derived for the string is dependent on the  material ($\rho$) and the tension applied to the string: $v=\sqrt{\tau/\rho}$.  
- In lecture we saw that two solutions to the wave equation came in the form of travelling waves, and standing waves.


```python
# Feel free to import any other packages you need here.
import numpy as np
import matplotlib.pyplot as plt
from geos604funcs import wave_eq_preset, normal_mode_seismogram, normal_mode_snapshot
%matplotlib widget  
```

---------------
## Problem 1: Travelling Waves [2.0 points]

Remember in lecture we described on solution to the wave equation as the harmonic wave solution.  Let's build some understanding about the harmonic wave solution, the parameters we use to describe it, and how it defines waves which travel in both space and time.  

Consider below the real part of the harmonic wave solution that we derived in class, which describes a wave traveling in the $+\hat{x}$ direction (see Euler's formula to go from the general solution to this version):

$$ u(x, t) = Acos(\omega t - kx)  $$

You are free to use the following function to define the harmonic wave solution for your solutions. Remember that there are some plotting examples in HW1 that will help you solve the problems below.
```python
def u(x, t, A, w, k):
    """
    Harmonic Wave Solution

    .. usage::
        x = np.linspace(0, 2 * np.pi, 1000)
        y = u(x, t=0, A=1, w=1, k=1)

    :type x: list
    :param x: distance spanned by the soution
    :type t: float
    :param t: time at which you want to look at u(x,t)
    :type A: float
    :param A: maximum amplitude of the wave
    :type w: float
    :param w: omega, angular frequency
    :type k: float
    :param k: wavenumber 
    """
    return A * np.cos(w * t - k * x)
```

-----------

1. [0.50] Plot $u(x, t)$ at time $t=0$ and for distance $0 \le X \le 2\pi$. Choose a value of $k$ so that there are four wavelengths contained in $X$ (set the other variables to $A=1$ and $\omega=2\pi$).

   a. What value of $k$ was required?  
   b. On the same figure, pick a point $x_0$ where $u$ is at a maximum value (where $u(x_0, t) = A$). Plot a marker on your plot at $u(x_{0}, 0)$

3. [0.75] Plot $u(x,t_i)$ for times $t_i=0.0, 0.25, 0.5, 0.75$ (please make one plot per time stamp, $t_i$). Also plot the a marker on each figure at $u(x_0, t_i)$.

   a. From these plots can you tell which way this wave is travelling?  
   b. Imagine you are an observer at the location of this marker. Write a few words describing what you are experiencing. Can you think of a real-world analogy for this situation?  
   c. What happens if you double the wavenumber $k$? What would an observer at $u(x_0, t_i)$ notice if all of a sudden the wavenumber doubled? (you don't need to plot this, but you can if it helps).  
   d. What would happen if you double the angular frequency $\omega$? From the point of view of the observer on the fixed marker, what will they notice if the angular frequency suddenly doubled? (you don't need to plot this, but you can if it helps).

   
3. [0.50] In the last part we considered the situation of a stationary observer/fixed point in space. Now we want to travel with the wave.  

   a. Plot a series of 4 figures for $u(x,t_i)$ for times $t_i=0.0, 0.25, 0.5, 0.75$
   - For each of the figures, place a marker at some location $x_i$ where $u(x_i, t_i) = u(x_{i+1}, t_{i+1})$. In other words, the marker should travel with the same part of the wave for each of the four figures.
   - **Note:** Don't do this manually, but find a some general relationship for $u(x_i, t_i)$.
    
   b. Same as in the previous situation, imagine an observer at the location of this marker. Write a few words describing what the observer is experiencing. Can you think of a real-world analogy for this situation?

   
4. [0.25] Please explain why the solution $cos(\omega t - kx)$ gives you a wave travelling in the direction you stated in Part 3, whereas $cos(\omega t + kx)$ gives you waves travelling in the other direction? You may find it helpful to explain this in terms of the math.  


```python
# Your answers here
```

---------------
## Problem 2: Waves in Space and Time  [2.0 points]

The solution $u(x, t) = Acos(\omega t - kx)$ defines a wave that travels in both space in time. In the last problem we looked at space continuously, but only snapshots in time. But this wave exists in both space and time, as shown in Stein and Wysession's Fig 2.2-3.  

Matplotlib has a nice example for plotting a surface here that can be modified to for our purposes (last accessed 1/17/25): https://matplotlib.org/stable/gallery/mplot3d/surface3d.html  

The following code block takes the important parts of that example and provides the setup for the problem, you will modify this function to make your own plot.

```python
# Import statements are already taken care of
fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(8,8))

# Set the data, fill in the '?' parts yourself
x = # ?
t = # ?
X, T = np.meshgrid(x, t)

# Wave constants
w = # ?
A = # ?
k = # ?

# Harmonic wave solution
U = # u(x,t) = ?

# Plot the surface, use a nice colormap
surf = ax.plot_surface(X, T, U, cmap="viridis")

# You can add anything else here to help you visualize

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=10, label="u(x,t)")

plt.show()
```
-------

1. [0.25] Using the surface plotting function above, plot the harmonic wave solution $u(x, t) = A cos(\pi t - 2 \pi x)$ in both space and time for $0 \le X \le 4$ and $0 \le t \le 4$ with $A=100$
2. [0.50] Think about a line that follows on continuous part of this wave in space and time (e.g., a continuous peak or continuous valley), what does this represent?
3. [0.25] If you defined the solution as $u(x, t) = A cos(\pi t + 2 \pi x)$, how would this change the quantity you identified in Part 2? (plotting this may help)
4. [1.00] Make a new plot where you recreate the wave defined in Problem 1.3. I.e., modify the ranges of $X$ and $T$, and the wave constants defined in Problem 1.3. Can you describe each of the figures you made in Problem 1.3 with this 3D plot. Some examples of what I'm looking for
   - "If I draw a line through this part of the 3D plot, I get the waveform shown at $t=t_0$"
   - "If I was to plot the marker shown in Problem 1.3, it would follow this path..."


```python
# Your answers here
```

---------------
## Problem 3: Reflection and Transmission [2.0 points]

Now we want to try and mimic the situation shown in Fig. 2.2-6 in Stein and Wysession to build some intuition about how solutions to the wave equation behave when they cross impedance contrasts (boundaries between different materials). This means we need to understand how the wave equation behaves at boundaires, and dictates the tradeoff of space and time for valid solutions. Let's look at one valid solution to the wave equation, the Gaussian.  

Consider a Guassian pulse moving in the $+\hat{x}$ direction, that satisfies the 1D wave equation:

$$ u(x, t) = Ae^{-((x-vt)-b)^2/ 2c^2} $$

where $A$ is the amplitude, $b$ is the starting position. I found a *sweet* 1D wave equation solver online that I have modified for use in this notebook. We will call the function `wave_eq_preset()` to play a few variations of this animation in our notebook.

#### Some notes on this animation:
> - If you re-run the cell containing the animation, it will restart. If you let it run for a while, it will also restart.
> - `Kernel` -> `Interrupt Kernel` will stop the animation
> - If the animation does not play, try 1) re-running the cell, 2) restarting the kernel, 3) restart JupyterLab. If it still does not work, please email/Slack message me!
> - Under the hood, the function uses a finite-difference scheme to solve the partial differential equation numerically.
> - All credit for the original function goes to the author Sacha Binder ([GitHub repository](https://github.com/sachabinder/wave_equation_simulations/tree/main)).

-----------
1. First, let's just make sure this thing works. Run the function `wave_eq_preset()` in a cell.
2. [1.0] Explain what you see in the animation by answering the following questions:  
   a. Is this string that is hosting the wave homogeneous across the entire domain? If not where is the junction?  
   b. Which side ($v_1$ or $v_2$) has the higher velocity?    
   c. Starting at $t\approx1.3$s, what happens at the left-hand boundary ($x=0$)? Explain in terms of acoustic impedance.  
   d. Starting at $t\approx2.1$s, what happens at the right-hand boundary ($x=1.5$)? Explain in terms of acoustic impedance.  
3. [0.5] Can you determine the velocities $v_1$ and $v_2$ of the string on either side of the junction? What are they?
4. [0.5] Now run `wave_eq_preset(choice=1)`. Is $v_2 > v_1$ or is $v_1 > v_2$? Describe the features of the animation that justify your answer.


```python
# Your answers here
```

---------------
## Problem 4: Standing Waves [4.0 points]

Another solution to the wave equation that we saw in class was the standing wave, a wave that oscillates in place. We want to build some intution on standing waves, and try and recreate Stein and Wysession Fig 2.2-8 to understand how summations of normal modes/standing waves allow us to generate traveling waves.    

Consider the standing wave solution we looked at in class (SW Eq. 2.2.40). 

$$ u(x, t) = \sin(n\pi x/L)\cos(\omega t) $$

Where $L$ is the length of the string. A traveling wave can be expressed as the weighted sum of a string's normal modes, which can be written (SW Eq. 2.2.42).

$$ u(x, t) = \sum_{n=0}^{\infty} sin(n\pi x_s/L)F(\omega_n)\sin(n\pi x/L)\cos(\omega_n t) $$

Where:
- $n$ gives integer values for each mode of the string
- $x_s$ is the position of a source which excites these modes (or the travelling wave)
- $\omega_n$ is the angular frequency ($\omega_n = \pi v n / L$, which comes from $\omega = kv$ and $\lambda = 2L/n$)
- $F(\omega_n)$ is a weight factor that describes how different frequencies contribute to the time history of the source

We can set the weighting term as a source term:

$$ F(\omega_n) = exp({-(\omega_n \tau)^2 / 4)}) $$

To make things easier I coded up the normal mode summation problem following the example in Stein and Wysession Section A8.1. I will ask you to play around with the function to build understanding of normal modes and how they can be summed to generate travelling waves, and similar to the 

> #### Some notes on these functions:
> - Both functions describe a string of finite length. An initial source term is "plucked" at $t=0$.
> - You can choose how many modes you want to sum using the parameter `number_modes`. Default value is 40
> - By default the functions plot all the modes (as in SW Fig 2.2-8). If your value of `number_modes` is too large this may be difficult to visualize. You can set `show_modes=False` to only plot the sum (snapshot or seismogram).

------------

1. Let's first run these functions
   - In one cell, run function `normal_mode_seismogram`
   - In the following cell, run `normal_mode_snapshot`    
   - For both functions set the parameters `number_modes=40` and `show_modes=True`
2. [1.0] Looking at the output figures:  
   a. Describe what you see in 'Seismogram'  
   b. What do you think the three pulses are in 'Seismogram'?  
   c. Describe what you see in 'Snapshot'  
   d. What is the yellow star in 'Snapshot'  
   e. What timestep $t$ is shown in 'Snapshot'  
3. [0.5] Decrease the number of modes by half and re-run the functions. What happens? What happens if you double the number of modes?   
    a. In 'Snapshot', we might expect to plot a smooth travelling wave like we saw in Problem 3, but that is not the case. Can you write a few words describing why we do not see a smooth wave?  
    b. What we would need to do to get a smooth wave? Try to make a smooth wave in 'Seismogram', what did you have to do?  
    
4. [1.0] In `normal_mode_snapshot`, you can choose the timestep to plot by setting parameter `time=1` (e.g., to get time 1, you can change the value to get later times). Play around with this and convince yourself you can see the travelling wave propagating.  
   a. Referencing these snapshots, where is the receiver location in 'Seismogram'?  
   b. In terms of acoustic impedance, describe the boundary conditions at the two sides of the string.  
   c. What is the velocity of the string?
6. [0.5] If you look at the code (you don't have to), I did not explicitely code up any boundary conditions, I simply plotted the equation for $u(x,t)$ (above). If that is the case, how does the travelling wave know what to do at the boundaries? (hint: this has to do with the lecture derivation)
7. [1.0] The source term is described by $F(\omega_n)$ described above.  
   a. Plot $F(\omega_n)$ for the first 200 modes with $\tau=0.4$. You should have all the variables available to you from information in this problem.  
   b. What function does $\tau$ describe? (I'm looking for a name)  
   c. Play around with the value of $\tau$, what does $\tau$ control?
   d. You can set $\tau$ in the functions `normal_mode_seismogram` and `normal_mode_snapshot` with the input parameter `tau=0.2` (or whatever value you want). How does changing $\tau$ affect your 'Seismogram'? Does this match your expectation from your plot of $\tau$?  
   e. What do you think would happen if you change the function that controls $\tau$? What would change in the plots from `normal_mode_seismogram` and `normal_mode_snapshot`? (you do not need to plot this, just describe)


```python
# Your answers here
```

-----------
## Problem N 
1. Approximately how much time did you spend on this homework assignment?
2. Did you find this homework particularly easy, adequate, or difficult?
3. Any feedback on the homework?



```python
# Your answers here
```
