# Homework 4: Ray Theory
- **Course:** UAF GEOS604 &mdash; Seismology  
- **Instructor:** Bryant Chow ([bhchow@alaska.edu](bhchow@alaska.edu))
- **Course Website:** [https://bryantchow.com/teaching/geos604](https://bryantchow.com/teaching/geos604) 
- **Last Modified:**  02/10/25

## Semester: Spring 2025
- **Total Points**: 10
- **Assigned**: February 11, 2025
- **Due Date and Time**: February 18, 2025 at the beginning of class

## Submitting Homework

1. If you have handwritten solutions to any problems:  
   a. include them directly in the notebook as a scan/picture in the appropriate section, OR  
   b. scan/photograph and include them separately. please also reference them in the notebook ("see handwritten notes submitted separately")
2. Please write any text-based answers in Markdown cells.
3. Please keep code cells to ~80 character line-width as exporting to PDF often cuts off code blocks.
4. Save your homework as a PDF file using `File` -> `Print` -> `Save to PDF`
5. Double check that the output PDF contains all your code, text, and any images you have included
6. Upload your PDF to the class Google Drive (see course website above) before the **Due Date and Time** listed above.

**See the class syllabus section on "Homework Policy" for late policy and information on academic integrity.**

---------------

## Motivation



```python
# Import cell, please run this cell before proceeding
# You may import any other packages you may need
import numpy as np
import matplotlib.pyplot as plt
from obspy.taup import TauPyModel, plot_travel_times, plot_ray_paths
from obspy.geodetics import locations2degrees
%matplotlib widget
```

# Problem 1: Ray Theory [2.5]

> **NOTE: No part of Problem 1 requires code. Problem 1 can be solved with pen and paper, or LaTex in a markdown cell, whatever is easiest for you**

1. (0.5) For a medium composed of upper, middle, and lower layers with velocities 6, 8, and 10 km/s, calculate the angle of incidence in the 8 and 10 km/s layers for a ray with an incident angle of 10$^\circ$ in the 6 km/s layer. What is the smallest angle of incidence in the 6 km/s layer that causes total internal reflection at the 8 km/s and 10 km/s interfaces?
 
2. (1.0) Consider two rays that originate from a source at $x=0$ and $z=0$, in a medium with velocity 1 km/s with angles of incidence 0$^\circ$ and 30$^\circ$ (see Figure 2.prob.23 below). Assume that these rays cross an interface at $z=2$ km into a medium with velocity 1.5 km/s and travel to the boundary at $z=4$ km.  For **each** of these ray paths:  
   a. Compute the angle of incidence in the upper layer, the ray path length in each layer, and the total travel time.  
   b. Compute the components and magnitude of the slowness vector $\mathbf{s} = (p, \eta)$ in each layer. Check that the magnitude is related to the velocity as expected.  
   c. Derive the total travel time from the scalar product of slowness and distance ($\mathbf{s}\cdot\mathbf{x}$). Remember to use the appropriate slowness components and horizontal and vertical distances in each layer. Check that these travel times agree with those from (a).

<img src="https://levee.wustl.edu/seismology/book/chapter2/chap2_fr/2_9_p2.jpg" width=700 height=700 class="center"/>
   

3. (1.0) For the `ScSP` conversion at the top of the downgoing slab in Fig. 2.6-15 (below), assume that `ScS` is traveling vertically in the slab, which dips at 30$^\circ$. Assume that the velocities in the slab are $\alpha=9.3$ km/s and $\beta=5.2$ km/s, and the overlying mantle between the slab and the surface has velocities $\alpha_2 = 8.0$ km/s and $\beta_2=4.6$ km/s.  
   a. Find the angle of incidence for `ScS` and `ScSP` at the top of the slab and at the earth's surface.  
   b. Use this result, and the seismograms shown in Fig. 2.6-15 (below) to estimate the depth of the slab. Bear in mind that the `ScSp` and `ScS` arrivals observed at a given station originated from different points on the slab.  

<img src="https://levee.wustl.edu/seismology/book/chapter2/chap2_fr/2_6_15.jpg" width=700 height=700 class="center"/>
 



---------

# Problem 2: Seismic Phases and TauP Background [2.5]

In this problem we will use [ObsPy's TauP](https://docs.obspy.org/packages/obspy.taup.html), a seismic ray tracing and travel-time program based on the [TauP Toolkit](https://github.com/crotwell/TauP) from [Crotwell et al. (1999)](https://pubs.geoscienceworld.org/ssa/srl/article/70/2/154/142385/The-TauP-Toolkit-Flexible-Seismic-Travel-time-and). This question will ask you to read through the TauP Documentation (https://docs.obspy.org/packages/obspy.taup.html) to answer the following questions. 

1. (0.25) The documentation lists a number of available models that define 1D Earth structure. We saw PREM in HW3. Some other models commonly referenced in literature include: Jeffreys-Bullens [jb], [IASP91](https://academic.oup.com/gji/article/105/2/429/705789?login=true), and [AK135](https://academic.oup.com/gji/article/122/1/108/575854?login=true). Referencing the [IASP91 paper (Kennett and Engdahl, 1991)](https://academic.oup.com/gji/article/105/2/429/705789?login=true) (available in the course Google Drive).  
  a. What measurements went into making this velocity model?   
  b. In units of seconds how much did IASP91 change the JB model? (*note that IASP91 came out 50 years after Jeffreys-Bullens, which was made using "hand-cranked calculators"*)  
  c. How was the model tested to ensure if accurately represented Earth structure?   
  d. Figures 7 and 8 in the paper show travel time curves in "reduced travel time". Describe in words how they convert from traveltime to reduced travel time, and what the point of this operation is.  


2. (0.25) The ObsPy TauP documentation provides phase naming convention, which follows mostly from seismological convention. These are also described in Stein and Wysession [SW] Table 3.5-2. The following questions will ask you to describe certain phases based on their phase names. Please use bullet points or numbers to denote each leg, e.g. `Scs` would be.
      > S: upgoing S wave  
      > c: topside CMB reflection  
      > s: upgoing s wave
      >
 Note that one of these phases is not physically possible  
   a. `PKIIKP`  
   b. `pPcPSKKP`  
   c. `SKPiKP`  
   d. Which phase is not physically possible? Why not?    

3. (0.5) One version of IASP91 can be found on the EarthScope (formerly IRIS) [Earth Model Collaboration (EMC)](https://ds.iris.edu/ds/products/emc-earthmodels/), a publicly hosted repository of Earth models. We will plot the model to get an understanding of what it looks like:  
  - Navigate to the IASP91 EMC page: https://ds.iris.edu/ds/products/emc-iasp91/  
  - Download the `IASP91.csv` file to somewhere accessible by this notebook.  
  - Read the file into this notebook (`np.loadtxt` and `np.genfromtxt` are good options, we did something similar in HW1)    
  - Plot the velocity data, trying to re-create the figure on the [IASP91 EMC page](https://ds.iris.edu/ds/products/emc-iasp91/ )  
    - **Important Difference** In your plot, make each data point is visible (e.g., https://matplotlib.org/stable/gallery/pyplots/pyplot_simple.html#sphx-glr-gallery-pyplots-pyplot-simple-py) 
    - Your figure does not have to be pixel for pixel the same, but it should have the same axis orientations etc.  

4. (1.0) Consulting the figure you made in the previous part:
  - If you zoom into your figure, you will see multiple, discrete velocity discontinuities/jumps (2 velocity values occupying the same depth).
  - **Identify** the approximate depth (in km) and associated layer name (e.g., Moho, Core-Mantle Boundary, etc.). You can consult your textbook, the internet, etc. to determine appropriate layer naming.
  - For each layer you have identified, **write a few words describing** why there is a velocity discontinuity here (or why we think there is). Your answer should be in terms of physical processes (i.e., is this a compositional change, mineralogical?).
  - There is one boundary that does not follow this trend right above the Core-Mantle Boundary (better seen in P velocities). **Write a few words describing** what this boundary is and what we think happens here.

5. (0.5) Dividing the earth into 4 distinct zones: 1) crust, 2) mantle, 3) outer core, 4) inner core, and assuming the continental Moho is the representative crust-mantle boundary. Calculate the mean value of $V_p$ and $V_s$ in **each** zone. Print out the values, we will use these for a later problem.




-----------


# Problem 3: TauP Arrival Times [1.5]

Consider the following M8.1 earthquake `11384712` that occurred in the Kermadec trench (https://ds.iris.edu/wilber3/find_stations/11384712), we want to think about the travel times and rays for the seismic phases recorded by seismometer `IU.COLA` in Fairbanks (https://earthquake.usgs.gov/monitoring/operations/stations/IU/COLA/).

> **Notes for this problem:**
> - Truncate all your distances, depths, and coordinates at 2 decimal places. e.g. 177.2791 $\rightarrow$ 177.27
> - The TauP documentation provides example code blocks for running all functions so you should be able to find examples of what you need to input there.


1. (0.25) Determine the distance (in degrees) between source and station, you may find the functions [locations2degrees](https://docs.obspy.org/packages/autogen/obspy.geodetics.base.kilometer2degrees.html) useful (already imported). Using the distance you just calculated, and the real source depth, acquire a list of arrival times for this earthquake using `get_travel_times()` and the `IASP91` model.  Set parameter `phase_list=['ttbasic']`

2. (0.25) Print the arrivals object returned by `get_travel_times()` (e.g., print(arrivals)).
    - Phase `SKS` arrives before phase `S`, why is that?
    - Why is there no phase `PKP`?

3. (0.5) Look at the Wilber3 waveforms for `IU.COLA` (https://ds.iris.edu/wilber3/find_stations/11384712):  
    a. Does Wilber3 use `IASP91` to define their predicted traveltimes?  
    b. Explain whether or not we might expect to see the latest arriving phase in the list of arrivals. Justify your answer based on the timing of the arrival.  

4. (0.5) Consider the great-circle path between the earthquake in Kermadec and receiver in Fairbanks. There is a large-scale "thing" somewhere in the Earth, along or very close to this ray path, that cannot be accounted for by ray theory. We discussed it in lecture and it is also discussed in the relevant chapters in the book.
    - Describe what this "thing" is, and why ray theory cannot account for it.
    - If, somehow, our model could account for this "thing", would your arrival times change? If so, in what direction (faster, slower)?
    - What would we need to include in our computer programs so that we could account for this "thing"? 



-------------

# Problem 4: Raypaths and Traveltimes [3.5]

A picture is worth a thousand words (and a Python plot is worth two thousand). Ray-tracers, TauP notwithsatnding, are great for visualizing ray paths and incidence angles, we can "shoot" rays to any location in our model to understand how they interact with velocity discontinuities.

## Plotting Tips
- If you get a `RuntimeWarning` from Matplotlib while making plots, include `plt.close("all")` at the top of your plotting code cell.
- You can make the TauP figures larger by feeding in a bigger figure object, e.g.,
```python
f, ax = plt.subplots(figsize=(10,10))
plt.axis('off')  # command for hiding the axis so that ObsPy can plot over
# ...
arrivals.plot_rays(fig=f)
```
- Plotting ray paths in Cartesian coordinates may make it easier to visualize shallow or nearby rays
```python
arrivals.plot_rays(plot_type="cartesian")
```

----------

1. Let's set up a similar source-receiver geometry as we had in Problem 2 (i.e., `11384712` $\rightarrow$ `IU.COLA`).
    - Plot the ray paths using `plot_rays()`
    - Set your source depth to 30 km
    - Set your distance to 98$^\circ$
    - Set your phase list to ["P"]  
    - Set `legend=True` in `plot_rays()` so that the phases will be labelled (see example code below).
```python
arrivals.plot_rays(legend=True, phase_list=["P"])
```

2. (1.0) Using **only** the following pieces of information: (1) the ray parameter equation $p=r\sin i/v$, (2) the fact that in TauP the Earth is a sphere with radius=6371 km; (3) any information stored in the `P` ray `Arrival` object (see note below), calculate the following (*show the process of how you attained your answer*):

    a. $\alpha$ (P-wave velocity) in the vicinity of the earthquake  
    b. $\alpha$ in the vicinity of the receiver   
    c. $\alpha$ where the ray bottoms out

```python
# Note: prints out a dictionary of all available Object variables
vars(arrivals[0])
```
4. (1.0) Now you are free to modify the inputs for `phase_list` and `distance_in_degree` (keep source depth at 30km). Play around with `plot_rays()` to identify the following phases. Confirm you found this phase by providing (1) the phase name and (2) a distance at which the ray is a valid solution (valid means TauP will allow you to plot it).  
    a. A crustal head wave that propagates at mantle velocities  
    b. A direct P-wave arriving at outer-core shadow zone distances, that does not reflect and does not enter the mantle    
    c. The `J` phase (*you'll have to look this up*)    
    d. A reflected S-wave that exactly reaches the antipode by circling the earth twice  
    e. An S-wave that enters the outer core, bounces multiple times inside the outer core, does not enter the inner core, and ends up at a distance $<60^\circ$ away from the source  

5. (1.5) Let's take a look at a traveltime curve. You can plot travel times using `plot_travel_times()`.  
- Set your source depth at 0 km
- Set your phase list to be `["P", "PcP", "ScP", "S", "PKiKP", "PP", "PKP", "SKP", "SKS", "SKKS", "PS"]`.  
- Referencing **only** this figure **and** the average layer velocity values you obtained in Problem 2.5 (this problem set), answer the following questions.  
- Please justify your answers, providing only values is not enough, you will need to provide justification based on the travel time curve and velocity values.  
- You may also plot ray paths if it helps you visualize the situation.  
  a. Some of the travel time phases are plotted with negative slopes (they arrive sooner at longer distances). How is that possible?    
  b. Approximately how deep is the core-mantle boundary?  
  c. Approximately how deep is the inner-core boundary?  
  d. A seismometer records a direct P-wave arrival at time $t$, and then a direct S-wave arrival 8 minutes after time $t$. Approximately how far away is the earthquake?  
  e. Between 0 and 20$^\circ$, there are a number of `P` and `S` curves that stop abruptly. Explain why this is.  
  f. Identify a portion of the traveltime curve where triplification occurs (e.g., at *X* degrees and *Y* minutes, the *Z* phase shows triplification).



-----------

## Problem N 
1. Approximately how much time did you spend on this homework assignment?
2. Did you find this homework particularly easy, adequate, or difficult?
3. Any feedback on the homework?



```python

```
