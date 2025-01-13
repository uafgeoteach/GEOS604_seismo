# Homework 1: Welcome to Seismology
- **Course:** UAF GEOS604 &mdash; Seismology  
- **Instructor:** Bryant Chow ([bhchow@alaska.edu](bhchow@alaska.edu))
- **Course Website:** [https://bryantchow.com/teaching/geos604](https://bryantchow.com/teaching/geos604) 
- **Last Modified:**  01/12/25

## Semester: Spring 2025
- **Total Points**: 6
- **Assigned**: Thursday Jan. 16, 2025
- **Due Date and Time**: Thursday Jan. 23, 2025 at 14:00:00 AK


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

This intro homework is meant to familiarize you with the Jupyter notebook format that most homework assignements will use, and to give you an intuitive sense of earthquakes, global seismicity, and Alaska seismicity, so that the next time your friend asks you when "The Big One" is happening, you can confidently say, "I don't know".



```python
# Feel free to import any other packages you need here
import numpy as np
import matplotlib.pyplot as plt
```

## Problem 1: Global and Alaska Seismicity [3.0 points] 
- `Wilber 3` ([https://ds.iris.edu/wilber3/find_event](https://ds.iris.edu/wilber3/find_event)) is a nifty Earthscope (formerly IRIS) tool for browsing earthquakes from a web browser.     
- We will be using `Wilber 3` to develop a feeling for global and regional (Alaska) seismicity rates.  
- The point of this Problem is to build some intuition on how many earthquakes are happening around us every day.


1. Navigate to Wikipedia's `Lists of Earthquakes` page and find the list of `Deadliest earthquakes`  
  a. [0.25] Since the year 2000, what earthquake on this list has the lowest magnitude earthquake?  
  b. [0.25] Since the year 2000, what earthquake on this list has the deepest depth?
2. Setting the magnitude you found in (1a) as $M_{min}$, and the depth you found in (1b) as $Z_{max}$  
   a. [0.25] Go to Wilber3 and determine how many earthquakes with $M_{min}<M<10$ and $0<Z<Z_{max}$ have occurred in the last year (starting from 01/01 of this year).  
   b. [0.25] What about in the last 1, 2, 3, 4 and 5 decades?
3. [0.25] Plot the values you found in (2) on a graph of number of earthquakes vs. time. Include (0, 0) as a data point.
4. [0.25] Plot a linear regression of the data plotted in (3). Based on the linear regression, how many earthquakes do we expect per year for this magnitude and depth range?
5. [0.25] Based on the value you found in (4), what is the chance on any given day that we would expect an earthquake that fits in this magnitude and depth range?
6. [0.25] Repeat steps (3) through (4) but only for the last 20 years. How does your "expected earthquakes per year" value change?
7. [0.25] Based on yearly seismicity for the past 3 years, does your answer in (3) or (6) make more sense? Write a few words on how this compares to what you previously expected regarding global seismicity.
8. [0.5] On Wilber 3, define a bounding box around Alaska with longitude $\approx$ [-172, -140] and latitude = [72, 49] (doesn't have to  be exact), repeat steps (2)-(5) for this region.
9. [0.25] Does the value you obtained in 8 differ from your expectation of seismicity in Alaska? Write a few words discussing.


```python

```

---------------------
## Problem 2: Body Waves at a Glance [1.0 point]

Let's build on Problem 1 and develop some general intuition on how far we expect to record seismic signals from events of certain magnitudes. This is useful for understanding what stations can and cannot be used for detailed analyses of earthquakes. 

1. [0.25] Using the same bounding box for Alaska as in Problem 1, search for events with $6.5<M<7.0$ and $0<Z<30$km in the past 10 years. Choose one of the events to use for the rest of the problem (if working together, select different events). What event did you pick? (Click through to the event page of the earthquake and copy paste the red text at the top, e.g,. 2022-01-11 Mww6.6 Fox Islands, Aleutian Islands).
2. On the event page, click on individual stations to see their recording of the event (leave the Networks and Channels as default values, you can adjust the waveform time range using the box in the bottom right of the pop-up window)   
  a. [0.25] At approximately what distance, $D$ (in degrees), can you no longer discern an identifiable P and S arrival?"
> By discernable I mean you can convince yourself you see a pulse at the expected P-wave or S-wave arrival time on any of the 3 channels  
> Only look for direct arrivals P and S, not PP or SS or other phases.

  b. [0.25] Do body waves re-appear at any distance? Can you think of a reason why we no longer see body wave arrivals at these distances?  
  c. [0.25] Include a screenshot of the waveform pop-up window for the furthest station in which you still can discern an identifiable P and S arrival.   


```python

```

------------------
## Problem 3: Surface waves at a Glance [0.5 point]

Let's do the same analysis but for larger magnitude events. Choosing any earthquake in the last years, with $0<Z<10 km$, what is the smallest magnitude earthquake whose **surface waves** are visible at most (or all) available stations? You can make use of the nifty `Show Record Section` button to do this quickly. (Try to make sure you pick an event that has atleast recordings at $\ge150^\circ$ epicentral distance). Include a screenshot of the record section for the event you end up choosing.



```python

```

------
## Problem 4: Earthquake Catalog [1.5  points]

- Now we'll look at events en-masse to get a more quantitative understanding of global seismicity.
- `Wilber 3` has a `Download events` button that generates a text file of all events in a given search.  

1. On Wilber 3, do a global search for events with $M_{min}<M<10$ and $0<Z<Z_{max}$ since the year 2000. Click the `Download events` and save the corresponding text as a `.txt` file.
2. [0.25] Read in the text file (Note: NumPy's `loadtxt` and `genfromtxt` are two options for this), then parse out the origin time, depth, and magnitude into separate lists
4. [0.25] Make a histogram of all magnitudes with bins of size 1. Approximately how many more earthquakes do we have in each bin?
6. [0.25] Now make a histogram of all depths. At what depth do earthquakes primarily occur?
7. You can convert the time strings we have into `datetime` objects and then plot them with Matplotlib (see example code below).  
   a. [0.5] Convert your time array into datetime objects and make a plot of time vs. cumulative magnitude for each magnitude bin (e.g., 6$\le$M$<$7)  
   b. [0.25] Play around with the data and note anything peculiar you observe. Pay attention to areas where the plots deviate from linear, or when plots coincide with one another. 


```python

```

-----------
## Problem N 
1. Approximately how much time did you spend on this homework assignment?
2. Did you find this homework particularly easy, adequate, or difficult?
3. Did you find this homework relevant, useful 
4. Any feedback on the homework?



<!-- your answers here -->
