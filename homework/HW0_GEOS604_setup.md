# Homework 0
UAF GEOS604 -- Seismology  
Bryant Chow [bhchow@alaska.edu](bhchow@alaska.edu)  

## Semester: Spring 2025
- **Assigned**: 01/14/2025 
- **Due Date and Time**: 01/21/2025 at 14:00:00 (no late penalties will be applied!)
- Total points: 1.0 

## Problem 0: Conda environment setup [1.0 Point]

On your personal computer, or the computer that you will be doing homework on:

>**NOTE** Steps 1&ndash;5 are one-time instructions, once you do this for HW1, you will not have to re-run these commands. Steps 6&ndash;?? will need to be run each time you start a new homework

1. Download and install Miniconda if you don't already have it: https://docs.anaconda.com/miniconda/

2. In a Terminal, clone the course repository to a directory where you will do homework 
    ```bash
        git clone https://github.com/UAFGEOTEACH/GEOS604_seismo.git
        cd GEOS604_seismo
    ```

3. Create a new Conda environment from the `environment.yml` file
    ```bash
        conda create --file=environment.yml
    ```
4. Activate the Conda environment (you will need to do this for all homework)
    ```bash
        conda activate geos604
    ```
5. Change directories into the `homeworks/` directory, you will find the homework assignments here
    ```bash
        cd homework/
    ``` 
6. Jupyter was installed in this environment, start a new Jupyter environment
    ```bash
        jupyter lab
    ```
    This should open up a web browser with the Jupyter interface (if it doesn't click on the link that comes at the end of the startup logs `http://localhost:8888/tree...`) 

7. In the web browser interface, double click the Jupyter notebook for Homework 1.

8. Open the `Text File` app from the launcher and in the created text file, answer the following questions
	1. Approximately how much time did you spend on this homework assignment?
	2. Did you find this homework particularly easy, adequate, or difficult?
	3. Any feedback on the homework?

9. Take a screenshot of your web browser showing the Jupyter Lab interface and the answers your wrote in part (8).

10. Upload your screenshot to your homework submission directory in the class Google Drive.



