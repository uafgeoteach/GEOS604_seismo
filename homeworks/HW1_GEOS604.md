# Homework 1: Stress and Strain
UAF GEOS604 -- Seismology  
Bryant Chow [bhchow@alaska.edu](bhchow@alaska.edu)  
Last Tested: DD/MM/YYYY

## Semester: Spring 2025
- **Assigned**: DD/MM/YYYY 
- **Due**: DD/MM/YYYY at 22:00:00 

## Instructions:
1. Homeworks will always be out of 10 points, each problem will say how many points it is worth. 
2. Late homework will be deducted 10\% of the overall grade, so the max total you can get for late homework is 9/10 points
3. The late policy is {\bf not} cumulative, but please try to turn in homework in a timely matter
4. Some problems may have multiple sub-parts that will divide the point total for a given problem
5. Please do the homework in order, some later problems may relate to previous parts
6. You are allowed (and encouraged) to work as partners or in groups, but each individual student must turn in their own assignment, written in their own words (e.g., don't copy-paste code blocks from each other. Your code may look similar but should not look exactly the same!
7. Instructions may sometimes include code snippets that you can copy-paste. Syntax highlighting should denote what you're looking at
8. You can also look at the instruction source code to see
    ```bash
        echo "this is an example terminal command"
    ```
    ```python
        # this is an example Python code block
        from numpy arange, pi, cos

        x = arange(0, 2 * pi, pi / 4)
        y = cos(x)
        
    ```



## Problem 0: Conda environment setup [0.25 Points]

On your personal computer, or the computer that you will be doing homework on:

>**NOTE** Steps 1&ndash;5 are one-time instructions, once you do this for HW1, you will not have to re-run these commands. Steps 6&ndash;?? will need to be run each time you start a new homework

1. [Create a GitHub account](https://github.com/signup?ref_cta=Sign+up&ref_loc=header+logged+out&ref_page=%2F&source=header-home) and email me (bhchow@alaska.edu) your GitHub username

1. Download and install Miniconda if you don't already have it: https://docs.anaconda.com/miniconda/

1. Fork the GEOS604_seismo repository (https://github.com/uafgeoteach/GEOS604_seismo/fork). Keep the default s ettings and hit `Create fork`

1. Clone your Git repository to a directory where you will do homework 
    ```bash
        git clone https://github.com/<YOUR_GITHUB_USERNAME>/GEOS604_seismo.git
        cd GEOS604_seismo
    ```
    >**NOTE**: Be sure to replace `<YOUR_GITHUB_USERNAME>` with your GitHub username to access your own forked repository.

1. Create a new Conda environment from the `environment.yml` file
    ```bash
        conda create --file=environment.yml
    ```
1. Activate the Conda environment (you will need to do this for all homework)
    ```bash
        conda activate geos604
    ```
1. Change directories into the `notebooks/` directory, you will do your homework here
    ```bash
        cd notebooks/
    ``` 
1. Jupyter was installed in this environment, start a new Jupyter environment
    ```bash
        jupyter notebook
    ```
    This should open up a web browser with the Jupyter interface (if it doesn't click on the link that comes at the end of the startup logs `http://localhost:8888/tree...`) 

1. Start a new Jupyter notebook, and name it: 
`last_name`_GEOS604\_HW1.ipynb (where `last_name` is your last name)

1. In a new terminal window (because your old one is now occupied with running the notebook) push your newly created, empty, Jupyter notebook to your cloned repository 

    ```bash
        git add -A  # add the new file to tracking
        git commit -m "first commit!"  # or whatever commit message you wanbt
        git push  # 
    ```

1. Please Complete your homeworks in these Jupyter notebook. If text answers are required please write them in markdown. If you need to show math, you can use latex-style math mode to do so 
    >**NOTE**: If you are not familiar with Latex, this is a great opportunity to learn, it is an incredibly powerful tool for document preparation. But since this is not a class on Latex, you are welcome to hand write math/derivations. Please email me if this is your preferred method and I can provide instructions on how to upload these in conjunction with your Juypter notebooks.


>**PLEASE READ**: It would behoove you to commit and push your work often, so that there is a timestamped, version controlled backup incase your computer fails, or you are unable to complete your assignment. If you tell me you can't turn in your homework because "my dog ate my harddrive!", I will likely take into account the commit history for your given assignment when making any decisions/exceptions.



# Problem N [0.25 Points]

1. Approximately how much time did you spend on this homework assignment?
2. Did you find this homework particularly easy, adequate, or difficult?
3. Any feedback on the homework?



