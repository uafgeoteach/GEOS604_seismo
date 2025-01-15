## Startup Procedures for GEOS604 Homework

These instructions are used anything you begin a new homework assignment 
(except for HW0) for GEOS 604. Please make sure you have followed the setup 
instructions in HW0 to install Conda and setup for Conda environment.

1. Open your Terminal and navigate to the directory where this repository is:
```bash
cd <path_to_GEOS604_repo>
```
2. Run `git pull` to download the latest repository will contain the most up to
   date version of the homeworks
```bash
git pull
```
3. Copy the weeks homework assignment into the directory `my_homework/`. (Note: 
   This will prevent you from updating the master copy in homework/, and 
   hopefully keep Git from complaining about differing versions)
```bash
# e.g., if you are doing HW1
cp -r homework/HW1*.ipynb my_homework/
```
4. Optional: It's always a good idea to back up your work. The simplest way is 
   to make a copy of `my_homework/` anytime you make significant progress, e.g.,
```bash
cp -r my_homework/ my_homework_backup_1-14-24
```
  You may also choose to backup in other ways. If you are ever asked to overwrite
  a file, think carefully about what you are about to do!
5. Change directory into the directory where you'll be doing homework 
   (my\_homework/)
```bash
cd my_homework/
```
6. Activate the course Conda environment
```bash
conda activate geos604
```
7. Start Jupyter Lab
```bash
jupyter lab
```
8. The Jupyter interface will start up and open a new browser window. If you 
   do not see a new window, copy-paste the link at the bottom that starts with
   `http://localhost:8888/tree...`
9. On the navigation bar on the left, double-click on the homework you just copied
10. Proceed with the Jupyter notebook. Remember to save often! And create backups
    as needed.
