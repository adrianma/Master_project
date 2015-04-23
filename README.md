# Master project
Adrian Martinez Gomez (adrianma@student.ethz.ch)
April 2015
#
This are the files for my Master thesis project.
The Matlab codes used during this project are to find under the folder "Codes".
The subfolder are to understand as follows:

* ADRIAN - main scripts and functions developed by myself during my project.
Here there are a couple of further subfolders and numerous scripts to perform
System ID as well as Control.
* archived_data - generated ARMAX models and simulations.
Note that a lot of the simulations were left out due to their big file size. However you can find the models and processed iddata.
* BLKTRIDIAG - helping function that can generate tridiagonal matrices (not written by me)
* VAGGELIS - main simulation of the EWH population. Contains some speed-up and vectorizatoin done by me, to accelerate the simulations. To further accelerate them generate a MEX file for the function "mixing.m" and edit the function "buoyCor.m" accordingly (as this is a very system-specific task it is left to the reader)

Now some instructions follow to perform some of the tasks from the report:

(1) Open-loop autonomous simulations:
This is all contained into the function "main_No_Ext_Control"
You have to basically generate the population; then construct the water draws; and lastly simulate with the selected number of EWHs n_app
(Note that you have to generate these files, since they are not provided)  
An example function call:
% number of EWHs and number of simulations
n_app = 1000; % EWHs
N_experiments = 20; % simulations
main_No_Ext_Control(n_app,N_experiments)

(2) Generate Set-point variation experiments for ARMAX System ID:
Use the file "main_experiments" (which also supports Probability Switching):
% Example (copy-paste these 4 lines into your Matlab prompt)
% get the data for 1000EWHs at hour 16
my_hour = 16;
main_experiments('SetPointVariation',my_hour,1000,1);

(3) Process the Data for System ID:
Use the function "Data_system_ID.m", where the call could be for example as follows:
% substract the averaged baseline for output
% choose as exogenous input the averaged water draws
my_Data = Data_system_ID('SetPointVariation',[3,9],1000,2,3,1);

(4) Assess the models (only works for data of hours 4,7,12,16,18,20,22):
Use the file "Assess_houly.m": This function splits the data into their hours, separates by sign and generates the ARMAX (or state-space) models.

Note: the models computed by me and presented in Table 4.1 of the written thesis are stored in "archived_data/SPV/Data_same/15032015_same_SPV.mat"

(5) Closed loop MPC control:
(This assumes that you have the GUROBI solver installed. Please see www.gurobi.com for installation)
I generated a script that runs the examples for artifical signals for the MPC controller for the given parameters of the thesis. You just simply have to run the script "experiments_ClosedLoop.m"
Take into account that there are LOTS of options to be added into the closed loop simulation (saving file options, changing the cost function, changing horizon, etc). 
It should be all documented in that script or within the file "main_CL_simulation".
