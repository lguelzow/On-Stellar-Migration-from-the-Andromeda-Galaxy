# On-Stellar-Migration-from-the-Andromeda-Galaxy
Simulation of hypervelocity stars (HVSs) trajectories originating in Andromeda for my paper "On Stellar Migration from the Andromeda Galaxy".

The base of this code is the simulation I developed for my Master's thesis project. For the paper based on it, I improved the simulation.

The purpose of the simulation is to calculate the trajectories of HVSs that are ejected from the Andromeda galaxy. In my project,
the fraction of these stars that come close to the Sun are then compared to data from the GAIA star survey, specifically the GAIA eDR3 catalogue.

At the beginning of each file, there is a short statement about what it does. The main simulation file "HVS-trajectory-simulation.py" has a full table of contents.
In addition, each file has comments explaing the purpose of each step.

Here is a textual description of what each file does:


M31-trajectory.py OR Calculation of Andromeda trajectory

The calculation of the trajectory of Andromeda in the last few billion years is done separately from the main calculation in the file "M31-trajectory.py".
We use the observational values for the present position and velocity of Andromeda from van der Marel (2012) as initial conditions. Since we use 
measurements of Andromeda at present time, the time interval for the calculation runs backwards from $t_0=13.8 Gyr$ to $t=10.0 Gyr$ to enable us to calculate 
the position and velocity of Andromeda in the past.
Next, the differential equations are defined as a function. As previously mentioned, we use the system of ordinary differential equations given by Eq. (\ref{DEs}).
For the acceleration we use what we established in Eq. (\ref{a_and}). The magnitude of the position vector is calculated every time the function is called since 
it is needed for the calculation of the acceleration. After the function is defined, ODEint is called and given the initial conditions, time interval and the 
function containing the DEs. ODEint calls the DE function in every step of the time interval. It integrates the velocity and acceleration components 
calculated from the position and velocity vector from the previous step to find the position and velocity vector at the current time step. For the first step,
it uses the values calculated from the initial vectors.
At the end, \textit{ODEint} returns an array containing the position and velocity vectors of Andromeda at every step of the time interval. The results are written
into the data file "M31-trajectory.txt". The main file of the simulation uses this data file to interpolate functions of the position and velocity of Andromeda 
depending on time.


RNG-initial-conditions.py OR Generation of random initial conditions

In addition to the trajectory of Andromeda, we also generate the initial conditions separately in the file "RNG-initial-conditions.py". At first, we establish
the amount of initial conditions we generate. This determines how long the resulting data file is and how many trajectories are calculated in the simulation.
At the same time, we set the constants and parameters we need for the generation of the HVSs initial conditions. Namely, these are the parameters used in the
Plummer profile (\ref{Plummer}) and the minimum and maximum values of the initial velocity. Given these parameters, we generate random position and velocity
vectors with the steps detailed in Section \ref{Initial}. Since we generate them in spherical coordinates, we subsequently transform them to Cartesian coordinates.
In this form, the initial conditions are written into the data file \texttt{HVS-initial-conditions.txt}. The main simulation file calls this file to read the
HVS initial conditions and adds the position and velocity of Andromeda correct initial conditions in the simulation coordinate systems.


HVS-trajectory-simulation.py OR Simulation of HVS trajectories

The main file of the simulation is "HVS-trajectory-simulation.py". At the beginning, we establish the amount of time steps in the calculation, the constants
discussed in Section \ref{LG} as well as the vectors, rotation matrices and angles discussed in Section \ref{trafo_sec}. Further, we create the data files
the results are written into at the end. In the next section we call on the previously generated data file "M31-trajectory.txt" for the position and
velocity of Andromeda over the last few billions years. We use this data to define time dependent functions for each component of the position and vector of
Andromeda. Following this, the data file "HVS-initial-conditions.txt" containing the initial conditions for the HVSs is called and read to determine
the amount of HVS trajectories that are calculated.
The for-loop that contains solving the differential equations \ref{DEs} is repeated for each trajectory. At first, we define the time interval for
the ODEint function according to the initial send-off time of the HVS.  Subsequently, we calculate the HVS initial position and velocity vectors 
in the simulation coordinate system by adding the position and velocity of Andromeda at the send-off time to the initial conditions read from the data file.
Now, we define the differential equations as a function. This is done similarly as previously described for the trajectory of Andromeda. We use (\ref{a_tot})
for the total acceleration of a HVS in the gravitational potential of both galaxies. The HVS position vector with respect to the centre of each galaxy is
calculated every time the functions is called. They are necessary to calculate the respective acceleration of each galaxy.
As before, ODEint calls the DE function in every time step and finds the position and velocity vectors by integration. After running through the time interval,
ODEint returns the HVS trajectory as a table with the time, position components and velocity components in each time step.
The distance of the HVS to the Milky Way centre as well as the Sun at $t_0$ is used to select trajectories that pass through the Milky Way and close to the Sun.
The selection process is discussed at the beginning of Section \ref{sec:results}. After the first steps of the selection process in the simulation coordinate
frame, we transform the position and velocity vectors of the remaining trajectories at $t_0$ to Galactic coordinates in the heliocentric rest frame with the
transformation established in Section \ref{trafo_sec}. The remaining steps of the selection process are performed to find the HVSs that come within $13 kpc$ of
the Sun.
Finally, we write the the position and velocity vectors of the final HVS selection into the data file "smallest-min-distances.txt". Further, we record secondary
results into two additional data files. The velocity data of the HVS selection within $50 kpc$ of the GC is written into "MW-centre-vel.txt" to analyse their
velocity distribution in comparison to the final HVS selection. Further, we write the distance to the Sun at $t_0$ of every calculated trajectory into 
"distances-to-sun.txt" in order to be able to find the initial conditions of a specific HVS trajectory.

