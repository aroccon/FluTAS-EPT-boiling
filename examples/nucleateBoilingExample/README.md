## Simulating nucleate boiling in a turbulent subcooled flow

In this simulation, nuclate boiling in a turbulent flow is simulated on a heated wall. The physical properties and the flow conditions are taken from [1]. For the turbulent inflow bondary condition the file **`sem_input_1.txt`** has to be present. This file enables the use of the divergence free synthetic eddy method [2] to produce isotropic turbulence. Within this file, the mean flow in stream-wise direction, the turbulent kintik energy and the dissipation rate - all as a function of the wall distance - is specified as input. Here, the turbulent inflow is located at the lower boundary in the y-direction. The turbulent input data is taken from [3], since a ductflow was used in [1] too.   <br>

The bubble nucleation is triggered by using random nucleation sites based on the algorithm proposed in [4]. The locations and activation temperatures of each nucleation site is listed in the file **`randomNucleationSites.in`**. The adaptation of this method for the use with the ACDI method is described in chapter 7.2 [6].

The surface of the solid - fluid interface is heated. The spatial extent and power of the heater is specified in **`ibm.in`** (line 6). More on the information regarding the implementation is found in chapter 7.4 in [6].<br>

This setup uses the outflow boundary condition from [6]. This feature is controlled by line 10 in file **`boiling.in`**. Using 'F' as the third argument in line 10, enables the use of the method in [6] chapter 4.6.<br>

Simplified LES model was implemented (see [5] and chapter 7.5 in [6]) to model the subgrid scale turbulence with a standard Smagorinsky model and Van Driest wall dampening. This feature is controlled by line 12 in the file **`boiling.in`**.<br>

Relavant literature for this simulation:

[1] *Artyom Kossolapov. Experimental Investigation of Subcooled Flow Boiling and CHF at Proto typical Pressures of Light Water Reactors. PhD thesis, Massachusetts Institute of Technology. Department of Nuclear Science and Engineering, 9 2021. URL* https://hdl.handle.net/1721.1/151891

[2] *Ruggero Poletto. Divergence free development of the sythetic eddy method in oder to improve syn thetic turbulence for embedded LES simulations. PhD thesis, School of Mechanical, Aerospace and Civil Engineering, The University of Manchester, 3 2015*

[3] *Sergio Pirozzoli, Davide Modesti, Paolo Orlandi, and Francesco Grasso. Turbulence and secondary motions in square duct flow. Journal of Fluid Mechanics, 840:631–655, 2018. doi:10.1017/jfm.2018.66*

[4] *Yohei Sato and Bojan Niceno. Nucleate pool boiling simulations using the interface tracking method: Boilingregimefromdiscretebubbletovapormushroomregion. InternationalJournal of Heat and Mass Transfer, 105:505–524, 2017. ISSN 0017-9310. doi:* https://doi.org/10.1016/j.ijheatmasstransfer.2016.10.018

[5] *Sato, Y., Kossolapov, A., Niceno, B., & Bucci, M. (2026). Interface tracking simulation for subcooled flow boiling of water at 10 bar. Journal of Fluid Mechanics, 1030, A40. doi:10.1017/jfm.2026.11253*

[6] *Lorenz Weber. Towards the Statistical Evaluation of Phase-Changing Flows: The Development of a Framework for Large-Scale Interface-Resolved Simulations (PhD Thesis), Karlsruher Institut für Technologie (KIT), 02 2026, URL https://doi.org/10.5445/IR/1000190080*

## How to use:
- Copy the FluTAS executable in this folder (may be done by editing and executing the script **`copy_apps.bash`**)
- edit **`dns.in`** to adjust number of processes for MPI execution (line 32)
- edit **`run.sh`** to adjust number of processes for MPI execution or modify file for batch execution on HPC clusters
- edit **`wallTemps.bash`** to adjust location of virtual enviornment for Python
- execute **`wallTemps.bash`**
- Maximum wall temperatures for each time step are extracted from the log file **`log.out`**
- The post process script created a file **`Evolution.pdf`** that visualizes the evolution of the max. wall temperature over the simulation runtime as an indicator for wether the simulation has reached its steady state.

Note that the input file **`scalar.in`** is not required for this case, but kept in the folder to avoid error messages.