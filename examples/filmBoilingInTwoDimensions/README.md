## Simulation of film boiling in two dimensions

Here the ACDI method is used to simulate film boiling in two dimensions. The computational setup is insprired by [1] but has a smaller wall superheat here (i.e. 1K above saturation temperature). The simulated Nusselt numbers are compared to the emprical correlation from [2].<br>


This setup uses the outflow boundary condition from [3]. This feature is controlled by line 10 in file **`boiling.in`**. Using 'F' as the third argument in line 10, enables the use of the original method in [3].<br>

Details on the computational setup using the ACDI method is found in [4] section 4.6 and in [5] chapter 6.7.<br>

Further and example script **`extract_VTK_data.py`** is used to convert the binary FluTAS files into VTK files for easy sharing of simulation data.<br>

Relavant literature for two timensional film boiling in 2D are:

[1] *Jinliang Xu Dongliang Sun and Qicheng Chen. Modeling of the evaporation and condensation phase-change problems with FLUENT. Numerical Heat Transfer, Part B: Fundamentals, 66(4): 326–342, 2014. doi: 10.1080/10407790.2014.915681. URL https://doi.org/10.1080/10407790.2014.915681*

[2] *V.V. Klimenko. Film boiling on a horizontal plate — new correlation. International Journal of Heat and Mass Transfer, 24(1):69–79, 1981. ISSN 0017-9310. doi: https://doi.org/10.1016/0017-9310(81)90094-6.*

[3] *Akash Dhruv. A vortex damping outflow forcing for multiphase flows with sharp interfacial jumps. Journal of Computational Physics, 511:113122, 2024. ISSN 0021-9991. doi: https://doi.org/10.1016/j.jcp.2024.113122*

[4] *Lorenz Weber, Aritra Mukherjee, Andreas G. Class, and Luca Brandt. A consistent and scalable framework suitable for boiling flows using the conservative diffuse interface method. Journal of Computational Physics, 551:114680, 2026. ISSN 0021-9991. doi: https://doi.org/10.1016/j.jcp.2026.114680*


[5] *Lorenz Weber. Towards the Statistical Evaluation of Phase-Changing Flows: The Development of a Framework for Large-Scale Interface-Resolved Simulations (PhD Thesis), Karlsruher Institut für Technologie (KIT), 02 2026, URL https://doi.org/10.5445/IR/1000190080*

## How to use:
- Copy the FluTAS executable in this folder (may be done by editing and executing the script **`copy_apps.bash`**)
- edit **`dns.in`** to adjust number of processes for MPI execution (line 32)
- edit **`run.sh`** to adjust number of processes for MPI execution or modify file for batch execution on HPC clusters
- edit **`do_postprocess.sh`** to adjust location of virtual enviornment
- execute **`do_postprocess.sh`**
- Nusselt numbers are extracted from the log file **`log.out`** and VTK files created from the binary files in the 'data' folder
- The post process script created the files **`NusseltNumbers.pdf`** and **`BubbleDetachment.pdf`** that visualizes the simulation results

Note that the input file **`scalar.in`** and **`nuc_sites_1.in`** is not required for this case, but kept in the folder to avoid error messages.