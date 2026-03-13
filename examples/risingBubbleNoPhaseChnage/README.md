## Simulation of a rising bubble without phase change

In this test case the ACDI method is used so simulate the dynamics of a three-dimensional rising bubble in the absence of phase change. The results are compare to the reference data obtained from the publication:

[1] *Jutta Adelsberger, Patrick Esser, Michael Griebel, Sven Gross, Margrit Klitz, and Alexander Rüttgers. 3D incompressible two-phase flow benchmark computations for rising droplets. In Proceedings of the 11th World Congress on Computational Mechanics (WCCM XI), 01 2014*

Details on the computational settings are found in Chapter 6.1 in

[2] *Lorenz Weber. Towards the Statistical Evaluation of Phase-Changing Flows: The Development of a Framework for Large-Scale Interface-Resolved Simulations (PhD Thesis), Karlsruher Institut für Technologie (KIT), 02 2026, URL https://doi.org/10.5445/IR/1000190080*

and in section 4.2 in

[3] *Lorenz Weber, Aritra Mukherjee, Andreas G. Class, and Luca Brandt. A consistent and scalable framework suitable for boiling flows using the conservative diffuse interface method. Journal of Computational Physics, 551:114680, 2026. ISSN 0021-9991. doi: https://doi.org/10.1016/j.jcp.2026.114680*

This example also demonstrates the use of the bubble tagging functionality of FluTAS. The bubble is tagged and evaluated during the runtime of the simulation. The fractional bubble information that is distributed across different processors for MPI execution is afterwards concatenated by executing **`concat_tags.bash`** (using the Python script **`concatenate_tags.py`**). In order for the tagging feature to work, the switch 'DO_POSTPROC ?= 1' needs to be set in the file **`FluTAS/src/apps/two_phase_acdi_boiling/app.two_phase_acdi_boiling`**.

## How to use:
- Copy the FluTAS executable in this folder (may be done by editing and executing the script **`copy_apps.bash`**)
- edit **`dns.in`** to adjust number of processes for MPI execution (line 32)
- edit **`run.sh`** to adjust number of processes for MPI execution or modify file for batch execution on HPC clusters
- edit **`concat_tags.bash`** and **`do_postprocess.sh`** to adjust location of virtual enviornment
- execute **`do_postprocess.sh`**
- This creates a file **`Bubble_rise_comparison.pdf`** that compares the results with [1]


Note that the input file **`scalar.in`** is not required for this case, but kept in the folder to avoid error messages.