# Some HPC projects I worked on at university

In the *MPI Parallelisation* folder is my work parallelising a random field Ising model [1] across multiple threads in a multi core CPU.
In the *3D Spanning Cluster* folder is my work on building a parallelised model to find a spanning cluster in 3D space by merging clusters found concurrently within smaller cubes; see **wifi.c** for the details of how it works. (The build here contains a bug when recombining partitions of the workload)
Reports of all my projects are included, explaining the approach taken and the effective reduction in time.

I plan to add the rest of the work I completed for this module in the future.

[1] Ising model wiki - https://en.wikipedia.org/wiki/Ising_model
