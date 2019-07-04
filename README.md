# Nociception project

### Firstly, compile mod files
nrnivmodl ./mod

### One tread one C-fiber simulation
nrniv -python onefibersimulation.py 2(number of model)

### Parallel simulation of сomplex activity (several C-fibers)
mpiexec -n 2(number of tread) nrniv -mpi -python parallelsimulation.py 
