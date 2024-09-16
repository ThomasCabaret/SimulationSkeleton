# Simulator of a minimalistic "emergence of life" event
- The simulated world exhibits spontaneous formation of self reproducting structures.
- This is a ready code blocks project
- https://www.codeblocks.org/

# Usage of the simulation
* See details in the smfl gui within the program
* At any time, in case of need, the R key removes every particles

# Experiment 1, vesicles formation:
* Click anywhere to start filling the simulation with green particles
* 600 are probably enough (you can see the count in the GUI)
* Click again to stop green particles creation
* Use the A key once, it will introduce a red catalyser particle and start the chain reaction hopefully leading to vesicles formation

# Experiment 2, vesicle reproduction:
* Check "destroy at boundary" in the GUI
* Still in the GUI put "brownian speed to 0"
* Use the S key roughly 30 time to created yellow particles
* You can toggle "centerize" with the C key
Do not forget to disable centerize after, with C key (it's only an help to build a vesicle)
* Use the A key to create a red particle (that will be in the vesicle if this one is centered)
* Click far from the vesicle to start a source of green particle.
If the source looks like a beam it means you forgot to disable centerize with the C key.
(Clicking again stops the green particles production in case of need)
* The vesicle will hopefully feed, grow, and reproduce.

# What's next ?
The best course of action is to migrate this project to a similar system implemented in alien project here:
https://github.com/chrxh/alien because it contains already a lot of features that will be required anyway at some point and a lot of performance optimisations to allow several orders of magnitude larger simulations. Whether or not both projects are merged I think it's best to centralize the communication in the alien-project discord to keep the A-LIFE community together, maybe in a dedicated room: https://discord.gg/xWYByB7sJ4.
After that the steps could be:
* Set up an end to end simulation of vesicle formation and their duplication. In this demonstration there are 2 distincts experiments essentially because vesicles need precursor to reproduce but free catalysts consume precursor and free catalysts are required for vesicule formation. In a nutshell different environnements, or change in environment over time is required. This can be integrated by crafting different areas in the simulated world as well as cyclic change like periodic flow mimicing sea current or tides. The aim is to have vesicle production at some location and some date in the simulation, and reproduction of them at another location later, without human intervention in between. Addition of new chemicals and new reactions might be necessary to mitigate the action of free catalysts in some areas.
* Add new chemicals and new reactions openning on some fitness improvment to see a phase of evolution after the initial viable replicators emergence.
* Add parametrized chemicals and a more multifunction physics of those polymer (kind of similar to protein job in the real world) to transition towards more open darwinian evolution. See what this "parametrized" concept is about here: https://docs.google.com/document/d/1i6MqmgbaOxabFZPrLGpODqssTl1e5lWWDPObg_nZyQo/edit?usp=sharing


# How to build the sources
The build process is mostly automated using the cross-platform CMake build system and the vcpkg package manager.

### Getting the sources
To obtain the sources, please open a command prompt in a suitable directory (which should not contain whitespace characters) and enter the following command:
```
git clone https://github.com/ThomasCabaret/SimulationSkeleton.git
```


### Build instructions (tested on windows only): 
- [Docker desktop](https://www.docker.com/products/docker-desktop/)
- [X server vcxsrv](https://sourceforge.net/projects/vcxsrv/)

1. Open vcxsrv and keep everything default
2. Open Docker desktop
3. Open project in vscode and click `Re-open in container`

Build steps:
```bash
# TODO: Find a way to have this command done in Dockerfile
cd /usr/local/vcpkg && git fetch --unshallow
# Clean up files if previous files exists :
rm -rf  /workspaces/SimulationSkeleton/build/*
cd /workspaces//SimulationSkeleton/ && mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . --config Release -j8
```

Run the program `./ParticleLife`


