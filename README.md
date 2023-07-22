# Particle-Life
- This is a ready code blocks project
- https://www.codeblocks.org/


At any time, in case of need, the R key removes every particles

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
https://github.com/chrxh/alien because it contains already a lot of features that will be required anyway at some point and a lot of performance optimisations to allow several orders of magnitude larger simulations.
After that the steps could be:
* Set up an end to end simulation of vesicle formation and their duplication. In this demonstration there are 2 distincts experiments essentially because vesicles need precursor to reproduce but free catalysts consume precursor and free catalysts are required for vesicule formation. In a nutshell different environnements, or change in environment over time is required. This can be implemented but crafting different areas in the simulated world as well as cyclic change like periodic flow mimicing sea current or tides. The aim is to have production vesicle at some location and some date in the simulation, and reproduction of them at another location later, without human intervention. Addition of new chemicals and new reactions might be necessary to mitigate the action of free catalists in some areas.
* Add new chemicals and new reaction openning on some fitness improvment to see a phase of evolution after the initial viable replicators emergence.
* Add parametrized chemicals and a more multifunction physics of those to transition towards more open darwinian evolution. See what this "parametrized" concept is about here: https://docs.google.com/document/d/1i6MqmgbaOxabFZPrLGpODqssTl1e5lWWDPObg_nZyQo/edit?usp=sharing
