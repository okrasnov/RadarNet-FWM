# RadarNet-FWM
## Forward Model of Radar Signals in the Distributed Radar System

This is a simple MATLAB simulator of the set of non-directive radars that are working simultaneously but independently and observes the same set of targets.  

The simulator generate a timeline of coherently processed bursts of sensing waveforms with specified parameters, which are presented as a set of received signals on the range-Doppler plane.  

It consists of two files:

- RadarNet_Sim.m - the simulator. It includes the full description of the
  simulating scene: the details of the radar set configuration, targets setup, including the highway as a moving in a straight line set of objects (can be a few independently moving lanes), and a free moving target that described by the set of waypoints. As soon as all parameters of the scene selected and the simulation timeline (duration and periodicity) are defined, the script generate a sequence of data and write them in spacified mat file.

- visualize_model.m - function that visualizes the input scene's geometry and the time sequence of range-Doppler planes for the radar sequence.

As an example, in '_output_' directory arev presented two records for 1.3 GHz 1W radars and for 9.3GHz 10W radars.

The directory '_doc_' includes the PDF file with a presentation that describe the basic idea of the simulator's design.