# Model Predictive Control
Self-Driving Car Engineer Nanodegree Program

In this project I implemented a MPC Controller to maneuver the vehicle in autonomous mode around the track. 

![alt text][image1]

## Parameters of MPC

> Model predictive control (MPC) is an advanced method of process control that is used to control a process while satisfying a set of constraints. The main advantage of MPC is the fact that it doesn't allow the current timeslot to be optimized, while keeping future timeslots in account. This is achieved by optimizing a finite time-horizon, but only implementing the current timeslot and then optimizing again

[Wikipedia](https://en.wikipedia.org/wiki/Model_predictive_control)

## Kinematic model

We implemented a kinematic model to control the vehicle. Kinematic model is a simplification of dynamic models that ignore tire forces, gravity and mass.

This model contains:
- x and y coordinates
- orientation angle
- velocity
- cross-track error 
- orientation error 

Actuators are acceleration and steering angle.

![alt text][image2]

## Polynomial fitting

To fit waypoints I use third-order polynomial as it was mentioned in lectures.

Also I transform coordinates of car from world's to local (main.cpp, lines: 117-123)

## Timestep length, duration and latency 

The idea is: we need to maximize N * dt to increase horizont of predictions. But we can't choose very big value for N because calculation will take a lot of time and also we can't choose very big value for dt because in this case you will not apply actuators well-timed.

I made a lot of expreiments and decidet to use **[N:10, dt:0.1]**.

Also in that case dt equals latency (100ms). It allows to fix problem with latency in easy way because we have gap exactly for one step (MPC.cpp, lines: 108-117)

## Improvement 

I'm planning to improve this project by tunning weight more accurate and by adding another parts of cost function.

[//]: # (Image References)
[image1]: ./mpc_car.png
[image2]: ./mpc_equations.png




