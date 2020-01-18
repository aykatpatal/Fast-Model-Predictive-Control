# Fast-Model-Predictive-Control
Fast Model Predictive Control using the Primal Barrier Method


One of the major areas of focus in the field of
autonomous vehicles is safe navigation and control. While there
are a plethora of techniques available for control, one of the most
studied technique is Model Predictive Control (MPC) of vehicle.
An even faster way is the Fast MPC implementation which
improves the speed of standard MPC by exploiting the structure
of the problem that suits the particular problem. In this paper,
we discuss the different aspects of design and implementation
of a fast MPC model on the 1/10th scale F1 race car. A nonlinear
bicycle model is linearized and used. We used
 some techniques that can be incorporated on MPC to
reach real-time constraints by reducing the computation cost
of optimization. Fast MPC provides the suboptimal solution yet
satisfies the inequality constrains. This way we can be sure that
the control strategy will never produce the control input to the
system which can’t be implementable.

The error in x and y was : 

![Alt text](https://github.com/aykatpatal/Fast-Model-Predictive-Control/blob/master/error_x.jpg)
