# Analysis

This page describes the different aspects of the MPC implementation for path following in this project, which is the 5th project in the 2nd term of the Udacity Self Driving Car Nanodegree.

The names of the headings below correspond to the criterea under the 'Implementation' part of the rubric.

## The Model

The state follows that introduced in the lectures. This includes x, y, psi, v, cte, and epsi, which are described in the table below:


| Variable | Description |
|:--------:|:-----------:|
| x        | The x position of the vehicle reference point, defined relative to the ego vehicle coordinate frame at the start of each control iteration. |
| y        | The y position of the vehicle reference point, defined relative to the ego vehicle coordinate frame at the start of each control iteration. |
| psi      | The angle measured between the ego vehicle longitudinal axis and the x axis of the ego vehicle coordinate frame at the start of each control iteration. Positive for counter clock-wise rotations when looking down at the xy plane from a coordinate with a positive z component (the coordinate system follows the RHR, i.e. x cross y is z) |
| v        | The speed of the vehicle reference point relative to the ground. |
| cte      | The cross track error of the planned path relative to the reference path |
| epsi     | The heading error of the planned path relative to the reference path |

The update equations are:

let k = j+1

x_k = x_j + v_j * sin(psi_j)
y_k = y_j + v_j * sin(psi_j)
psi_k = psi_j + v_j * (delta_j / LF) * dt
v_k = v_j + a_j * dt
cte_k = (f_j - y_j) + (v_j * CppAD::sin(epsi_j) * dt)
epsi_k = (psi_j - psides_j) + v0 * delta0 / Lf * dt

Note that the fg vector is filled using the constraints, which for x is x_k - (x_j + v_j * sin(psi_j)). The others follow a similar pattern.

## Timestep Length and Elasped Duration (N & dt)

I started with N 100 and dt 0.1, but found that this was projecting too far into the distance. I settled on N 100 and dt 0.1, since this provided a reasonable lookahead distance and planning resolution. Also, 10 prediction points is much easier to debug using text printouts than 30 or even 20 prediction points. I also wanted to ensure not to project too far because the projection distance is proportional to speed. In case I wanted to adjust the speed to 50 or 60 mph, I didn't want to have any problem with fitting a polynomial to a data set that was not one to one, or with projected too far so that the model becomes too inaccurate.

Note that in general, we want to choose values that are as small as possible for dt and as large as possible for N, such that the desired horizon is covered and we do not induce too much compuatational load.

## Polynomial Fitting and MPC Preprocessing

A third order polynomial is fitted to the waypoints to enable easier calculation of CTE and Heading Error.

No modification is made to the vehicle state, actuators, or waypoints prior to applying the MPC, other than the adjustment for Latency described in the next section.

After MPC, the steering angle is adjusted from an angle to a range [-1, +1] by dividing by the maximum possible angle of +/-25 degrees. The acceleration value is also converted to a throttle value by multiplying by a scalar based on data collected near the indended speed limit. This is in addition to the cost based adjustments that happen in the MPC cost calculation step.

## Model Predictive Control with Latency

In order to compensate for the 100 ms latency, I first made the controller work without the latency. Upon adding the latency back in, I attempted to adjust how the calculation of the constraints was done in FG_eval, but this did not seem to work well (I used the inputs from 2 time steps ago instead of 1 to calculate the values for the current time step). This did not seem to work well, so I decided to try the suggestion from the course, which was to advance the state using the model by the latency amount before processing using the MPC. This did work quite well, although the performance was still worse that the system without the latency injected (which is to be expected, since the model does not predict perfectly what will happen).