# CarND-Path-Planning-Project
Self-Driving Car Engineer Nanodegree Program
##### Authon: Hsin-Cheng Chao (olala7846@gmail.com)

##### Note:
This project was forked from Udacity and has been modified heavily, for setup and build detail, please reference the [Udacity's REPO](https://github.com/udacity/CarND-Path-Planning-Project). for the project requirement, please reference the [Rebrics](https://review.udacity.com/#!/rubrics/1020/view)

## Reflection
My implementation was based on the [BMW/Bosch paper](http://ieeexplore.ieee.org/abstract/document/5509799/) and modified to fit the simulator environment.

### System Architecture:
The program communicate to the Udacity simulator through web socket. The simulator sends some messages to the c++ program:

1. The current car state (Global Cartesian coordinate, Frenet Coordinate, car speed ...) 
2. Previouly unexecuted trajectory in Cartesian coordinates.
3. Sensor fusion data (position and speed of other vehicles)

The main.cpp should returns the updated future trajectories in Cartesian coordinate.

#### Conversion between coordinate system.
Most of my work (FSM, trajectory generation) was done under the Frenet coordinate system. So I use the `generate_splines()` function to generate the spline once and store it in memory before the sockets starts, for every waypoint I use six other waypints to generate a spline and use the spline to convert the Frenet coordinate to Cartisian coordinate (See the `getXY` function. 

#### Path Plaining
The Path plaining algorithm was done using a FSM (finite state machine) and JMT (jerk minimization trajectory) with some cost functions. more details below. 

#### Finite State Machine
I created 5 states for the FSM (See the `update_states` function for more details):

1. Velocity Keeping: In this state, the vehicle will try to keep it's speed as high as possible (I use speed_limit * 0.9) as a safe speed. Whenver the car detects vehicle close ahead, it will transition the state to Vehicle Following.
2. Vehicle Following: The car will try to keep a safe distance and speed from the vehicle ahead (5.0 + leading_vehicle_speed * reaction_time) and try to change lane when it's safe to do so.
3. Lane Change Left / Lane Change Right: Car will try to move to the target lane in a safe s and d trajectory.

#### Jerk Minimization Trajectory
I use the jerk minimization algorithm and use the end of the previously trajectory as the start config state of the JMT and generate many possible end states according to the current FSM states.

1. Velocity Keeping / Land change Left / Lane Change Right: The vehicle will try to end go to the target position (30 meters ahead) with different duration and end velocity as end config and choose the one with min cost.

```
// psudo code for sampling different end config
double target_duration = 2.0;
double target_speed = speed_limit * 0.9;
for d in range(target_duration - 1.0, target_duration + 1.0, 0.1)
	for v in range(target_speed - 10.0, target_speed + 10.0, 1.0) {
		end_config = {car_s + 30.0, v, 0.0}
		trajectory = JMT(start_config, end_config, d)
	}
}
```

2. Vehicle Following: In this state the vehicle will try to keep the same speed and a safe distance to the leading vehicle, different from other states, since the end config speed and distance are basically decided, I sample different duration (time to finish) as the end config for the JMT.

```
// psudo code for vehicle following JMT sampling
target_duration = 2.0;
reaction_time = 0.5;
safe_distance = 5.0 + lv_s * reaction_time;  // lv: leading vehicle
for d in range(target_duration - 1.0, d <= target_duration +5.0, 0.5) {
	lv_future_s = lv_s + lv_speed * d  // calculate the position of lv after duration
	end_config = {lv_future_s - safe_distance, lv_speed, 0.0}
	trajectory = JMT(start_config, end_config, d)
}
```

#### Cost function and trajectory validation
Before calculating the cost, I first validate the trajectory speed and acceleration by calculating the derivative of the generate trajectoy. The speed should not break the speed limit and should not go below 0.0, and the acceleration should not exceed the project rubric requirements (10m/s*s). 

After that I calculates the cost of each sampled valid trajectories, by checking the collision, min/max speed, whether it stays in the center of lane ...

```
// psudo code for cost function
min_cost = 1e10;
best_traj = None;
for traj in possible_trajectories:
	cost = 0.0
	cost += 1000 * collision_cost(traj, sensor_fustion)
	cost += 50 * max_min_speed_cost(traj)
	cost += 25 * center_lane_cost(traj)
	if cost < min_cost:
	    min_cost = cost
	    best_traj = traj
```

### Conclusion
The algorithm itself is quite simple, but the coding parts are complicated. There are many edge cases needs to handle (e.g. when the s coordinates when from max_s to 0.0, many thing could break). Also it took me lots of time trial and error trying decide how much end configurations to sample for the JMT. The less trajectory sampled, the less time it take to find the best trajectory (The system will be more responsive but it is more possible to find a less optimal trajectory).
