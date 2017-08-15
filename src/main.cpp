#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include <deque>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "Eigen-3.3/Eigen/Dense"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

const double TIME_STEP = 0.02;
const double MPH_TO_MPS = 0.44704;  // miles per hour to meters per second
const double SPEED_LIMIT = 49.9 * MPH_TO_MPS;
const double MAX_ACC = 9.0;
const double MAX_JERK = 9.0;

enum CarState {
  velocity_keeping,
  vehicle_following
};
CarState current_car_state = velocity_keeping;


// save previous frenet trajectory globally
deque<double> prev_s_traj, prev_d_traj;
vector<deque<double>> prev_frenet_trajectory = {prev_s_traj, prev_d_traj};

vector<tk::spline> global_splines;

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

void generate_splines(vector<double> maps_x, vector<double> maps_y) {
  int num_points = maps_x.size();
  for(int i = 0; i< num_points; i++) {
    // fit spline with 6 points
    // the target point is p2
    int p0 = (i - 2 + num_points) % num_points;
    int p1 = (i - 1 + num_points) % num_points;
    int p2 = (i + num_points) % num_points;
    int p3 = (i + 1 + num_points) % num_points;
    int p4 = (i + 2 + num_points) % num_points;
    int p5 = (i + 3 + num_points) % num_points;

    vector<double> X = {maps_x[p0], maps_x[p1], maps_x[p2], maps_x[p3], maps_x[p4], maps_x[p5]};
    vector<double> Y = {maps_y[p0], maps_y[p1], maps_y[p2], maps_y[p3], maps_y[p4], maps_y[p5]};

    // affine transformation
    double x_shift = X[2];
    double y_shift = Y[2];
    double theta = atan2(Y[3] - Y[2], X[3] - X[2]);

    int num_spline_points = X.size();
    vector<double> _X(num_spline_points), _Y(num_spline_points);
    for(int i = 0; i < num_spline_points; i++) {
      // translate P0 to origin
      double x_t = X[i] - x_shift;
      double y_t = Y[i] - y_shift;
      _X[i] = x_t * cos(-theta) - y_t * sin(-theta);
      _Y[i] = x_t * sin(-theta) + y_t * cos(-theta);
    }

    tk::spline spline_func;
    spline_func.set_points(_X, _Y);
    global_splines.push_back(spline_func);
  }
}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s,
    vector<double> maps_x, vector<double> maps_y, vector<double> maps_dx, vector<double> maps_dy)
{
  int num_points = maps_x.size();
  // should generate splines before getXY;
  assert(num_points == global_splines.size());

	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

  // fit spline
  auto spline_func = global_splines[prev_wp];
  double ratio = (s - maps_s[prev_wp]) / (maps_s[wp2]-maps_s[prev_wp]);

  // Points in car coordinates on prev_wp
  double x0 = maps_x[prev_wp];
  double x1 = maps_x[wp2];
  double y0 = maps_y[prev_wp];
  double y1 = maps_y[wp2];
  double dx = x1 - x0;
  double dy = y1 - y0;
  double theta = atan2(dy, dx);

  double _x = ratio * sqrt(dx * dx + dy * dy);
  double _y = spline_func(_x);

  double x, y;
  // revert affine transformation
  x = x0 + _x * cos(theta) - _y * sin(theta);
  y = y0 + _x * sin(theta) + _y * cos(theta);

  // add d * unit norm vector
  double nx = maps_dx[prev_wp] + ratio * (maps_dx[wp2] - maps_dx[prev_wp]);
  double ny = maps_dy[prev_wp] + ratio * (maps_dy[wp2] - maps_dy[prev_wp]);
  x = x + d * nx;
  y = y + d * ny;

	return {x, y};

}

/* JMT: jerk minimized trajectory */
vector<double> JMT(vector< double> start, vector <double> end, double T)
{
  /*
  Calculate the Jerk Minimizing Trajectory that connects the initial state
  to the final state in time T.

  INPUTS

  start - the vehicles start location given as a length three array
      corresponding to initial values of [s, s_dot, s_double_dot]

  end   - the desired end state for vehicle. Like "start" this is a
      length three array.

  T     - The duration, in seconds, over which this maneuver should occur.

  OUTPUT
  an array of length 6, each value corresponding to a coefficent in the polynomial
  s(t) = a_0 + a_1 * t + a_2 * t**2 + a_3 * t**3 + a_4 * t**4 + a_5 * t**5

  EXAMPLE

  > JMT( [0, 10, 0], [10, 10, 0], 1)
  [0.0, 10.0, 0.0, 0.0, 0.0, 0.0]
  */

  MatrixXd A = MatrixXd(3, 3);
	A << T*T*T, T*T*T*T, T*T*T*T*T,
      3*T*T, 4*T*T*T,5*T*T*T*T,
      6*T, 12*T*T, 20*T*T*T;

	MatrixXd B = MatrixXd(3,1);
	B << end[0]-(start[0]+start[1]*T+.5*start[2]*T*T),
      end[1]-(start[1]+start[2]*T),
      end[2]-start[2];

	MatrixXd Ai = A.inverse();

	MatrixXd C = Ai*B;

	vector <double> result = {start[0], start[1], .5*start[2]};
	for(int i = 0; i < C.size(); i++)
	{
    result.push_back(C.data()[i]);
	}

	return result;
}

int get_lane(double car_d) {
  int lane_id = 0;
  if (4.0 < car_d && car_d <= 8.0) {
    lane_id = 1;
  } else if (8.0 < car_d) {
    lane_id = 2;
  }
  return lane_id;
}

double sensor_fusion_speed(json car_data) {
    int car_id = car_data[0];
    double car_x = car_data[1];
    double car_y = car_data[2];
    double car_vx = car_data[3];
    double car_vy = car_data[4];
    double car_s = car_data[5];
    double car_d = car_data[6];
    return sqrt(car_vx * car_vx + car_vy * car_vy);
}


/* (Olala): update car state machine */
void update_car_state(const vector<double> car, const json &sensor_fusion) {
  double car_s = car[0];
  double car_d = car[1];
  double car_speed = car[2];
  int current_lane = get_lane(car_d);

  int closest_vehicle_idx = -1;
  double closest_distance = 1e5; // large number
  for (int vehicle_idx = 0; vehicle_idx < sensor_fusion.size(); vehicle_idx++) {
    json car_data = sensor_fusion[vehicle_idx];
    int vehicle_id = car_data[0];
    double vehicle_x = car_data[1];
    double vehicle_y = car_data[2];
    double vehicle_vx = car_data[3];
    double vehicle_vy = car_data[4];
    double vehicle_s = car_data[5];
    double vehicle_d = car_data[6];

    // skip car behind or car in different lanes
    if (vehicle_s < car_s || abs(vehicle_d - car_d) > 2.0) {
      continue;
    }
    double distance = vehicle_s - car_s;
    if (distance < closest_distance) {
      closest_vehicle_idx = vehicle_idx;
      closest_distance = distance;
    }
  }

  if (true || closest_vehicle_idx == -1 || closest_distance > 100.0) {
    // std::cout << "velocity keeping\n";
    current_car_state = velocity_keeping;

  } else {
    current_car_state = vehicle_following;
    // std::cout << "vehicle following\n";
  }

}

double poly_eval(double x, vector<double> coeffs) {
  double result = 0.0;
  double t = 1.0;
  for(int i=0; i<coeffs.size(); i++){
    result += coeffs[i] * t;
    t *= x;
  }
  return result;
}

// helper function to get x, x_d, x_dd using the last 3 trajectory
vector<double> get_traj_end_config(deque<double> traj) {
  int traj_size = traj.size();
  assert(traj_size >= 3);
  double x0 = traj[traj_size - 3];
  double x1 = traj[traj_size - 2];
  double x2 = traj[traj_size - 1];

  double v1 = (x1 - x0) / TIME_STEP;
  double v2 = (x2 - x0) / TIME_STEP;

  double a2 = (v2 - v1) / TIME_STEP;
  return {x2, v2, a2};
}

// probability density function
double pdf(double mean, double sigma, double x) {
  double sig2 = sigma * sigma;
  double delta = x - mean;
  return 1.0 / (2.0 * M_PI * sig2) * exp(delta * delta / 2.0 / sig2);
}

// calculate s cost
double s_diff_cost(vector<deque<double>> traj, const vector<double> end_config) {
  auto traj_s = traj[0];
  auto traj_d = traj[1];

  auto traj_end_config = get_traj_end_config(traj_s);
  vector<double> sigma_s = {10.0, 4.0, 2.0};

  double cost = 0.0;
  for(int i = 0 ; i < 3; i++) {
    double actual = traj_end_config[i];
    double expected = end_config[i];
    double sigma = sigma_s[i];
    double diff = abs(expected - actual);
    cost += erf(diff / sigma);
  }

  return cost;
}

// Binary cost funciton for max acceleration
double max_acceleration_cost(vector<deque<double>> traj) {
  auto traj_s = traj[0];
  auto traj_d = traj[1];
  int traj_size = traj_s.size();
  for(int i = 2; i < traj_size; i++) {
    double s2 = traj_s[i];
    double s1 = traj_s[i - 1];
    double s0 = traj_s[i - 2];
    double v1 = (s1 - s0) / TIME_STEP;
    double v2 = (s2 - s1) / TIME_STEP;
    double a2 = (v2 / v1) / TIME_STEP;
    if (a2 > MAX_ACC) {
      std::cout << "max accleration exceed\n";
      return 1.0;
    }
  }
  return 0.0;
}

// Binary cost funciton for car collision detection
double collision_cost(vector<deque<double>> traj, const json &sensor_fusion) {
  auto traj_s = traj[0];
  auto traj_d = traj[1];
  for (int i = 0; i < traj_s.size(); i++) {

    double car_s = traj_s[i];
    double car_d = traj_d[i];
    double delta_t = i * TIME_STEP;

    for (int car_idx = 0; car_idx < sensor_fusion.size(); car_idx++) {
      json car_data = sensor_fusion[car_idx];
      int id = car_data[0];
      double target_x = car_data[1];
      double target_y = car_data[2];
      double target_vx = car_data[3];
      double target_vy = car_data[4];
      double target_s = car_data[5];
      double target_d = car_data[6];

      double target_speed_s = sqrt(target_vx * target_vx + target_vy * target_vy);
      double target_future_s = target_s + target_speed_s * delta_t;
      if (car_d - 3.0 < target_d && target_d < car_d + 3.0 && target_future_s > car_s) {
        double distance = abs(target_future_s - car_s);
        if (distance < 5.0) {
          std::cout << "will kiss\n";
          return 1.0;
        }
      }
    }
  }
  return 0.0;
}

// Binary cost function that checks whether the car breaks the speed limit
double speed_limit_cost(vector<deque<double>> traj) {
  auto traj_s = traj[0];
  auto traj_d = traj[1];
  int traj_size = traj_s.size();
  double SAFE_MARGIN = SPEED_LIMIT * 0.2;
  for(int i = 1; i < traj_size; i++) {
    double s1 = traj_s[i];
    double s0 = traj_s[i - 1];
    double v1 = (s1 - s0) / TIME_STEP;
    // add v1 < 0.0 because we don't want to spot on the highway
    if (v1 >= SPEED_LIMIT) {
      return 1.0;
    }
  }
  return 0.0;
}


vector<vector<double>> very_constraints() {
  vector<vector<double>> all_constraints;
  for(double dt = 0.0; dt <= 2.0; dt += 0.1) {
    for(double dv = -10.0; dv <= 0.0; dv += 1.0) {
      vector<double> a_constraint = {dt, dv};
      all_constraints.push_back(a_constraint);
    }
  }

  // for(double dt = 0.0; dt <= 2.0; dt += 0.1) {
  //   all_constraints.push_back({dt, 0.0});
  // }
  return all_constraints;
}

vector<vector<double>> generate_trajectory(
    const json &car, const json &sensor_fusion,
    vector<double> maps_s, vector<double> maps_x, vector<double> maps_y,
    vector<double> maps_dx, vector<double> maps_dy) {

  const double PREDICT_HORIZON = 2.5;  // time to predict ahead

  auto previous_path_x = car["previous_path_x"];
  auto previous_path_y = car["previous_path_y"];

  // Previous path's end s and d values
  double end_path_s = car["end_path_s"];
  double end_path_d = car["end_path_d"];

  // Save all possible frenet trajectories here cad select
  // the one with min cost
  vector<vector<deque<double>>> frenet_trajectories;

  // setup start config
  double start_s, start_s_d, start_s_dd;
  double start_d, start_d_d, start_d_dd;

  auto prev_s_vals = prev_frenet_trajectory[0];
  auto prev_d_vals = prev_frenet_trajectory[1];

  if (previous_path_x.size() <= 2) {
    std::cout << "generate traj from scratch\n";

    start_s = car["s"];
    start_s_d = car["speed"];
    start_s_dd = 0.0;
    start_d = car["d"];
    start_d_d = 0.0;
    start_d_dd = 0.0;
  } else {
    std::cout << "generate path with previously return traj\n";

    // remove trajectories that already been comsumed
    while (prev_s_vals.size() > previous_path_x.size()) {
      prev_s_vals.pop_front();
      prev_d_vals.pop_front();
    }

    // TODO(Olala): save the previously end config directly
    assert(prev_s_vals.size() == previous_path_x.size());
    int prev_step_size = prev_s_vals.size();
    double s0 = prev_s_vals[prev_step_size - 3];
    double s1 = prev_s_vals[prev_step_size - 2];
    double s2 = prev_s_vals[prev_step_size - 1];
    double d0 = prev_d_vals[prev_step_size - 3];
    double d1 = prev_d_vals[prev_step_size - 2];
    double d2 = prev_d_vals[prev_step_size - 1];
    double vs1 = (s1 - s0) / TIME_STEP;
    double vs2 = (s2 - s1) / TIME_STEP;
    double vd1 = (d1 - d0) / TIME_STEP;
    double vd2 = (d2 - d1) / TIME_STEP;
    double as2 = (vs2 - vs1) / TIME_STEP;
    double ad2 = (vd2 - vd1) / TIME_STEP;

    start_s = s2;
    start_s_d = vs2;
    start_s_dd = as2;
    start_d = d2;
    start_d_d = vd2;
    start_d_dd = ad2;
  }
  vector<double> start_s_config = {start_s, start_s_d, start_s_dd};
  vector<double> start_d_config = {start_d, start_d_d, start_d_dd};
  // std::cout << "Start config:" << start_s << "," << start_s_d << ", " << start_s_dd << std::endl;

  // Setup current target config
  double time_to_goal = 3.0;
  double current_speed = start_s_d;
  double target_speed = min(SPEED_LIMIT, (current_speed + MAX_ACC * time_to_goal));
  // double target_s = start_s + 0.5 * (current_speed + target_speed) * time_to_goal * 0.5;
  // TODO(Olala): what is the correct target_s?
  double target_s = start_s + 50;
  int target_lane = 1;
  double target_d = 2.0 + 4.0 * target_lane;
  std::cout << "target:" << target_s << "," << target_speed << "," << time_to_goal << std::endl;

  // Try different end configs with very_constraints()
  auto all_constraints = very_constraints();
  for (int c_idx = 0; c_idx < all_constraints.size(); c_idx++) {

    auto a_constraint = all_constraints[c_idx];
    double dt = a_constraint[0];
    double dv = a_constraint[1];

    vector<double> try_end_s = {target_s, target_speed + dv, 0.0};
    vector<double> try_end_d = {target_d, 0.0, 0.0};
    double try_time_to_goal = time_to_goal + dt;

    deque<double> next_s_vals;
    deque<double> next_d_vals;

    // Use some previous trajectories left
    for (int i = 0; i < prev_s_vals.size(); i++) {
      next_s_vals.push_back(prev_s_vals[i]);
      next_d_vals.push_back(prev_d_vals[i]);
    }

    // JMT
    auto s_trajectory_coeff = JMT(start_s_config, try_end_s, try_time_to_goal);
    auto d_trajectory_coeff = JMT(start_d_config, try_end_d, try_time_to_goal);

    int prev_step_size = prev_s_vals.size();
    int total_stpes = PREDICT_HORIZON / TIME_STEP;

    for (int i = 0; (i + prev_step_size) <= total_stpes; i++) {
      double t = (i + 1) * TIME_STEP;
      double s = poly_eval(t, s_trajectory_coeff);
      double d = poly_eval(t, d_trajectory_coeff);

      next_s_vals.push_back(s);
      next_d_vals.push_back(d);
    }

    vector<deque<double>> a_trajectory = {next_s_vals, next_d_vals};
    frenet_trajectories.push_back(a_trajectory);
  }

  // Evaluate trajectory cost
  int best_trajectory_idx = 0;
  double min_cost = 1e10;
  for (int i = 0; i < frenet_trajectories.size(); i++) {
    auto trajectory = frenet_trajectories[i];
    double cost = 0.0;
    cost += 1000.0 * collision_cost(trajectory, sensor_fusion);
    cost += 1000.0 * max_acceleration_cost(trajectory);
    cost += 500.0 * speed_limit_cost(trajectory);
    // cost += 1.0 * s_diff_cost(trajectory, end_s);

    if(cost < min_cost) {
      best_trajectory_idx = i;
      min_cost = cost;
    }
  }

  std::cout << "best traj:" << best_trajectory_idx << " cost:" << min_cost << std::endl;
  auto best_trajectory = frenet_trajectories[best_trajectory_idx];

  // save frenet trajectory for next time
  auto next_s_vals = best_trajectory[0];
  auto next_d_vals = best_trajectory[1];
  prev_frenet_trajectory.clear();
  prev_frenet_trajectory.push_back(next_s_vals);
  prev_frenet_trajectory.push_back(next_d_vals);

  // Convert Frenet coordinate to Cartesian coordinate
  vector<double> next_x_vals;
  vector<double> next_y_vals;

  // std::cout << "frenet to cartesian\n";
  for(int i=0; i < next_s_vals.size(); i++) {
    double s = next_s_vals[i];
    double d = next_d_vals[i];
    // std::cout << "Frenet[" << i << "]:" << s << "," << d << "\n";
    vector<double> point = getXY(s, d, maps_s, maps_x, maps_y, maps_dx, maps_dy);
    next_x_vals.push_back(point[0]);
    next_y_vals.push_back(point[1]);
    // std::cout << "Cartesian[" << i << "] " << point[0] << ", " << point[1] << "\n";
  }
  vector<vector<double>> trajectory = {next_x_vals, next_y_vals};

  return trajectory;
}


int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  // Generate spline and save it for later use
  generate_splines(map_waypoints_x, map_waypoints_y);

  h.onMessage([
      &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](
        uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);

        string event = j[0].get<string>();

        if (event == "telemetry") {
          // j[1] is the data JSON object

        	// Main car's localization Data
          double car_x = j[1]["x"];
          double car_y = j[1]["y"];
          double car_s = j[1]["s"];
          double car_d = j[1]["d"];
          double car_yaw = j[1]["yaw"];
          double car_speed = j[1]["speed"];

          // Previous path data given to the Planner
          auto previous_path_x = j[1]["previous_path_x"];
          auto previous_path_y = j[1]["previous_path_y"];

          // Previous path's end s and d values
          double end_path_s = j[1]["end_path_s"];
          double end_path_d = j[1]["end_path_d"];

          // Sensor Fusion Data, a list of all other cars on the same side of the road.
          auto sensor_fusion = j[1]["sensor_fusion"];

          json msgJson;

          vector<double> next_x_vals;
          vector<double> next_y_vals;

          // (Olala): My code below

          // Behavior plaining
          vector<double> car = {car_s, car_d, car_speed};
          update_car_state(car, sensor_fusion);


          // Trajectory Genereation
          auto best_trajectory = generate_trajectory(
              j[1], sensor_fusion,
              map_waypoints_s, map_waypoints_x, map_waypoints_y,
              map_waypoints_dx, map_waypoints_dy);

          std::cout << "traj_size:" << best_trajectory[0].size() << std::endl;

          next_x_vals = best_trajectory[0];
          next_y_vals = best_trajectory[1];

          // TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;

          auto msg = "42[\"control\","+ msgJson.dump()+"]";

          //this_thread::sleep_for(chrono::milliseconds(1000));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}

