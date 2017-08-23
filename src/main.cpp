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
const double SPEED_LIMIT = 49.5 * MPH_TO_MPS;
const double MAX_ACC = 9.0;
const double MAX_JERK = 9.0;

// The max s value before wrapping around the track back to 0
const double max_s = 6945.554;
double half_max_s = max_s / 2.0;

enum CarState {
  velocity_keeping,
  vehicle_following,
  lane_change_left,
  lane_change_right,
};
CarState current_car_state = velocity_keeping;


// save previous frenet trajectory globally
deque<double> prev_s_traj, prev_d_traj;
vector<deque<double>> prev_frenet_trajectory = {prev_s_traj, prev_d_traj};

vector<tk::spline> global_splines;
int lane_change_target_lane = 1;

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

vector<double> derivative(vector<double> coeff) {
  vector<double> derivative_coeffs;
  for (int i = 1; i < coeff.size(); i++) {
    derivative_coeffs.push_back(i * coeff[i]);
  }
  return derivative_coeffs;
}

double s_dist(double car_s, double target_s) {
  if (target_s < (car_s - half_max_s)) target_s += max_s;
  if (target_s > (car_s + half_max_s)) target_s -= max_s;
  return target_s - car_s;
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
  assert(0.0 <= s && s <= max_s);
  int num_points = maps_x.size();
  // should generate splines before getXY;
  assert(num_points == global_splines.size());

	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}
  // handle last waypoin
  if (prev_wp == -1) {
    prev_wp = maps_x.size() - 1;
  }

	int wp2 = (prev_wp+1)%maps_x.size();

  // fit spline
  auto spline_func = global_splines[prev_wp];
  // handle next_wp == 0 (s[0] will be 0.0)
  double next_wp_s = maps_s[wp2];
  if(next_wp_s == 0.0) {
    next_wp_s = max_s;
  }
  double ratio = (s - maps_s[prev_wp]) / (next_wp_s - maps_s[prev_wp]);

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

int get_leading_vehicle_id(double car_s, double car_d, const json &sensor_fusion) {
  int current_lane = get_lane(car_d);

  int leading_vehicle_id = -1;
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

    // distance to target, positive means ahead, negative means behind
    car_s = fmod(car_s, max_s);
    double distance = s_dist(car_s, vehicle_s);

    if (distance < 0 || abs(vehicle_d - car_d) > 2.5) {
      continue;
    }
    if (distance < 100.0 && distance < closest_distance) {
      leading_vehicle_id = vehicle_id;
      closest_distance = distance;
    }
  }
  return leading_vehicle_id;
}

json get_leading_vehicle_by_id(int id, const json &sensor_fusion) {
  for (int i = 0; i < sensor_fusion.size(); i++) {
    json car_data = sensor_fusion[i];
    if (car_data[0] == id) {
      return car_data;
    }
  }
  assert(false); // vehicle id not found
}

/* (Olala): update car state machine */
void update_car_state(const vector<double> car, const json &sensor_fusion) {

  double car_s = car[0];
  double car_d = car[1];
  double car_speed = car[2];
  int current_lane = get_lane(car_d);

  int leading_vehicle_id = get_leading_vehicle_id(car_s, car_d, sensor_fusion);

  json leading_vehicle_data;
  double lv_s, lv_d, lv_dist;
  if (leading_vehicle_id != -1) {
    leading_vehicle_data = get_leading_vehicle_by_id(leading_vehicle_id, sensor_fusion);
    lv_s = leading_vehicle_data[5];
    lv_d = leading_vehicle_data[6];
    lv_dist = lv_s - car_s;
  }

  double KEEP_SPEED_DIST = 50.0;
  bool left_is_clear = true;
  bool right_is_clear = true;
  for (int vid = 0; vid < sensor_fusion.size(); vid++) {
    json car_data = sensor_fusion[vid];
    int vehicle_id = car_data[0];
    double vehicle_x = car_data[1];
    double vehicle_y = car_data[2];
    double vehicle_vx = car_data[3];
    double vehicle_vy = car_data[4];
    double vehicle_s = car_data[5];
    double vehicle_d = car_data[6];
    double vehicle_dist = s_dist(car_s, vehicle_s);

    int vehicle_lane = get_lane(vehicle_d);
    if (-20.0 <= vehicle_dist && vehicle_dist <= (KEEP_SPEED_DIST + 10.0)) {
      if (vehicle_lane == (current_lane - 1))
        left_is_clear = false;
      if (vehicle_lane == (current_lane + 1))
        right_is_clear = false;
    }
    if (current_lane == 0)
      left_is_clear = false;
    if (current_lane == 2)
      right_is_clear = false;
  }

  if (current_car_state == velocity_keeping) {
    double target_d = 2.0 + current_lane * 4.0;
    double target_d_diff = abs(target_d - car_d);
    if (leading_vehicle_id == -1 || lv_dist > KEEP_SPEED_DIST) {
      // keep current state
    } else {
      current_car_state = vehicle_following;
    }
  } else if (current_car_state == vehicle_following) {
    double target_d = 2.0 + current_lane * 4.0;
    double target_d_diff = abs(target_d - car_d);
    if (leading_vehicle_id == -1 || lv_dist > (KEEP_SPEED_DIST + 10.0)) {
      current_car_state = velocity_keeping;
    } else if(target_d_diff > 0.1) {
      // don't change lane before cetering current lane
      std::cout << "Keep in center of lane\n";
    }else if (left_is_clear) {
      current_car_state = lane_change_left;
      lane_change_target_lane = current_lane - 1;
    } else if (right_is_clear) {
      current_car_state = lane_change_right;
      lane_change_target_lane = current_lane + 1;
    }
  } else if (current_car_state == lane_change_left) {
    double target_d = 2.0 + lane_change_target_lane * 4.0;
    double target_d_diff = abs(target_d - car_d);
    if (target_d_diff >= 0.1) {
      // center current lane
    } else {
      current_car_state = velocity_keeping;
    }
  } else if (current_car_state == lane_change_right) {
    double target_d = 2.0 + lane_change_target_lane * 4.0;
    double target_d_diff = abs(target_d - car_d);
    if (target_d_diff >= 0.1) {
      // cener curerrent lane
    } else {
      current_car_state = velocity_keeping;
    }
  }
  std::cout << "current state:" << current_car_state << " target lane:" << lane_change_target_lane<< std::endl;
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

double collision_cost(vector<deque<double>> traj, vector<vector<vector<double>>> predictions) {
  auto traj_s = traj[0];
  auto traj_d = traj[1];
  assert(predictions[0].size() == traj_s.size());

  for (int v = 0; v < predictions.size(); v++) {
    auto target_prediciton = predictions[v];
    for (int i = 0; i < traj_s.size(); i++) {
      // our car is named "eagle"
      double eagle_s = traj_s[i];
      double eagle_d = traj_d[i];

      auto target_pos = target_prediciton[i];
      auto target_s = target_pos[0];
      auto target_d = target_pos[1];

      if (abs(eagle_d - target_d) < 3.0) {
        double s_distance = s_dist(eagle_s, target_s);
        if (abs(s_distance) < 5.0) { // will collide
          // std::cout << "Will kiss: \n";
          // std::cout << "i:" << i << ", eagle_s:" << eagle_s << ", eagle_d" << eagle_d;
          // std::cout << ", target_s:" << target_s << ", target_d" << target_d << std::endl;
          return 1.0;
        }
      }
    }
  }
  return 0.0;
}

// Binary cost function that checks whether the car breaks the speed limit
double speed_range_cost(vector<deque<double>> traj) {
  auto traj_s = traj[0];
  auto traj_d = traj[1];
  int traj_size = traj_s.size();

  double min_speed = SPEED_LIMIT;
  double max_speed = 0.0;

  for(int i = 1; i < traj_size; i++) {
    double s1 = traj_s[i];
    double s0 = traj_s[i - 1];
    double v = (s1 - s0) / TIME_STEP;

    min_speed = min_speed > v? v: min_speed;
    max_speed = max_speed < v? v: max_speed;
  }

  double TARGET_RATIO = 0.9;
  double target_speed = SPEED_LIMIT * TARGET_RATIO;
  double min_speed_cost = 0.0;
  double max_speed_cost = 0.0;
  if (min_speed < target_speed) {
    min_speed_cost = 1.0 - min_speed / target_speed;
  }

  if (max_speed > target_speed) {
    max_speed_cost = (max_speed - target_speed) / (SPEED_LIMIT - target_speed);
  }

  return (max_speed_cost + min_speed_cost) / 2.0;
}

double d_offset_cost(vector<deque<double>> traj) {
  auto traj_s = traj[0];
  auto traj_d = traj[1];
  int traj_size = traj_s.size();
  double total_offset = 0.0;
  for (int i = 0; i < traj_size; i++) {
    double d = traj_d[i];
    total_offset += abs(fmod(d, 4.0) - 2.0) / 2.0;
  }
  return total_offset / traj_size;
}

vector<vector<vector<double>>> enumerate_coeffs_combs(
    vector<vector<double>> possible_s_coeffs,
    vector<vector<double>> possible_d_coeffs) {

  std::cout << "s:" << possible_s_coeffs.size() << " d:" << possible_d_coeffs.size() << std::endl;
  vector<vector<vector<double>>> combinations;
  for (int i = 0; i < possible_s_coeffs.size(); i++) {
    for (int j = 0; j < possible_d_coeffs.size(); j++) {
      auto s_coeffs = possible_s_coeffs[i];
      auto d_coeffs = possible_d_coeffs[j];
      combinations.push_back({s_coeffs, d_coeffs});
    }
  }
  return combinations;
}

// check whether the JMT coefficients will break the car limits
// within check_duration
bool check_is_d_JMT_good(vector<double> jmt_coeffs, double check_duration) {
  auto v_coeffs = derivative(jmt_coeffs);
  auto a_coeffs = derivative(v_coeffs);
  for (double t = TIME_STEP; t <= check_duration; t += TIME_STEP) {
    double v = poly_eval(t, v_coeffs);
    double a = abs(poly_eval(t, a_coeffs));
    if (v > SPEED_LIMIT || a > MAX_ACC) {
      return false;
    }
  }
  return true;
}

// check whether the JMT coefficients will break the car limits
// within check_duration
bool check_is_s_JMT_good(vector<double> jmt_coeffs, double check_duration) {
  auto v_coeffs = derivative(jmt_coeffs);
  auto a_coeffs = derivative(v_coeffs);
  auto jerk_coeffs = derivative(a_coeffs);
  for (double t = TIME_STEP; t <= check_duration; t += TIME_STEP) {
    double v = poly_eval(t, v_coeffs);
    double a = abs(poly_eval(t, a_coeffs));
    double jerk = abs(poly_eval(t, jerk_coeffs));
    if (v > SPEED_LIMIT || v < 0.0 || a > MAX_ACC || jerk > MAX_JERK) {
      return false;
    }
  }
  return true;
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


  // setup start config
  double start_s, start_s_d, start_s_dd;
  double start_d, start_d_d, start_d_dd;
  double start_t_offset = 0.0;

  auto prev_s_vals = prev_frenet_trajectory[0];
  auto prev_d_vals = prev_frenet_trajectory[1];

  vector<vector<double>> possible_s_coeffs;
  vector<vector<double>> possible_d_coeffs;

  if (previous_path_x.size() <= 2) {
    std::cout << "generate traj from scratch\n";

    start_s = car["s"];
    start_s_d = car["speed"];
    start_s_dd = 0.0;
    start_d = car["d"];
    start_d_d = 0.0;
    start_d_dd = 0.0;
  } else {

    // remove trajectories that already been comsumed
    while (prev_s_vals.size() > previous_path_x.size()) {
      prev_s_vals.pop_front();
      prev_d_vals.pop_front();
    }

    // TODO(Olala): save the previously end config directly
    assert(prev_s_vals.size() == previous_path_x.size());
    int prev_s_size = prev_s_vals.size();
    double s0 = prev_s_vals[prev_s_size - 3];
    double s1 = prev_s_vals[prev_s_size - 2];
    double s2 = prev_s_vals[prev_s_size - 1];
    start_t_offset = prev_s_size * TIME_STEP;

    if (s2 < s1) s1 -= max_s;
    if (s1 < s0) s0 -= max_s;
    std::cout << "s0,1,2:" << s0 <<", " << s1 << ", " << s2 << std::endl;

    double d0 = prev_d_vals[prev_s_size - 3];
    double d1 = prev_d_vals[prev_s_size - 2];
    double d2 = prev_d_vals[prev_s_size - 1];
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

  // Velocity keeping
  if (current_car_state == velocity_keeping ||
      current_car_state == lane_change_left ||
      current_car_state == lane_change_right) {

    double end_s = start_s + 100.0;
    double min_duration = 0.1;
    double max_duration = 20.0;
    double target_speed = SPEED_LIMIT * 0.9;
    for (double duration = min_duration; duration <= max_duration; duration += 0.25) {
      auto try_end_s = {end_s, target_speed, 0.0};
      auto s_coeffs = JMT(start_s_config, try_end_s, duration);
      if (check_is_s_JMT_good(s_coeffs, duration)){
        possible_s_coeffs.push_back(s_coeffs);
      }
    }

    int target_lane = get_lane(start_d);
    if (current_car_state == lane_change_left ||
        current_car_state == lane_change_right) {
      target_lane = lane_change_target_lane;
    }
    double target_d = 2.0 + 4.0 * target_lane;
    std::cout << "target lane:" << target_lane << ", target_d:" << target_d << std::endl;
    double duration = 4.0;
    double target_d_speed = 0.0;

    for (double d_offset = -0.4; d_offset <= 0.4; d_offset += 0.2) {
      vector<double> try_end_d = {target_d + d_offset, target_d_speed, 0.0};
      auto d_coeffs = JMT(start_d_config, try_end_d, duration);
      if (check_is_d_JMT_good(d_coeffs, duration)) {
        possible_d_coeffs.push_back(d_coeffs);
      }
    }
  } else if (current_car_state == vehicle_following) {

    double car_s = car["s"];
    double car_d = car["d"];
    int leading_vehicle_id = get_leading_vehicle_id(car_s, car_d, sensor_fusion);
    json leading_vehicle = get_leading_vehicle_by_id(leading_vehicle_id, sensor_fusion);

    double lv_x = leading_vehicle[1];
    double lv_y = leading_vehicle[2];
    double lv_vx = leading_vehicle[3];
    double lv_vy = leading_vehicle[4];
    double lv_s = leading_vehicle[5];
    double lv_d = leading_vehicle[6];

    double lv_speed = sqrt(lv_vx * lv_vx + lv_vy * lv_vy);
    double target_v = min(SPEED_LIMIT * 0.9, lv_speed);
    double REACTION_TIME = 0.5;
    double safe_dist_to_lv = target_v * REACTION_TIME + 10.0;

    for (double duration = 1.0; duration <= 10.0; duration += 0.5) {
      for (double target_dist_ahead = safe_dist_to_lv; target_dist_ahead <= safe_dist_to_lv + 5.0; target_dist_ahead += 2.5) {
        double target_s = lv_s + lv_speed * (duration + start_t_offset) - target_dist_ahead;
        // handle target become very small when over max_s
        double dist_ahead = s_dist(start_s, target_s);

        vector<double> try_end_s = {start_s + dist_ahead, target_v, 0.0};
        auto s_coeffs = JMT(start_s_config, try_end_s, duration);
        bool is_traj_good = check_is_s_JMT_good(s_coeffs, duration);
        if (is_traj_good) {
          possible_s_coeffs.push_back(s_coeffs);
        }
      }
    }
    if (possible_s_coeffs.size() == 0) {
      std::cout << "failed to generate trajectory\n";
      std::cout << "start config:\n";
      std::cout << "start_s" << start_s << " start_v:" << start_s_d << std::endl;
    }

    int target_lane = get_lane(start_d);
    double target_d = 2.0 + 4.0 * target_lane;
    double duration = 4.0;
    double target_d_speed = 0.0;

    for (double d_offset = -0.4; d_offset <= 0.4; d_offset += 0.2) {
      vector<double> try_end_d = {target_d + d_offset, target_d_speed, 0.0};
      auto d_coeffs = JMT(start_d_config, try_end_d, duration);
      if (check_is_d_JMT_good(d_coeffs, duration)) {
        possible_d_coeffs.push_back(d_coeffs);
      }
    }
  }

  int total_stpes = PREDICT_HORIZON / TIME_STEP;

  // Calculate future sensor fusion for costs
  // prediction[vehicle_id][i] = {s, d};
  vector<vector<vector<double>>> predictions;
  for (int v = 0; v < sensor_fusion.size(); v++) {

    json vehicle_data = sensor_fusion[v];
    int vehicle_id = vehicle_data[0];
    double vehicle_x = vehicle_data[1];
    double vehicle_y = vehicle_data[2];
    double vehicle_vx = vehicle_data[3];
    double vehicle_vy = vehicle_data[4];
    double vehicle_s = vehicle_data[5];
    double vehicle_d = vehicle_data[6];

    double vehicle_speed = sqrt(vehicle_vx * vehicle_vx + vehicle_vy * vehicle_vy);
    vector<vector<double>> vehicle_prediction;
    for (int i = 0; i < total_stpes; i++) {
      vehicle_s = vehicle_s + vehicle_speed * TIME_STEP;
      vehicle_prediction.push_back({vehicle_s, vehicle_d});
    }
    predictions.push_back(vehicle_prediction);
  }


  // Try to evaluate the cost of different variations
  double min_cost = 1e10;
  vector<deque<double>> best_traj;

  auto combinations = enumerate_coeffs_combs(possible_s_coeffs, possible_d_coeffs);
  assert(combinations.size() > 0);
  for (int c = 0; c < combinations.size(); c++) {
    auto a_comb = combinations[c];
    auto s_coeffs = a_comb[0];
    auto d_coeffs = a_comb[1];


    deque<double> next_s_vals;
    deque<double> next_d_vals;

    // Use some previous trajectories left
    for (int i = 0; i < prev_s_vals.size(); i++) {
      next_s_vals.push_back(prev_s_vals[i]);
      next_d_vals.push_back(prev_d_vals[i]);
    }

    int num_prev_steps = prev_s_vals.size();
    int total_stpes = PREDICT_HORIZON / TIME_STEP;

    for (int i = 0; (i + num_prev_steps) < total_stpes; i++) {
      double t = (i + 1) * TIME_STEP;
      double s = poly_eval(t, s_coeffs);
      double d = poly_eval(t, d_coeffs);

      next_s_vals.push_back(fmod((s + max_s), max_s));
      next_d_vals.push_back(d);
    }

    double cost = 0.0;
    auto trajectory = {next_s_vals, next_d_vals};
    // TODO(Olala): calculate cost
    cost += 1000.0 * collision_cost(trajectory, predictions);
    cost += 50.0 * speed_range_cost(trajectory);
    cost += 20.0 * d_offset_cost(trajectory);

    if (cost < min_cost) {
      min_cost = cost;
      best_traj.clear();
      best_traj.push_back(next_s_vals);
      best_traj.push_back(next_d_vals);
    }
  }


  // std::cout << "best traj:" << best_trajectory_idx << " cost:" << min_cost << std::endl;
  auto best_trajectory = best_traj;

  // save frenet trajectory for next time
  auto next_s_vals = best_trajectory[0];
  auto next_d_vals = best_trajectory[1];
  prev_frenet_trajectory.clear();
  prev_frenet_trajectory.push_back(next_s_vals);
  prev_frenet_trajectory.push_back(next_d_vals);

  // Convert Frenet coordinate to Cartesian coordinate
  vector<double> next_x_vals;
  vector<double> next_y_vals;

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

          // std::cout << "traj_size:" << best_trajectory[0].size() << std::endl;

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

