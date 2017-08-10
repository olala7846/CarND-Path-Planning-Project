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


// save previous frenet trajectory globally
vector<deque<double>> prev_frenet_trajectory;

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

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s,
    vector<double> maps_x, vector<double> maps_y, vector<double> maps_dx, vector<double> maps_dy)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

  // calculate XY with spline
  int spline_points = 6;
  int map_size = maps_x.size();
  vector<double> X(spline_points), Y(spline_points);
  for(int i=0; i < spline_points; i++) {
    int wp = (prev_wp + i + map_size - 2) % map_size;
    X[i] = maps_x[wp];
    Y[i] = maps_y[wp];
  }

  // affine transformation
  double x0 = X[0];
  double y0 = Y[0];
  double theta = atan2(Y[1] - Y[0], X[1] - X[0]);

  vector<double> _X(spline_points), _Y(spline_points);
  for(int i = 0; i < spline_points; i++) {
    // translate P0 to origin
    double x_t = X[i] - x0;
    double y_t = Y[i] - y0;
    _X[i] = x_t * cos(-theta) - y_t * sin(-theta);
    _Y[i] = x_t * sin(-theta) + y_t * cos(-theta);
  }

  // fit spline
  tk::spline spline_func;
  spline_func.set_points(_X, _Y);
  double ratio = (s - maps_s[prev_wp]) / (maps_s[wp2]-maps_s[prev_wp]);
  double _x = _X[2] + (_X[3] - _X[2]) * ratio;
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


/* (Olala): PTG: polynomial trajectory generating */
string get_behavior(string current_behavior, const json &sensor_fusion) {
  string behavior ("full_speed_ahead");
  return behavior;
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

vector<vector<double>> generate_trajectories(
    double target_s, double target_speed, int target_lane,
    double time_to_goal, const json &car, const json &sensor_fusion,
    vector<double> maps_s, vector<double> maps_x, vector<double> maps_y,
    vector<double> maps_dx, vector<double> maps_dy) {

  const double time_step = 0.02;
  const double time_to_predict = 2.5;  // only predict 2 seconds ahead
  const double time_use_predicted = 1.0;  // keep 1 second of predicted path

  vector<vector<vector<double>>> trajectories;
  vector<vector<deque<double>>> frenet_trajectories;

  double car_x = car["x"];
  double car_y = car["y"];
  double car_s = car["s"];
  double car_d = car["d"];
  double car_yaw = car["yaw"];
  double car_speed = car["speed"];

  car_speed = car_speed * 0.44704;  // 1 mph = 0.44704 m/s

  auto previous_path_x = car["previous_path_x"];
  auto previous_path_y = car["previous_path_y"];

  // Previous path's end s and d values
  double end_path_s = car["end_path_s"];
  double end_path_d = car["end_path_d"];

  double target_d = target_lane * 4.0 + 2.0;

  deque<double> next_s_vals;
  deque<double> next_d_vals;

  vector<double> next_x_vals;
  vector<double> next_y_vals;

  for(double dt=-1.0; dt <= 1.0; dt += 0.1) {

    if (previous_path_x.size() <= 3){  // generate trajectory from scratch

      vector<double> start_s = {car_s, car_speed, 0.0};
      vector<double> end_s = {target_s, target_speed, 0.0};
      auto s_trajectory_coeff = JMT(start_s, end_s, time_to_goal);

      vector<double> start_d = {car_d, 0.0, 0.0};
      vector<double> end_d = {target_d, 0.0, 0.0};
      auto d_trajectory_coeff = JMT(start_d, end_d, time_to_goal);

      for(double t=time_step; t<= time_to_predict; t+= time_step){
        double s = poly_eval(t, s_trajectory_coeff);
        double d = poly_eval(t, d_trajectory_coeff);
        next_s_vals.push_back(s);
        next_d_vals.push_back(d);
        vector<double> point = getXY(s, d, maps_s, maps_x, maps_y, maps_dx, maps_dy);
        next_x_vals.push_back(point[0]);
        next_y_vals.push_back(point[1]);
      }


    } else {  // use some previous generated trajectory

      auto prev_s_vals = prev_frenet_trajectory[0];
      auto prev_d_vals = prev_frenet_trajectory[1];

      // remove trajectories that already been comsumed
      for(int i=0; i< previous_path_x.size() - prev_s_vals.size(); i++) {
          prev_s_vals.pop_front();
          prev_d_vals.pop_front();
      }

      int steps_to_keep = (time_use_predicted + 0.01) / time_step;

      for (int i=0; i<steps_to_keep; i++) {
        next_s_vals.push_back(prev_s_vals[i]);
        next_d_vals.push_back(prev_d_vals[i]);
        next_x_vals.push_back(previous_path_x[i]);
        next_y_vals.push_back(previous_path_y[i]);
      }

      // start JMT from the last config
      double s5 = prev_s_vals[steps_to_keep-1];
      double s4 = prev_s_vals[steps_to_keep-2];
      double s3 = prev_s_vals[steps_to_keep-3];
      double d5 = prev_d_vals[steps_to_keep-1];
      double d4 = prev_d_vals[steps_to_keep-2];
      double d3 = prev_d_vals[steps_to_keep-3];

      double vs5 = (s5 - s4) / time_step;
      double vs4 = (s4 - s3) / time_step;

      double vd5 = (d5 - d4) / time_step;
      double vd4 = (d4 - d3) / time_step;

      double as5 = (vs5 - vs4) / time_step;
      double ad5 = (vd5 - vd4) / time_step;

      vector<double> start_s = {s5, vs5, as5};
      vector<double> end_s = {target_s, target_speed, 0.0};
      auto s_trajectory_coeff = JMT(start_s, end_s, time_to_goal);

      vector<double> start_d = {d5, vd5, ad5};
      vector<double> end_d = {target_d, 0.0, 0.0};
      auto d_trajectory_coeff = JMT(start_d, end_d, time_to_goal);

      for(double t=time_step; t <= ( time_to_predict - time_use_predicted) ; t += time_step){
        double s = poly_eval(t, s_trajectory_coeff);
        double d = poly_eval(t, d_trajectory_coeff);
        next_s_vals.push_back(s);
        next_d_vals.push_back(d);
        vector<double> point = getXY(s, d, maps_s, maps_x, maps_y, maps_dx, maps_dy);
        next_x_vals.push_back(point[0]);
        next_y_vals.push_back(point[1]);
      }

    }

    vector<deque<double>> frenet_trajectory = {next_s_vals, next_d_vals};
    frenet_trajectories.push_back(frenet_trajectory);

    vector<vector<double>> trajectory = {next_x_vals, next_y_vals};
    trajectories.push_back(trajectory);

  }

  auto best_trajectory = trajectories[0];
  auto best_frenet_trajectory = frenet_trajectories[0];

  prev_frenet_trajectory.clear();
  prev_frenet_trajectory.push_back(best_frenet_trajectory[0]);
  prev_frenet_trajectory.push_back(best_frenet_trajectory[1]);

  return trajectories[0];
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

  string current_behavior ("full_speed_ahead");

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

  h.onMessage([
      &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy,
      &current_behavior](
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
          // first skip this step
          string new_behavior = get_behavior(current_behavior, sensor_fusion);
          current_behavior = new_behavior;

          double target_speed = 49.0 * 0.44704; // to meters per second
          double T = 4.0;
          double target_s = car_s + T * (target_speed - 0.1);
          int target_lane = 1;

          // Trajectory Genereation
          auto best_trajectory = generate_trajectories(
              target_s, target_speed, target_lane, T, j[1], sensor_fusion,
              map_waypoints_s, map_waypoints_x, map_waypoints_y, map_waypoints_dx, map_waypoints_dy);
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

