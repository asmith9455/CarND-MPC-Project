#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <string>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "helpers.h"
#include "json.hpp"
#include "MPC.h"
#include <numeric>
#include <cmath>

// for convenience
using nlohmann::json;
using std::string;
using std::vector;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

int main()
{
  uWS::Hub h;

  // MPC is initialized here!
  MPC mpc{2, 1.0};

  h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
    std::cout << sdata << std::endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2')
    {
      string s = hasData(sdata);
      if (s != "")
      {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry")
        {
          // j[1] is the data JSON object
          vector<double> ptsx = j[1]["ptsx"];
          vector<double> ptsy = j[1]["ptsy"];
          double px = j[1]["x"];
          double py = j[1]["y"];
          double psi = j[1]["psi"];
          double v = j[1]["speed"];

          /**
           * TODO: Calculate steering angle and throttle using MPC.
           * Both are in between [-1, 1].
           */
          double steer_value;
          double throttle_value;

          json msgJson;
          // NOTE: Remember to divide by deg2rad(25) before you send the
          //   steering value back. Otherwise the values will be in between
          //   [-deg2rad(25), deg2rad(25] instead of [-1, 1].
          msgJson["steering_angle"] = steer_value;
          msgJson["throttle"] = throttle_value;

          // Display the MPC predicted trajectory
          vector<double> mpc_x_vals;
          vector<double> mpc_y_vals;

          /**
           * TODO: add (x,y) points to list here, points are in reference to 
           *   the vehicle's coordinate system the points in the simulator are 
           *   connected by a Green line
           */

          mpc_x_vals = {0.0, 10.0};
          mpc_y_vals = {0.0, 0.0};

          msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;

          // Display the waypoints/reference line
          vector<double> next_x_vals;
          vector<double> next_y_vals;

          int closest_index_to_ego = -1;
          double best_dist_to_ego = ::std::numeric_limits<double>::max();

          if (ptsx.size() != ptsy.size())
          {
            throw ::std::runtime_error("Waypoint components vectors do not have the same size.");
          }

          // for (int i = 0; i < ptsx.size(); ++i)
          // {
          //   using ::std::pow;
          //   const auto dist_to_ego = pow(ptsx[i] - px, 2) + pow(ptsy[i] - py, 2);

          //   if (dist_to_ego < best_dist_to_ego)
          //   {
          //     best_dist_to_ego = dist_to_ego;
          //     closest_index_to_ego = i;
          //   }
          // }

          next_x_vals.clear();
          next_y_vals.clear();

          ::Eigen::MatrixXd reference_points(2, ptsx.size());

          const ::Eigen::Matrix<double, 2, 1> ego_long{::std::cos(psi), ::std::sin(psi)};
          const ::Eigen::Matrix<double, 2, 1> ego_lat{::std::cos(psi+M_PI_2), ::std::sin(psi+M_PI_2)};
          const ::Eigen::Matrix<double, 2, 1> ego_pose{px, py};

          // for(int i = closest_index_to_ego; i < 10; ++i)
          // {
          //   const int mod_index{i % ptsx.size()};
          //   next_x_vals.push_back(ptsx.at(mod_index));
          //   next_y_vals.push_back(ptsy.at(mod_index));
          // }

          const auto p1 = ego_long * 0.0 + ego_pose;
          const auto p2 = ego_long * 20.0 + ego_pose;
          const auto p3 = ego_long * 20.0 + ego_pose + ego_lat * 2.0;


          const ::Eigen::Matrix<double,2,1> p1ec{(p1 - ego_pose).dot(ego_long), (p1 - ego_pose).dot(ego_lat)};
          const ::Eigen::Matrix<double,2,1> p2ec{(p2 - ego_pose).dot(ego_long), (p2 - ego_pose).dot(ego_lat)};
          const ::Eigen::Matrix<double,2,1> p3ec{(p3 - ego_pose).dot(ego_long), (p3 - ego_pose).dot(ego_lat)};

          /**
           * TODO: add (x,y) points to list here, points are in reference to 
           *   the vehicle's coordinate system the points in the simulator are 
           *   connected by a Yellow line
           */

          for(int i = 0; i < ptsx.size(); ++i)
          {
            const ::Eigen::Matrix<double,2,1> p{ptsx[i], ptsy[i]};
            const ::Eigen::Matrix<double,2,1> p_ego_coords{(p - ego_pose).dot(ego_long), (p - ego_pose).dot(ego_lat)};

            next_x_vals.push_back(p_ego_coords(0));
            next_y_vals.push_back(p_ego_coords(1));
          }

          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;

          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          std::cout << msg << std::endl;
          // Latency
          // The purpose is to mimic real driving conditions where
          //   the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          //   around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE SUBMITTING.
          std::this_thread::sleep_for(std::chrono::milliseconds(100));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        } // end "telemetry" if
      }
      else
      {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    } // end websocket if
  }); // end h.onMessage

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port))
  {
    std::cout << "Listening to port " << port << std::endl;
  }
  else
  {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }

  h.run();
}