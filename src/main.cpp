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

          ::Eigen::Matrix<double, 6, 1> state;

          if (ptsx.size() != ptsy.size())
          {
            throw ::std::runtime_error("Waypoint components vectors do not have the same size.");
          }

          const auto initial_cte = ClosestDistanceFromPointToPath(px, py, ptsx, ptsy);
          const auto initial_epsi = AngleErrorAtClosestSegment(px, py, psi, ptsx, ptsy);

          state(0) = px;
          state(1) = py;
          state(2) = psi;
          state(3) = v;
          state(4) = initial_cte;
          state(5) = initial_epsi;

          // Display the MPC predicted trajectory
          vector<double> mpc_x_vals;
          vector<double> mpc_y_vals;

          ::Eigen::VectorXd x_vals_eigen(ptsx.size()), y_vals_eigen(ptsy.size());

          for (int i = 0; i < ptsx.size(); ++i)
          {
            x_vals_eigen(i) = ptsx[i];
            y_vals_eigen(i) = ptsy[i];
          }

          const auto actuator_command = mpc.Solve(state, polyfit(x_vals_eigen, y_vals_eigen, 3), mpc_x_vals, mpc_y_vals);

          steer_value = -1.0 * actuator_command.first / deg2rad(25.0);
          throttle_value = actuator_command.second;

          // Display the waypoints/reference line
          vector<double> next_x_vals;
          vector<double> next_y_vals;

          int closest_index_to_ego = -1;
          double best_dist_to_ego = ::std::numeric_limits<double>::max();

          next_y_vals.clear();

          ::Eigen::MatrixXd reference_points(2, ptsx.size());

          const ::Eigen::Matrix<double, 2, 1> ego_long{::std::cos(psi), ::std::sin(psi)};
          const ::Eigen::Matrix<double, 2, 1> ego_lat{::std::cos(psi + M_PI_2), ::std::sin(psi + M_PI_2)};
          const ::Eigen::Matrix<double, 2, 1> ego_pose{px, py};

          for (int i = 0; i < ptsx.size(); ++i)
          {
            const ::Eigen::Matrix<double, 2, 1> p{ptsx[i], ptsy[i]};
            const ::Eigen::Matrix<double, 2, 1> p_ego_coords{(p - ego_pose).dot(ego_long), (p - ego_pose).dot(ego_lat)};

            next_x_vals.push_back(p_ego_coords(0));
            next_y_vals.push_back(p_ego_coords(1));
          }

          if(mpc_x_vals.size() != mpc_y_vals.size())
          {
            throw ::std::runtime_error("mpc x and y sizes don't match");
          }

          std::vector<double> mpc_x_vals_ego_coords, mpc_y_vals_ego_coords;

          for (int i = 0; i < mpc_x_vals.size(); ++i)
          {
            const ::Eigen::Matrix<double, 2, 1> p{mpc_x_vals[i], mpc_x_vals[i]};
            const ::Eigen::Matrix<double, 2, 1> p_ego_coords{(p - ego_pose).dot(ego_long), (p - ego_pose).dot(ego_lat)};

            mpc_x_vals_ego_coords.push_back(p_ego_coords(0));
            mpc_y_vals_ego_coords.push_back(p_ego_coords(1));
          }

          std::cout << "--------------------------------------------------" << std::endl;

          std::cout << "ego pose (x,y,psi) is: (" << px << ", " << py << ", " << psi << ")" << std::endl;
          std::cout << "steering command is: " << steer_value << std::endl;
          std::cout << "throttle command is: " << throttle_value << std::endl << std::endl;
          std::cout << "waypoints (global frame / ego frame) are: " << std::endl;
          for(int i = 0; i < ptsx.size(); ++i)
          {
            std::cout << "-- (" << ptsx[i] << ", " << ptsy[i] << ")" << " ... ( " << next_x_vals[i] << ", " << next_y_vals[i] << ")" << std::endl;
          }
          std::cout << "mpc best path is: " << std::endl;
          for(int i = 0; i < mpc_x_vals.size(); ++i)
          {
            std::cout << "-- (" << mpc_x_vals[i] << ", " << mpc_y_vals[i] << ")" << " ... ( " << mpc_x_vals_ego_coords[i] << ", " << mpc_y_vals_ego_coords[i] << ")" << std::endl;
          }

          std::cout << std::endl;
          std::cout << "--------------------------------------------------" << std::endl;
          std::cout << std::endl << std::endl << std::endl;

          //convert mpc points into ego vehicle space

          // NOTE: Remember to divide by deg2rad(25) before you send the
          //   steering value back. Otherwise the values will be in between
          //   [-deg2rad(25), deg2rad(25] instead of [-1, 1].
          // multiply steering angle by -1.0 since in the simulator, positive
          // steering angle means 'right turn', while in the MPC update step,
          // positive steering angle means 'left turn'
          msgJson["steering_angle"] = steer_value;
          msgJson["throttle"] = throttle_value;

          msgJson["mpc_x"] = mpc_x_vals_ego_coords;
          msgJson["mpc_y"] = mpc_y_vals_ego_coords;

          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;

          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          // std::cout << msg << std::endl;
          // Latency
          // The purpose is to mimic real driving conditions where
          //   the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          //   around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE SUBMITTING.
          // std::this_thread::sleep_for(std::chrono::milliseconds(100)); ///@todo: be sure to uncomment and handle this before submitting
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