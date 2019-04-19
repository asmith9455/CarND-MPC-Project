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
  MPC mpc{30, 0.1};

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

          // prepare data ----------------------------------------
          // -----------------------------------------------------

          using Vec2D = ::Eigen::Matrix<double, 2, 1>;

          const Vec2D ego_longitudinal_axis__global{::std::cos(psi), ::std::sin(psi)};
          const Vec2D ego_lateral_axis__global{::std::cos(psi + M_PI_2), ::std::sin(psi + M_PI_2)};
          const Vec2D ego_pos__global{px, py};

          const auto ref_line_x__global = ptsx;
          const auto ref_line_y__global = ptsy;

          std::vector<double> ref_line_x__ego, ref_line_y__ego;

          TransformFromGlobalToEgo(
              ego_pos__global, ego_longitudinal_axis__global, ego_lateral_axis__global,
              ref_line_x__global, ref_line_y__global,
              ref_line_x__ego, ref_line_y__ego);

          const auto ref_line_polynomial__ego = polyfit_vecs(ref_line_x__ego, ref_line_y__ego, 3);

          ::Eigen::Matrix<double, 6, 1> state;
          state(0) = 0;
          state(1) = 0;
          state(2) = 0;
          state(3) = v;
          state(4) = ref_line_polynomial__ego(0);
          state(5) = -atan(ref_line_polynomial__ego[1]);

          // -----------------------------------------------------
          // -----------------------------------------------------

          // run mpc----------------------------------------------
          // -----------------------------------------------------
          vector<double> mpc_x_vals__ego;
          vector<double> mpc_y_vals__ego;

          const auto commands = mpc.Solve(state, ref_line_polynomial__ego, mpc_x_vals__ego, mpc_y_vals__ego);

          const double delta = commands.first;
          const double steer_cmd = -1.0 * delta / deg2rad(25.0);
          const double accel = commands.second;
          const double throttle_cmd = accel;

          // -----------------------------------------------------
          // -----------------------------------------------------

          // print debug info-------------------------------------
          // -----------------------------------------------------

          std::cout << "--------------------------------------------------" << std::endl;

          std::cout << "ego pose (x,y,psi) is: (" << px << ", " << py << ", " << psi << ")" << std::endl;
          std::cout << "steering command is: " << steer_cmd << std::endl;
          std::cout << "throttle command is: " << throttle_cmd << std::endl
                    << std::endl;
          std::cout << "waypoints (global frame / ego frame) are: " << std::endl;
          for (int i = 0; i < ref_line_x__global.size(); ++i)
          {
            std::cout << "-- (" << ref_line_x__global[i] << ", " << ref_line_y__global[i] << ")"
                      << " ... ( " << ref_line_x__ego[i] << ", " << ref_line_y__ego[i] << ")" << std::endl;
          }
          std::cout << "mpc best path is: " << std::endl;
          for (int i = 0; i < mpc_x_vals__ego.size(); ++i)
          {
            std::cout << " ( " << mpc_x_vals__ego[i] << ", " << mpc_y_vals__ego[i] << ")" << std::endl;
          }

          std::cout << std::endl;
          std::cout << "--------------------------------------------------" << std::endl;
          std::cout << std::endl
                    << std::endl
                    << std::endl;

          // -----------------------------------------------------
          // -----------------------------------------------------

          //convert mpc points into ego vehicle space

          // NOTE: Remember to divide by deg2rad(25) before you send the
          //   steering value back. Otherwise the values will be in between
          //   [-deg2rad(25), deg2rad(25] instead of [-1, 1].
          // multiply steering angle by -1.0 since in the simulator, positive
          // steering angle means 'right turn', while in the MPC update step,
          // positive steering angle means 'left turn'
          json msgJson;
          msgJson["steering_angle"] = steer_cmd;
          msgJson["throttle"] = throttle_cmd;

          msgJson["mpc_x"] = mpc_x_vals__ego;
          msgJson["mpc_y"] = mpc_y_vals__ego;

          msgJson["next_x"] = ref_line_x__ego;
          msgJson["next_y"] = ref_line_y__ego;

          static auto current = std::chrono::high_resolution_clock::now();
          static auto last = std::chrono::high_resolution_clock::now();

          current = std::chrono::high_resolution_clock::now();

          ::std::cout << "current time: " << std::chrono::system_clock::to_time_t(current) << ::std::endl;
          ::std::cout << "last time: " << std::chrono::system_clock::to_time_t(last) << ::std::endl;
          ::std::cout << "time diff: " << std::chrono::duration<double, std::milli>(current - last).count() << ::std::endl;

          last = current;

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