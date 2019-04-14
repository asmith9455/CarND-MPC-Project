#ifndef HELPERS_H
#define HELPERS_H

#include <string>
#include "Eigen-3.3/Eigen/Core"

using Eigen::VectorXd;
using std::string;

double clamp(const double v, const double min, const double max)
{
  return std::min(std::max(v, min), max);
}

double ClosestDistanceFromPointToPath(
    const double px,
    const double py,
    const ::std::vector<double> global_path_x_values,
    const ::std::vector<double> global_path_y_values)
{
  if (global_path_x_values.size() != global_path_y_values.size())
  {
    throw ::std::runtime_error("Path component vector sizes are not equivalent.");
  }

  double best_mag = ::std::numeric_limits<double>::max();

  for (int i = 0; i < global_path_x_values.size() - 1; ++i)
  {
    using Vec2D = ::Eigen::Matrix<double, 2, 1>;

    Vec2D v{px, py};
    Vec2D p1{global_path_x_values[i], global_path_y_values[i]};
    Vec2D p2{global_path_x_values[i+1], global_path_y_values[i+1]};

    const auto seg = p2 - p1;
    const auto join = p2 - v;

    if (seg.norm() < 1e-10)
    {
      throw ::std::runtime_error("~0 length segment");
    }

    const auto seg_unit = seg / seg.norm();

    const auto join_along_seg = clamp(join.dot(seg_unit), 0.0, seg.norm()) * seg_unit;

    const auto point_on_seg = p2 + join_along_seg;

    const auto closest_dist = (point_on_seg - v).norm();

    if (closest_dist < best_mag)
    {
      best_mag = closest_dist;
    }
  }

  return best_mag;
}

double AngleErrorAtClosestSegment(
    const double px,
    const double py,
    const double psi,
    const ::std::vector<double> global_path_x_values,
    const ::std::vector<double> global_path_y_values)
{
  if (global_path_x_values.size() != global_path_y_values.size())
  {
    throw ::std::runtime_error("Path component vector sizes are not equivalent.");
  }

  using Vec2D = ::Eigen::Matrix<double, 2, 1>;

  double best_mag = ::std::numeric_limits<double>::max();
  Vec2D best_seg{};

  for (int i = 0; i < global_path_x_values.size() - 1; ++i)
  {
    Vec2D v{px, py};
    Vec2D p1{global_path_x_values[i], global_path_y_values[i]};
    Vec2D p2{global_path_x_values[i+1], global_path_y_values[i+1]};

    const auto seg = p2 - p1;
    const auto join = p2 - v;

    if (seg.norm() < 1e-10)
    {
      throw ::std::runtime_error("~0 length segment");
    }
    
    const auto seg_unit = seg / seg.norm();

    const auto join_along_seg = clamp(join.dot(seg_unit), 0.0, seg.norm()) * seg_unit;

    const auto point_on_seg = p2 + join_along_seg;

    const auto closest_dist = (point_on_seg - v).norm();

    if (closest_dist < best_mag)
    {
      best_mag = closest_dist;
      best_seg = seg;
    }
  }

  Vec2D vehicle_orientation{::std::cos(psi), ::std::sin(psi)};

  //compute angle
  const auto theta = ::std::acos(best_seg.dot(vehicle_orientation) / best_seg.norm());

  return theta;
}

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s)
{
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos)
  {
    return "";
  }
  else if (b1 != string::npos && b2 != string::npos)
  {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

//
// Helper functions to fit and evaluate polynomials.
//

// Evaluate a polynomial.
double polyeval(const VectorXd &coeffs, double x)
{
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); ++i)
  {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Fit a polynomial.
// Adapted from:
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
VectorXd polyfit(const VectorXd &xvals, const VectorXd &yvals, int order)
{
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);

  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); ++i)
  {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); ++j)
  {
    for (int i = 0; i < order; ++i)
    {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);

  return result;
}

#endif // HELPERS_H