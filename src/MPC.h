#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"
// #include <matplotlibcpp.h>

class MPC
{
public:
  MPC(::std::size_t N, ::std::double_t dt);

  virtual ~MPC();

  using VectorXd = ::Eigen::VectorXd;

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuations.
  ::std::pair<double, double> Solve(const VectorXd &state, const VectorXd &coeffs, ::std::vector<double> &x_vals, ::std::vector<double> &y_vals);

private:
  ::std::size_t N_;
  ::std::double_t dt_;
  ::std::size_t plot_history_size_;
  ::std::vector<double> x_vals_;
  // ::matplotlibcpp::Plot delta_plot_;
  // ::matplotlibcpp::Plot accel_plot_;
};

#endif // MPC_H
