#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include <iostream>
#include <string>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include <numeric>
#include <deque>
// #include <matplotlibcpp.h>

// namespace plt = matplotlibcpp;

using CppAD::AD;
using Eigen::VectorXd;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
//   simulator around in a circle with a constant steering angle and velocity on
//   a flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
//   presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;
const double reference_accel_rate = 1.0; //accelerate at 1 m/s/s to reach desired speed
const double reference_speed = 10.0;     //desired speed

struct TrajectoryCosts
{
  using CostT = AD<double>;

  TrajectoryCosts() : cte{0}, epsi{0}, velocity{0},
                      delta_magnitude{0}, accel_magnitude{0},
                      delta_smoothness{0}, accel_smoothness{0}
  {
  }

  CostT TotalCost() const
  {
    return cte + epsi + velocity + delta_magnitude + accel_magnitude + delta_smoothness + accel_smoothness;
  }

  CostT cte;
  CostT epsi;
  CostT velocity;
  CostT delta_magnitude;
  CostT accel_magnitude;
  CostT delta_smoothness;
  CostT accel_smoothness;
};

using ADvector = CPPAD_TESTVECTOR(AD<double>);

template <typename ContainerT>
TrajectoryCosts CalculateTrajectoryCosts(const ContainerT &vars, const ::std::size_t N, const std::double_t dt)
{
  TrajectoryCosts costs;

  const ::std::size_t x_start{0};
  const ::std::size_t y_start{N};
  const ::std::size_t psi_start{y_start + N};
  const ::std::size_t v_start{psi_start + N};
  const ::std::size_t cte_start{v_start + N};
  const ::std::size_t epsi_start{cte_start + N};
  const ::std::size_t delta_start{epsi_start + N};
  const ::std::size_t a_start{delta_start + (N - 1)};
  const ::std::double_t ref_v{10.0};

  for (int k = 1; k < N; ++k)
  {
    const double t = dt * k;

    //reconstruct current speed
    auto predicted_speed_at_t = vars[v_start];
    for (int i = 1; i <= k; ++i)
    {
      predicted_speed_at_t += vars[a_start + i - 1] * dt;
    }
    const auto current_speed = vars[v_start];

    const auto current_speed_error = (reference_speed - current_speed);

    const auto expected_time_to_achieve_ref_speed = current_speed_error / reference_accel_rate;
    const auto time_ratio = t / expected_time_to_achieve_ref_speed;
    const auto ref_speed_for_this_time = time_ratio * current_speed_error + current_speed;

    const auto d_v_cost = ::CppAD::pow(predicted_speed_at_t - reference_speed, 2);
    costs.velocity += d_v_cost;
  }

  for (int k = 0; k < N; ++k)
  {
    const auto d_cte_cost = ::CppAD::pow(vars[cte_start + k], 2);
    const auto d_epsi_cost = ::CppAD::pow(vars[epsi_start + k], 2);

    costs.cte += d_cte_cost;
    costs.epsi += d_epsi_cost;
  }

  for (int k = 0; k < (N - 1); ++k)
  {
    const auto d_delta_mag_cost = ::CppAD::pow(vars[delta_start + k], 2);
    const auto d_accel_mag_cost = ::CppAD::pow(vars[a_start + k], 2);
    costs.delta_magnitude += d_delta_mag_cost;
    costs.accel_magnitude += d_delta_mag_cost;
  }

  for (int k = 0; k < (N - 2); ++k)
  {
    const auto d_delta_smoothness = ::CppAD::pow(vars[delta_start + k + 1] - vars[delta_start + k], 2);
    const auto d_accel_smoothness = 1e4 * ::CppAD::pow(vars[a_start + k + 1] - vars[a_start + k], 2);

    costs.delta_smoothness += d_delta_smoothness;
    costs.accel_smoothness += d_accel_smoothness;
  }

  return costs;
}

void PrintCostDebugInfo(const TrajectoryCosts &costs)
{
  ::std::cout << "------------------------------------------------" << ::std::endl;
  ::std::cout << "Trajectory Costs Are: " << ::std::endl;
  ::std::cout << "-- cte: " << costs.cte << ::std::endl;
  ::std::cout << "-- epsi: " << costs.epsi << ::std::endl;
  ::std::cout << "-- velocity: " << costs.velocity << ::std::endl;
  ::std::cout << "-- delta_magnitude: " << costs.delta_magnitude << ::std::endl;
  ::std::cout << "-- accel_magnitude: " << costs.accel_magnitude << ::std::endl;
  ::std::cout << "-- delta_smoothness: " << costs.delta_smoothness << ::std::endl;
  ::std::cout << "-- accel_smoothness: " << costs.accel_smoothness << ::std::endl;
  ::std::cout << "------------------------------------------------" << ::std::endl;
}

class FG_eval
{
public:
  using ADvector = CPPAD_TESTVECTOR(AD<double>);
  // Fitted polynomial coefficients
  VectorXd coeffs;
  FG_eval(VectorXd coeffs, ::std::size_t N, ::std::double_t dt) : N_{N}, dt_{dt} { this->coeffs = coeffs; }

  void operator()(ADvector &fg, const ADvector &vars)
  {
    fg[0] = 0;

    const ::std::size_t x_start{0};
    const ::std::size_t y_start{N_};
    const ::std::size_t psi_start{y_start + N_};
    const ::std::size_t v_start{psi_start + N_};
    const ::std::size_t cte_start{v_start + N_};
    const ::std::size_t epsi_start{cte_start + N_};
    const ::std::size_t delta_start{epsi_start + N_};
    const ::std::size_t a_start{delta_start + (N_ - 1)};
    const ::std::double_t ref_v{10.0};

    AD<double> velocity_cost{0};
    AD<double> cte_cost{0};
    AD<double> epsi_cost{0};
    AD<double> smooth_input_cost{0};

    fg[1 + x_start] = vars[x_start];
    fg[1 + y_start] = vars[y_start];
    fg[1 + psi_start] = vars[psi_start];
    fg[1 + v_start] = vars[v_start];
    fg[1 + cte_start] = vars[cte_start];
    fg[1 + epsi_start] = vars[epsi_start];

    for (int k = 1; k < N_; ++k)
    {
      /// as in quiz, calculate predictions based on optimizer inputs accel and delta
      /// adapted from the code presented in the MPC lesson here: https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/f1820894-8322-4bb3-81aa-b26b3c6dcbaf/lessons/338b458f-7ebf-449c-9ad1-611eb933b076/concepts/ee21948d-7fad-4821-b61c-0d69bbfcc425
      const AD<double> x1 = vars[x_start + k];
      const AD<double> y1 = vars[y_start + k];
      const AD<double> psi1 = vars[psi_start + k];
      const AD<double> v1 = vars[v_start + k];
      const AD<double> cte1 = vars[cte_start + k];
      const AD<double> epsi1 = vars[epsi_start + k];

      const AD<double> x0 = vars[x_start + k - 1];
      const AD<double> y0 = vars[y_start + k - 1];
      const AD<double> psi0 = vars[psi_start + k - 1];
      const AD<double> v0 = vars[v_start + k - 1];
      const AD<double> cte0 = vars[cte_start + k - 1];
      const AD<double> epsi0 = vars[epsi_start + k - 1];

      const AD<double> delta0 = vars[delta_start + k - 1];
      const AD<double> a0 = vars[a_start + k - 1];

      const AD<double> f0 = coeffs[0] + coeffs[1] * x0;
      const AD<double> psides0 = CppAD::atan(coeffs[1]);

      fg[1 + x_start + k] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt_);
      fg[1 + y_start + k] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt_);
      fg[1 + psi_start + k] = psi1 - (psi0 + v0 * delta0 / Lf * dt_);
      fg[1 + v_start + k] = v1 - (v0 + a0 * dt_);
      fg[1 + cte_start + k] =
          cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt_));
      fg[1 + epsi_start + k] =
          epsi1 - ((psi0 - psides0) + v0 * delta0 / Lf * dt_);
    }

    const auto costs = CalculateTrajectoryCosts(vars, N_, dt_);

    // PrintCostDebugInfo(costs);

    fg[0] = costs.TotalCost();
  }

private:
  ::std::size_t N_;
  ::std::double_t dt_;
};

//
// MPC class definition implementation.
//
MPC::MPC(::std::size_t N, ::std::double_t dt) : N_{N}, dt_{dt} {
  // plt::figure_size(1200, 780);
}
MPC::~MPC() {}

::std::pair<double, double> MPC::Solve(const VectorXd &state, const VectorXd &coeffs, ::std::vector<double> &x_vals, ::std::vector<double> &y_vals)
{
  bool ok = true;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  /**
   * TODO: Set the number of model variables (includes both states and inputs).
   * For example: If the state is a 4 element vector, the actuators is a 2
   *   element vector and there are 10 timesteps. The number of variables is:
   *   4 * 10 + 2 * 9
   */
  const ::std::size_t state_size{6};
  const ::std::size_t actuator_count{2};

  ::std::size_t n_vars = state_size * N_ + actuator_count * (N_ - 1);
  /**
   * TODO: Set the number of constraints
   */
  size_t n_constraints = state_size * N_;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; ++i)
  {
    vars[i] = 0;
  }

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);

  const ::std::size_t x_start{0};
  const ::std::size_t y_start{N_};
  const ::std::size_t psi_start{y_start + N_};
  const ::std::size_t v_start{psi_start + N_};
  const ::std::size_t cte_start{v_start + N_};
  const ::std::size_t epsi_start{cte_start + N_};
  const ::std::size_t delta_start{epsi_start + N_};
  const ::std::size_t a_start{delta_start + (N_ - 1)};

  for (int i = 0; i < delta_start; ++i)
  {
    vars_lowerbound[i] = -std::numeric_limits<double>::max();
    vars_upperbound[i] = std::numeric_limits<double>::max();
  }

  for (int i = 0; i < N_ - 1; ++i)
  {
    vars_lowerbound[delta_start + i] = Lf * 0.436332;
    vars_upperbound[delta_start + i] = Lf * 0.436332;
  }

  vars_lowerbound[x_start] = state(0);
  vars_upperbound[x_start] = state(0);

  vars_lowerbound[y_start] = state(1);
  vars_upperbound[y_start] = state(1);

  vars_lowerbound[psi_start] = state(2);
  vars_upperbound[psi_start] = state(2);

  vars_lowerbound[v_start] = state(3);
  vars_upperbound[v_start] = state(3);

  vars_lowerbound[cte_start] = state(4);
  vars_upperbound[cte_start] = state(4);

  vars_lowerbound[epsi_start] = state(5);
  vars_upperbound[epsi_start] = state(5);

  for (int i = 0; i < N_ - 1; ++i)
  {
    vars_lowerbound[a_start + i] = -1.0;
    vars_upperbound[a_start + i] = 1.0;
  }
  /**
   * TODO: Set lower and upper limits for variables.
   */

  vars[x_start] = state(0);
  vars[y_start] = state(1);
  vars[psi_start] = state(2);
  vars[v_start] = state(3);
  vars[cte_start] = state(4);
  vars[epsi_start] = state(5);

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; ++i)
  {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs, N_, dt_);

  // NOTE: You don't have to worry about these options
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  //   of sparse routines, this makes the computation MUCH FASTER. If you can
  //   uncomment 1 of these and see if it makes a difference or not but if you
  //   uncomment both the computation time should go up in orders of magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          5.0\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "BEST SOLUTION PROPERTIES: " << std::endl;
  std::cout << "--success : " << (solution.status == CppAD::ipopt::solve_result<Dvector>::success) << std::endl;
  std::cout << "--cost : " << cost << std::endl;
  std::cout << "--(delta, accel) : (" << solution.x[delta_start] << ", " << solution.x[a_start] << ")" << std::endl;

  /**
   * TODO: Return the first actuator values. The variables can be accessed with
   *   `solution.x[i]`.
   *
   * {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
   *   creates a 2 element double vector.
   */

  for (int i = 0; i < N_; ++i)
  {
    x_vals.push_back(solution.x[x_start + i]);
  }
  for (int i = 0; i < N_; ++i)
  {
    y_vals.push_back(solution.x[y_start + i]);
  }

  assert(x_vals.size() == y_vals.size());

  const auto costs = CalculateTrajectoryCosts(solution.x, N_, dt_);

  PrintCostDebugInfo(costs);
  
  static std::deque<::std::double_t> vel_costs;

  vel_costs.push_back(CppAD::Value(costs.velocity));

  if (vel_costs.size() > 30)
  {
    vel_costs.pop_front();
  }

  std::cout << "velocity costs latest to oldest:" << std::endl;
  for (auto it = vel_costs.rbegin(); it != vel_costs.rend(); ++it)
  {
    std::cout << *it << ", ";
  }

  // plt::plot(vel_costs);

  return {solution.x[delta_start], solution.x[a_start]};
}