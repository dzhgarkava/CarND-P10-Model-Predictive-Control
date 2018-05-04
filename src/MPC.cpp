#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// TODO: Set the timestep length and duration
size_t N = 10;
double dt = 0.1;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

// Order:
// x1 - xN
// y1 - yN
// psi1 - psiN
// v1 - vN
// cte1 - cteN
// epsi1 - epsiN
// delta1 - delta(N-1)
// a1 - a(N-1)
int x_index = 0;
int y_index = x_index + N;
int psi_index = y_index + N;
int v_index = psi_index + N;
int cte_index = v_index + N;
int epsi_index = cte_index + N;
int delta_index = epsi_index + N;
int a_index = delta_index + N - 1;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.
    fg[0] = 0;

    // The part of the cost based on the reference state.
    for (int t = 0; t < N; t++)
    {
      fg[0] += CppAD::pow(vars[cte_index + t] , 2);
      fg[0] += CppAD::pow(vars[epsi_index + t], 2);
      fg[0] += CppAD::pow(vars[v_index + t], 2);
    }

    // Minimize change-rate.
    for (int t = 0; t < N - 1; t++)
    {
      fg[0] += CppAD::pow(vars[delta_index + t], 2);
      fg[0] += CppAD::pow(vars[a_index + t], 2);
    }

    // Minimize the value gap between sequential actuations.
    for (int t = 0; t < N - 2; t++) {
      fg[0] += CppAD::pow(vars[delta_index + t + 1] - vars[delta_index + t], 2);
      fg[0] += CppAD::pow(vars[a_index + t + 1] - vars[a_index + t], 2);
    }

    fg[1 + x_index] = vars[x_index];
    fg[1 + y_index] = vars[y_index];
    fg[1 + psi_index] = vars[psi_index];
    fg[1 + v_index] = vars[v_index];
    fg[1 + cte_index] = vars[cte_index];
    fg[1 + epsi_index] = vars[epsi_index];

    for (int t = 1; t < N ; t++)
    {
      AD<double> x1 = vars[x_index + t];
      AD<double> y1 = vars[y_index + t];
      AD<double> psi1 = vars[psi_index + t];
      AD<double> v1 = vars[v_index + t];
      AD<double> cte1 = vars[cte_index + t];
      AD<double> epsi1 = vars[epsi_index + t];

      // The state at time t.
      AD<double> x0 = vars[x_index + t - 1];
      AD<double> y0 = vars[y_index + t - 1];
      AD<double> psi0 = vars[psi_index + t - 1];
      AD<double> v0 = vars[v_index + t - 1];
      AD<double> cte0 = vars[cte_index + t - 1];
      AD<double> epsi0 = vars[epsi_index + t - 1];

      // Only consider the actuation at time t.
      AD<double> delta0 = vars[delta_index + t - 1];
      AD<double> a0 = vars[a_index + t - 1];

      AD<double> f0 = coeffs[0] + coeffs[1] * x0;
      AD<double> psides0 = CppAD::atan(coeffs[1]);


      // Here's `x` to get you started.
      // The idea here is to constraint this value to be 0.
      //
      // Recall the equations for the model:
      // x_[t] = x[t-1] + v[t-1] * cos(psi[t-1]) * dt
      // y_[t] = y[t-1] + v[t-1] * sin(psi[t-1]) * dt
      // psi_[t] = psi[t-1] + v[t-1] / Lf * delta[t-1] * dt
      // v_[t] = v[t-1] + a[t-1] * dt
      // cte[t] = f(x[t-1]) - y[t-1] + v[t-1] * sin(epsi[t-1]) * dt
      // epsi[t] = psi[t] - psides[t-1] + v[t-1] * delta[t-1] / Lf * dt
      fg[1 + x_index + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[1 + y_index + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[1 + psi_index + t] = psi1 - (psi0 + v0 * delta0 / Lf * dt);
      fg[1 + v_index + t] = v1 - (v0 + a0 * dt);
      fg[1 + cte_index + t] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
      fg[1 + epsi_index + t] = epsi1 - ((psi0 - psides0) + v0 * delta0 / Lf * dt);
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  size_t n_vars = 6 * N + 2 * (N -1 );
  // TODO: Set the number of constraints
  size_t n_constraints = 6 * N;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++)
  {
    vars[i] = 0;
  }

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);

    // TODO: Set lower and upper limits for variables.

  // Set all non-actuators upper and lowerlimits
  // to the max negative and positive values.
  for (int i = 0; i < delta_index; i++)
  {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  for (auto i = delta_index; i < a_index; i++)
  {
    vars_lowerbound[i] = -25 * M_PI / 180;
    vars_upperbound[i] = 25 * M_PI / 180;
  }

  // Acceleration upper and lower limits.
  for (auto i = a_index; i < n_vars; i++)
  {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++)
  {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

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
  std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
    vector<double> ret;

    ret.push_back(solution.x[delta_index]);
    ret.push_back(solution.x[a_index]);

    for (int i = 0; i < N - 1; ++i)
    {
        ret.push_back(solution.x[x_index + i]);
        ret.push_back(solution.x[y_index + i]);
    }

    return ret;
}
