//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NewtonIteration.h"

#include "InputParameters.h"
#include "Conversion.h"
#include "MooseEnum.h"
#include "Moose.h"
#include "MathUtils.h"

#include <limits>
#include <string>
#include <cmath>

// template <>
// InputParameters
// validParams<NewtonIteration>()
// {
//   InputParameters params = emptyInputParameters();

InputParameters
NewtonIteration::validParams()
{
  InputParameters params = emptyInputParameters();

  // Newton iteration control parameters
  params.addParam<Real>(
      "relative_tolerance", 1e-8, "Relative convergence tolerance for Newton iteration");
  params.addParam<Real>(
      "absolute_tolerance", 1e-14, "Absolute convergence tolerance for Newton iteration");
  params.addParam<Real>("acceptable_multiplier",
                        10,
                        "Factor applied to relative and absolute "
                        "tolerance for acceptable convergence if "
                        "iterations are no longer making progress");

  // diagnostic output parameters
  MooseEnum internal_solve_output_on_enum("never on_error always", "on_error");
  params.addParam<MooseEnum>("internal_solve_output_on",
                             internal_solve_output_on_enum,
                             "When to output internal Newton solve information");
  params.addParam<bool>("internal_solve_full_iteration_history",
                        false,
                        "Set true to output full internal Newton iteration history at times "
                        "determined by `internal_solve_output_on`. If false, only a summary is "
                        "output.");
  params.addParamNamesToGroup("internal_solve_output_on internal_solve_full_iteration_history",
                              "Debug");

  return params;
}

NewtonIteration::NewtonIteration(
    const InputParameters & parameters)
  : _check_range(true),
    _line_search(true),
    _bracket_solution(true),
    _internal_solve_output_on(
        parameters.get<MooseEnum>("internal_solve_output_on").getEnum<InternalSolveOutput>()),
    _max_its(1000), // Far larger than ever expected to be needed
    _internal_solve_full_iteration_history(
        parameters.get<bool>("internal_solve_full_iteration_history")),
    _internal_iteration_rel_tol(parameters.get<Real>("relative_tolerance")),
    _internal_iteration_abs_tol(parameters.get<Real>("absolute_tolerance")),
    _acceptable_multiplier(parameters.get<Real>("acceptable_multiplier")),
    _num_resids(30),
    _residual_history(_num_resids, std::numeric_limits<Real>::max()),
    _iteration(0),
    _initial_residual(0.0),
    _residual(0.0),
    _svrms_name(parameters.get<std::string>("_object_name"))
{
}

Real
NewtonIteration::minimumPermissibleValue(
    const Real /*trial_value*/) const
{
  return std::numeric_limits<Real>::lowest();
}

Real
NewtonIteration::maximumPermissibleValue(
    const Real /*trail_value*/) const
{
  return std::numeric_limits<Real>::max();
}

void
NewtonIteration::returnNewtonSolve(const Real trial_value,
                                                        Real & scalar,
                                                        const ConsoleStream & console)
{
  // construct the stringstream here only if the debug level is set to ALL
  std::stringstream * iter_output =
      (_internal_solve_output_on == InternalSolveOutput::ALWAYS) ? new std::stringstream : nullptr;

  // do the internal solve and capture iteration info during the first round
  // iff full history output is requested regardless of whether the solve failed or succeeded
  auto solve_state = internalSolve(trial_value,
                                   scalar,
                                   _internal_solve_full_iteration_history ? iter_output : nullptr);
  if (solve_state != SolveState::SUCCESS &&
      _internal_solve_output_on != InternalSolveOutput::ALWAYS)
  {
    // output suppressed by user, throw immediately
    if (_internal_solve_output_on == InternalSolveOutput::NEVER)
      throw MooseException("");

    // user expects some kind of output, if necessary setup output stream now
    if (!iter_output)
      iter_output = new std::stringstream;

    // add the appropriate error message to the output
    switch (solve_state)
    {
      case SolveState::NAN_INF:
        *iter_output << "Encountered inf or nan in material return mapping iterations.\n";
        break;

      case SolveState::EXCEEDED_ITERATIONS:
        *iter_output << "Exceeded maximum iterations in material return mapping iterations.\n";
        break;

      default:
        mooseError("Unhandled solver state");
    }

    // if full history output is only requested for failed solves we have to repeat
    // the solve a second time
    if (_internal_solve_full_iteration_history)
      internalSolve(trial_value, scalar, iter_output);

    // Append summary and throw exception
    outputIterationSummary(iter_output, _iteration);
    throw MooseException(iter_output->str());
  }

  if (_internal_solve_output_on == InternalSolveOutput::ALWAYS)
  {
    // the solve did not fail but the user requested debug output anyways
    outputIterationSummary(iter_output, _iteration);
    console << iter_output->str();
  }
}

NewtonIteration::SolveState
NewtonIteration::internalSolve(const Real trial_value,
                                                   Real & scalar,
                                                   std::stringstream * iter_output)
{
  scalar = initialGuess(trial_value);
  Real scalar_old = scalar;
  Real scalar_increment = 0.0;
  const Real min_permissible_scalar = minimumPermissibleValue(trial_value);
  const Real max_permissible_scalar = maximumPermissibleValue(trial_value);
  Real scalar_upper_bound = max_permissible_scalar;
  Real scalar_lower_bound = min_permissible_scalar;
  _iteration = 0;

  _initial_residual = _residual = computeResidual(trial_value, scalar);

  Real residual_old = _residual;
  Real init_resid_sign = MathUtils::sign(_residual);
  Real reference_residual = computeReferenceResidual(trial_value, scalar);

  if (converged(_residual, reference_residual))
  {
    iterationFinalize(scalar);
    outputIterationStep(
        iter_output, _iteration, trial_value, scalar, _residual, reference_residual);
    return SolveState::SUCCESS;
  }

  _residual_history.assign(_num_resids, std::numeric_limits<Real>::max());
  _residual_history[0] = _residual;

  while (_iteration < _max_its && !converged(_residual, reference_residual) &&
         !convergedAcceptable(_iteration, _residual, reference_residual))
  {
    scalar_increment = -_residual / computeDerivative(trial_value, scalar);
    scalar = scalar_old + scalar_increment;

    if (_check_range)
      checkPermissibleRange(scalar,
                            scalar_increment,
                            scalar_old,
                            min_permissible_scalar,
                            max_permissible_scalar,
                            iter_output);

    _residual = computeResidual(trial_value, scalar);
    reference_residual = computeReferenceResidual(trial_value, scalar);
    iterationFinalize(scalar);

    if (_bracket_solution)
      updateBounds(
          scalar, _residual, init_resid_sign, scalar_upper_bound, scalar_lower_bound, iter_output);

    if (converged(_residual, reference_residual))
    {
      outputIterationStep(
          iter_output, _iteration, trial_value, scalar, _residual, reference_residual);
      break;
    }
    else
    {
      bool modified_increment = false;

      // Line Search
      if (_line_search)
      {
        if (residual_old - _residual != 0.0)
        {
          Real alpha = residual_old / (residual_old - _residual);
          alpha = MathUtils::clamp(alpha, 1.0e-2, 1.0);

          if (alpha != 1.0)
          {
            modified_increment = true;
            scalar_increment *= alpha;
            if (iter_output)
              *iter_output << "  Line search alpha = " << alpha
                           << " increment = " << scalar_increment << std::endl;
          }
        }
      }

      if (_bracket_solution)
      {
        // Check to see whether trial scalar_increment is outside the bounds, and set it to a point
        // within the bounds if it is
        if (scalar_old + scalar_increment >= scalar_upper_bound ||
            scalar_old + scalar_increment <= scalar_lower_bound)
        {
          if (scalar_upper_bound != max_permissible_scalar &&
              scalar_lower_bound != min_permissible_scalar)
          {
            const Real frac = 0.5;
            scalar_increment =
                (1.0 - frac) * scalar_lower_bound + frac * scalar_upper_bound - scalar_old;
            modified_increment = true;
            if (iter_output)
              *iter_output << "  Trial scalar_increment exceeded bounds.  Setting between "
                              "lower/upper bounds. frac: "
                           << frac << std::endl;
          }
        }
      }

      // Update the trial scalar and recompute residual if the line search or bounds checking
      // modified the increment
      if (modified_increment)
      {
        scalar = scalar_old + scalar_increment;
        _residual = computeResidual(trial_value, scalar);
        reference_residual = computeReferenceResidual(trial_value, scalar);
        iterationFinalize(scalar);

        if (_bracket_solution)
          updateBounds(scalar,
                       _residual,
                       init_resid_sign,
                       scalar_upper_bound,
                       scalar_lower_bound,
                       iter_output);
      }
    }

    outputIterationStep(
        iter_output, _iteration, trial_value, scalar, _residual, reference_residual);

    ++_iteration;
    residual_old = _residual;
    scalar_old = scalar;
    _residual_history[_iteration % _num_resids] = _residual;
  }

  if (std::isnan(_residual) || std::isinf(_residual))
    return SolveState::NAN_INF;

  if (_iteration == _max_its)
    return SolveState::EXCEEDED_ITERATIONS;

  return SolveState::SUCCESS;
}

bool
NewtonIteration::converged(const Real residual, const Real reference)
{
  return (std::abs(residual) <= _internal_iteration_abs_tol ||
          std::abs(residual / reference) <= _internal_iteration_rel_tol);
}

bool
NewtonIteration::convergedAcceptable(const unsigned int it,
                                                         const Real residual,
                                                         const Real reference)
{
  // Require that we have at least done _num_resids evaluations before we allow for
  // acceptable convergence
  if (it < _num_resids)
    return false;

  // Check to see whether the residual has dropped by convergence_history_factor over
  // the last _num_resids iterations. If it has (which means it's still making progress),
  // don't consider it to be converged within the acceptable limits.
  const Real convergence_history_factor = 10.0;
  if (std::abs(residual * convergence_history_factor) <
      std::abs(_residual_history[(it + 1) % _num_resids]))
    return false;

  // Now that it's determined that progress is not being made, treat it as converged if
  // we're within the acceptable convergence limits
  return converged(residual / _acceptable_multiplier, reference);
}

void
NewtonIteration::outputIterationStep(std::stringstream * iter_output,
                                                         const unsigned int it,
                                                         const Real trial_value,
                                                         const Real scalar,
                                                         const Real residual,
                                                         const Real reference_residual)
{
  if (iter_output)
  {
    *iter_output << " iteration=" << it << " Trial Value =" << trial_value
                 << " scalar=" << scalar << " residual=" << residual
                 << " ref_res=" << reference_residual
                 << " rel_res=" << std::abs(residual) / reference_residual
                 << " rel_tol=" << _internal_iteration_rel_tol << " abs_res=" << std::abs(residual)
                 << " abs_tol=" << _internal_iteration_abs_tol << '\n';
  }
}

void
NewtonIteration::outputIterationSummary(std::stringstream * iter_output,
                                                            const unsigned int total_it)
{
  if (iter_output)
    *iter_output << "In " << total_it << " iterations the residual went from " << _initial_residual
                 << " to " << _residual << " in '" << _svrms_name << "'.\n";
}

void
NewtonIteration::checkPermissibleRange(Real & scalar,
                                                           Real & scalar_increment,
                                                           const Real scalar_old,
                                                           const Real min_permissible_scalar,
                                                           const Real max_permissible_scalar,
                                                           std::stringstream * iter_output)
{
  if (scalar > max_permissible_scalar)
  {
    scalar_increment = (max_permissible_scalar - scalar_old) / 2.0;
    scalar = scalar_old + scalar_increment;
    if (iter_output)
      *iter_output << "Scalar greater than maximum (" << max_permissible_scalar
                   << ") adjusted scalar=" << scalar << " scalar_increment=" << scalar_increment
                   << std::endl;
  }
  else if (scalar < min_permissible_scalar)
  {
    scalar_increment = (min_permissible_scalar - scalar_old) / 2.0;
    scalar = scalar_old + scalar_increment;
    if (iter_output)
      *iter_output << "Scalar less than minimum (" << min_permissible_scalar
                   << ") adjusted scalar=" << scalar << " scalar_increment=" << scalar_increment
                   << std::endl;
  }
}

void
NewtonIteration::updateBounds(const Real scalar,
                                                  const Real residual,
                                                  const Real init_resid_sign,
                                                  Real & scalar_upper_bound,
                                                  Real & scalar_lower_bound,
                                                  std::stringstream * iter_output)
{
  // Update upper/lower bounds as applicable
  if (residual * init_resid_sign < 0.0 && scalar < scalar_upper_bound)
  {
    scalar_upper_bound = scalar;
    if (scalar_upper_bound < scalar_lower_bound)
    {
      scalar_upper_bound = scalar_lower_bound;
      scalar_lower_bound = 0.0;
      if (iter_output)
        *iter_output << "  Corrected for scalar_upper_bound < scalar_lower_bound" << std::endl;
    }
  }
  // Don't permit setting scalar_lower_bound > scalar_upper_bound (but do permit the reverse).
  // This ensures that if we encounter multiple roots, we pick the lowest one.
  else if (residual * init_resid_sign > 0.0 && scalar > scalar_lower_bound &&
           scalar < scalar_upper_bound)
    scalar_lower_bound = scalar;
}
