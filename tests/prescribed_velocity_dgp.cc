#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/fe/fe_values.h>
#include <aspect/global.h>
#include <aspect/simulator_signals.h>

namespace aspect
{
  using namespace dealii;

  // Global variables (to be set by parameters)
  bool prescribe_internal_velocities;

  // Because we do not initially know what dimension we're in, we need
  // function parser objects for both 2d and 3d.
  Functions::ParsedFunction<2> prescribed_velocity_indicator_function_2d (2);
  Functions::ParsedFunction<3> prescribed_velocity_indicator_function_3d (3);
  Functions::ParsedFunction<2> prescribed_velocity_function_2d (2);
  Functions::ParsedFunction<3> prescribed_velocity_function_3d (3);

  /**
   * Declare additional parameters.
   */
  void declare_parameters(const unsigned int dim,
                          ParameterHandler &prm)
  {
    prm.declare_entry ("Prescribe internal velocities", "false",
                       Patterns::Bool (),
                       "Whether or not to use any prescribed internal velocities. "
                       "Locations in which to prescribe velocities are defined "
                       "in section ``Prescribed velocities/Indicator function'' "
                       "and the velocities are defined in section ``Prescribed "
                       "velocities/Velocity function''. Indicators are evaluated "
                       "at the center of each cell, and all DOFs associated with "
                       "the specified velocity component at the indicated cells "
                       "are constrained."
                      );

    prm.enter_subsection ("Prescribed velocities");
    {
      prm.enter_subsection ("Indicator function");
      {
        if (dim == 2)
          Functions::ParsedFunction<2>::declare_parameters (prm, 2);
        else
          Functions::ParsedFunction<3>::declare_parameters (prm, 3);
      }
      prm.leave_subsection ();

      prm.enter_subsection ("Velocity function");
      {
        if (dim == 2)
          Functions::ParsedFunction<2>::declare_parameters (prm, 2);
        else
          Functions::ParsedFunction<3>::declare_parameters (prm, 3);
      }
      prm.leave_subsection ();
    }
    prm.leave_subsection ();
  }

  template <int dim>
  void parse_parameters(const Parameters<dim> &,
                        ParameterHandler &prm)
  {
    prescribe_internal_velocities = prm.get_bool ("Prescribe internal velocities");
    prm.enter_subsection ("Prescribed velocities");
    {
      prm.enter_subsection("Indicator function");
      {
        try
          {
            if (dim == 2)
              prescribed_velocity_indicator_function_2d.parse_parameters (prm);
            else
              prescribed_velocity_indicator_function_3d.parse_parameters (prm);
          }
        catch (...)
          {
            std::cerr << "ERROR: FunctionParser failed to parse\n"
                      << "\t'Prescribed velocities.Indicator function'\n"
                      << "with expression\n"
                      << "\t'" << prm.get("Function expression") << "'";
            throw;
          }
      }
      prm.leave_subsection();

      prm.enter_subsection("Velocity function");
      {
        try
          {
            if (dim == 2)
              prescribed_velocity_function_2d.parse_parameters (prm);
            else
              prescribed_velocity_function_3d.parse_parameters (prm);
          }
        catch (...)
          {
            std::cerr << "ERROR: FunctionParser failed to parse\n"
                      << "\t'Prescribed velocities.Velocity function'\n"
                      << "with expression\n"
                      << "\t'" << prm.get("Function expression") << "'";
            throw;
          }
      }
      prm.leave_subsection();
    }
    prm.leave_subsection ();
  }

  /**
   * This function retrieves the unit support points (in the unit cell) for the current element.
   * The DGP element used when 'set Use locally conservative discretization = true' does not
   * have support points. If these elements are in use, a fictitious support point at the cell
   * center is returned for each shape function that corresponds to the pressure variable,
   * whereas the support points for the velocity are correct; the fictitious points don't matter
   * because we only use this function when interpolating the velocity variable, and ignore the
   * evaluation at the pressure support points.
   */
  template <int dim>
  std::vector< Point<dim> >
  get_unit_support_points_for_velocity(const SimulatorAccess<dim> &simulator_access)
  {
    std::vector< Point<dim> > unit_support_points;
    if ( !simulator_access.get_parameters().use_locally_conservative_discretization )
      {
        return simulator_access.get_fe().get_unit_support_points();
      }
    else
      {
        //special case for discontinuous pressure elements, which lack unit support points
        std::vector< Point<dim> > unit_support_points;
        const unsigned int dofs_per_cell = simulator_access.get_fe().dofs_per_cell;
        for (unsigned int dof=0; dof < dofs_per_cell; ++dof)
          {
            // base will hold element, base_index holds node/shape function within that element
            const unsigned int
            base       = simulator_access.get_fe().system_to_base_index(dof).first.first,
            base_index = simulator_access.get_fe().system_to_base_index(dof).second;
            // get the unit support points for the relevant element
            std::vector< Point<dim> > my_support_points = simulator_access.get_fe().base_element(base).get_unit_support_points();
            if ( my_support_points.size() == 0 )
              {
                //manufacture a support point, arbitrarily at cell center
                if ( dim==2 )
                  unit_support_points.push_back( Point<dim> (0.5,0.5) );
                if ( dim==3 )
                  unit_support_points.push_back( Point<dim> (0.5,0.5,0.5) );
              }
            else
              {
                unit_support_points.push_back(my_support_points[base_index]);
              }
          }
        return unit_support_points;
      }
  }

  /**
   * This function is called by a signal which is triggered after the other constraints
   * have been calculated. This enables us to define additional constraints in the mass
   * matrix on any arbitrary degree of freedom in the model space.
   */
  template <int dim>
  void constrain_internal_velocities (const SimulatorAccess<dim> &simulator_access,
                                      AffineConstraints<double> &current_constraints)
  {
    if (prescribe_internal_velocities)
      {
        const std::vector< Point<dim> > points = get_unit_support_points_for_velocity(simulator_access);
        const Quadrature<dim> quadrature (points);
        FEValues<dim> fe_values (simulator_access.get_fe(), quadrature, update_quadrature_points);
        typename DoFHandler<dim>::active_cell_iterator cell;

        // Loop over all cells
        for (cell = simulator_access.get_dof_handler().begin_active();
             cell != simulator_access.get_dof_handler().end();
             ++cell)
          if (! cell->is_artificial())
            {
              fe_values.reinit (cell);
              std::vector<types::global_dof_index> local_dof_indices(simulator_access.get_fe().dofs_per_cell);
              cell->get_dof_indices (local_dof_indices);

              for (unsigned int q=0; q<quadrature.size(); q++)
                // If it's okay to constrain this DOF
                if (current_constraints.can_store_line(local_dof_indices[q]) &&
                    !current_constraints.is_constrained(local_dof_indices[q]))
                  {
                    // Get the velocity component index
                    const unsigned int c_idx =
                      simulator_access.get_fe().system_to_component_index(q).first;

                    // If we're on one of the velocity DOFs
                    if ((c_idx >=
                         simulator_access.introspection().component_indices.velocities[0])
                        &&
                        (c_idx <=
                         simulator_access.introspection().component_indices.velocities[dim-1]))
                      {
                        // Which velocity component is this DOF associated with?
                        const unsigned int component_direction
                          = (c_idx
                             - simulator_access.introspection().component_indices.velocities[0]);

                        // we get time passed as seconds (always) but may want
                        // to reinterpret it in years
                        double time = simulator_access.get_time();
                        if (simulator_access.convert_output_to_years())
                          time /= year_in_seconds;

                        prescribed_velocity_indicator_function_2d.set_time (time);
                        prescribed_velocity_indicator_function_3d.set_time (time);
                        prescribed_velocity_function_2d.set_time (time);
                        prescribed_velocity_function_3d.set_time (time);

                        const Point<dim> p = fe_values.quadrature_point(q);

                        // Because we defined and parsed our parameter
                        // file differently for 2d and 3d we need to
                        // be sure to query the correct object for
                        // function values. The function parser
                        // objects expect points of a certain
                        // dimension, but Point p will be compiled for
                        // both 2d and 3d, so we need some careful
                        // casts to make this compile. Casting Point
                        // objects between 2d and 3d is somewhat
                        // meaningless but the offending casts will
                        // never actually be executed because they are
                        // protected by the if/else statements.
                        double indicator, u_i;
                        if (dim == 2)
                          {
                            indicator = prescribed_velocity_indicator_function_2d.value
                                        (reinterpret_cast<const Point<2>&>(p),
                                         component_direction);
                            u_i       = prescribed_velocity_function_2d.value
                                        (reinterpret_cast<const Point<2>&>(p),
                                         component_direction);
                          }
                        else
                          {
                            indicator = prescribed_velocity_indicator_function_3d.value
                                        (reinterpret_cast<const Point<3>&>(p),
                                         component_direction);
                            u_i       = prescribed_velocity_function_3d.value
                                        (reinterpret_cast<const Point<3>&>(p),
                                         component_direction);
                          }

                        if (indicator > 0.5)
                          {
                            // Add a constraint of the form dof[q] = u_i
                            // to the list of constraints.
                            current_constraints.add_line (local_dof_indices[q]);
                            current_constraints.set_inhomogeneity (local_dof_indices[q], u_i);
                          }
                      }
                  }
            }
      }
  }

  // Connect declare_parameters and parse_parameters to appropriate signals.
  void parameter_connector ()
  {
    SimulatorSignals<2>::declare_additional_parameters.connect (&declare_parameters);
    SimulatorSignals<3>::declare_additional_parameters.connect (&declare_parameters);

    SimulatorSignals<2>::parse_additional_parameters.connect (&parse_parameters<2>);
    SimulatorSignals<3>::parse_additional_parameters.connect (&parse_parameters<3>);
  }

  // Connect constraints function to correct signal.
  template <int dim>
  void signal_connector (SimulatorSignals<dim> &signals)
  {
    signals.post_constraints_creation.connect (&constrain_internal_velocities<dim>);
  }

  // Tell ASPECT to send signals to the connector functions
  ASPECT_REGISTER_SIGNALS_PARAMETER_CONNECTOR(parameter_connector)
  ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>, signal_connector<3>)
}
