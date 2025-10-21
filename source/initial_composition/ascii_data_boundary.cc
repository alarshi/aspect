/*
  Copyright (C) 2011 - 2025 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
 */


#include <aspect/global.h>
#include <aspect/initial_composition/ascii_data_boundary.h>

#include <deal.II/base/parameter_handler.h>



namespace aspect
{
  namespace InitialComposition
  {
    template <int dim>
    AsciiDataBoundary<dim>::AsciiDataBoundary ()
      = default;


    template <int dim>
    void
    AsciiDataBoundary<dim>::initialize ()
    {
      // Currently, the plugin is only implemented such that the boundary is computed from the top surface
      // of the geometry model
      surface_boundary_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");

      if (use_variable_crustal_density && use_lithosphere_depth)
        ascii_data_boundary->initialize({surface_boundary_id}, 3);
      else
        ascii_data_boundary->initialize({surface_boundary_id}, 1);
    }


    template <int dim>
    double
    AsciiDataBoundary<dim>::
    initial_composition (const Point<dim> &position,
                         const unsigned int n_comp) const
    {
      const double depth = this->get_geometry_model().depth(position);

      // Here, we assume that the first column after the boundary coordinates is the depth to the layer
      const double crustal_depth = ascii_data_boundary->get_data_component(surface_boundary_id, position, 0);

      if (depth <= crustal_depth && n_comp == 0)
        return ascii_data_boundary->get_data_component(surface_boundary_id, position, 1);

      // The following will only work if we have two compositions and here we assume that the second
      // composition is the lithosphere
      if (use_lithosphere_depth && n_comp == 1)
        {
          // compute the depth to the base of the lithosphere, this is the 3rd column after the boundary coordinates
          const double lithosphere_depth = ascii_data_boundary->get_data_component(surface_boundary_id, position, 2);

          if (depth <= lithosphere_depth && depth > crustal_depth)
            return compositional_value;
        }


      return 0.;
    }


    template <int dim>
    void
    AsciiDataBoundary<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("Ascii data boundary model");
        {
          Utilities::AsciiDataBoundary<dim>::declare_parameters(prm,
                                                                "$ASPECT_SOURCE_DIR/data/initial-temperature/adiabatic-boundary/",
                                                                "adiabatic_boundary.txt", "Crustal depths");

          prm.declare_entry ("Compositional value", "1",
                             Patterns::List (Patterns::Double(0)),
                             "The value assigned to the compositional field in between the surface of the model and the layer "
                             "defined by the ascii data boundary file."
                             "Units: none.");
          prm.declare_entry ("Use variable crustal densities", "false",
                             Patterns::Bool (),
                             "Whether to use variable crustal densities taken as an input from as "
                             "ascii data boundary file. ");
          prm.declare_entry ("Use lithosphere depth", "false",
                             Patterns::Bool (),
                             "Whether to use the lithosphere depth from the ascii data boundary file. "
                             "If true, the third column in the ascii data file is interpreted as the lithosphere depth. "
                             "If false, the third column is ignored and only the second column is used for compositional value.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    AsciiDataBoundary<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("Ascii data boundary model");
        {
          ascii_data_boundary = std::make_unique<Utilities::AsciiDataBoundary<dim>>();
          ascii_data_boundary->initialize_simulator(this->get_simulator());
          ascii_data_boundary->parse_parameters(prm, "Crustal depths");
          compositional_value = prm.get_double ("Compositional value");
          use_variable_crustal_density = prm.get_bool ("Use variable crustal densities");
          use_lithosphere_depth = prm.get_bool ("Use lithosphere depth");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace InitialComposition
  {
    ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(AsciiDataBoundary,
                                              "ascii data boundary",
                                              "Implementation of a model in which the composition "
                                              "is defined in a layer from the surface to the depth "
                                              "defined in the input ascii data. A constant compositional "
                                              "value is assigned to the layer, which is a user-defined "
                                              "parameter. "
                                              "Note the required format of the input data: The first lines "
                                              "may contain any number of comments if they begin with "
                                              "`#', but one of these lines needs to "
                                              "contain the number of grid points in each dimension as "
                                              "for example `# POINTS: 3 3'. "
                                              "The order of the data columns "
                                              "has to be `x', `Depth [m]' in a 2d model and "
                                              " `x', `y', `Depth [m]' in a 3d model, which means that "
                                              "there has to be a single column "
                                              "containing the depth. "
                                              "Note that the data in the input "
                                              "files need to be sorted in a specific order: "
                                              "the first coordinate needs to ascend first, "
                                              "followed by the second in order to "
                                              "assign the correct data to the prescribed coordinates. "
                                              "If you use a spherical model, "
                                              "then the assumed grid changes. `x' will be replaced by "
                                              "the radial distance of the point to the bottom of the model, "
                                              "`y' by the azimuth angle and `z' by the polar angle measured "
                                              "positive from the north pole. The grid will be assumed to be "
                                              "a latitude-longitude grid. Note that the order "
                                              "of spherical coordinates is `r', `phi', `theta' "
                                              "and not `r', `theta', `phi', since this allows "
                                              "for dimension independent expressions.")
  }
}
