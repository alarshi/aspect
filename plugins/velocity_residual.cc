/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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
#include <aspect/simulator.h>
#include "velocity_residual.h"

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      VelocityResidual<dim>::
      VelocityResidual ()
        :
        DataPostprocessorVector<dim> ("velocity_residual",
        			update_values | update_quadrature_points | update_gradients)
      {}

      template <int dim>
      void
	  VelocityResidual<dim>::initialize ()
	  {
	    std::set<types::boundary_id> surface_boundary_set;
	    surface_boundary_set.insert(1); //*** hard coded
	// The input ascii table contains one data column (velocity components) in addition to the coordinate columns.
	    Utilities::AsciiDataBoundary<dim>::initialize(surface_boundary_set, dim);
	  }


      template <int dim>
      void
      VelocityResidual<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double> > &computed_quantities) const
      {

    	Assert ((computed_quantities[0].size() == dim), ExcInternalError());
        auto cell = input_data.template get_cell<DoFHandler<dim> >();

        const double velocity_scaling_factor =
                      this->convert_output_to_years() ? year_in_seconds : 1.0;


        // We only want to output dynamic topography at the top and bottom
        // boundary, so only compute it if the current cell has
        // a face at the top or bottom boundary.
        bool cell_at_top_boundary = false;
        for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
          if (cell->at_boundary(f) &&
              (this->get_geometry_model().translate_id_to_symbol_name (cell->face(f)->boundary_id()) == "top"))
            cell_at_top_boundary = true;

        if (cell_at_top_boundary)
		{
		  for (unsigned int q=0; q<computed_quantities.size(); ++q)
		  {
		    for (unsigned int d = 0; d < dim; ++d)
			{
		      Tensor<1,dim> velocity;

			  velocity[d] = Utilities::AsciiDataBoundary<dim>::get_data_component(
										  this->get_geometry_model().translate_symbolic_boundary_name_to_id("top"),
										  input_data.evaluation_points[q], d);
			  computed_quantities[q](d) = velocity[d] - input_data.solution_values[q][d] * velocity_scaling_factor;
			}
		  }
		}
      }


      template <int dim>
      void
	  VelocityResidual<dim>::declare_parameters (ParameterHandler &prm)
	  {
    	prm.enter_subsection("Postprocess");
    	{
          prm.enter_subsection("Visualization");
    	  {
		    Utilities::AsciiDataBase<dim>::declare_parameters(prm,
						  "$ASPECT_SOURCE_DIR/data/initial-temperature/adiabatic-boundary/",
						  "adiabatic_boundary.txt", "Surface properties");
    	  }
    	  prm.leave_subsection();
    	  }
    	prm.leave_subsection();
	  }


      template <int dim>
      void
	  VelocityResidual<dim>::parse_parameters (ParameterHandler &prm)
	  {
	    prm.enter_subsection("Postprocess");
		  prm.enter_subsection("Visualization");
		  {
			Utilities::AsciiDataBase<dim>::parse_parameters(prm, "Surface properties");
		  }
		  prm.leave_subsection();
	    prm.leave_subsection();
	  }


    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(VelocityResidual,
                                                  "velocity residual",
                                                  "A visualization for caculating the surface velocity to compare"
												  "with the GPS velocities.")
    }
  }
}
