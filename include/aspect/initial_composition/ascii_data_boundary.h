/*
  Copyright (C) 2011 - 2023 by the authors of the ASPECT code.

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


#ifndef _aspect_initial_composition_ascii_data_boundary_h
#define _aspect_initial_composition_ascii_data_boundary_h

#include <aspect/initial_composition/interface.h>

#include <aspect/simulator_access.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace InitialComposition
  {
    /**
     * A class that implements the prescribed compositional fields determined
     * from a AsciiData input file.
     *
     * @ingroup InitialCompositionModels
     */
    template <int dim>
    class AsciiDataBoundary : public Interface<dim>, public aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Empty Constructor.
         */
        AsciiDataBoundary ();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        void
        initialize () override;

        /**
         * Return the initial composition as a function of position. For the
         * current class, this function returns value from the text files.
         */
        double
        initial_composition (const Point<dim> &position,
                             const unsigned int n_comp) const override;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * The dataset that contains the initial composition defining the
         * layer depths.
         */
        unsigned int surface_boundary_id;

        /**
         * The dataset that contains the initial boundary layer data.
         */
        std::unique_ptr<Utilities::AsciiDataBoundary<dim>> ascii_data_boundary;

        double compositional_value;

        /**
         * Whether to use variable crustal densities taken as an input from
         * the ascii data boundary file.
         */
        bool use_variable_crustal_density;

        /**
         * Whether to use lithosphere depth from the ascii data boundary file.
         */
        bool use_lithosphere_depth;
    };
  }
}


#endif
