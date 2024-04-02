/*
  Copyright (C) 2023 by the authors of the ASPECT code.

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

// Plugin to enable remote control of ASPECT using stdin to enable
// an Inversion.

#include "aspect/material_model/interface.h"
#include "aspect/material_model/visco_plastic.h"
#include "aspect/material_model/rheology/visco_plastic.h"
#include "aspect/time_stepping/interface.h"
#include <aspect/simulator.h>
#include <aspect/utilities.h>
#include <cstddef>
#include <deal.II/base/mpi.h>
#include <deal.II/base/utilities.h>
#include <world_builder/world.h>
#include <deal.II/grid/tria.h>

#include <iostream>
#include <typeinfo>

using namespace aspect;

template <int dim>
class Inversion : public TimeStepping::Interface<dim>, public aspect::SimulatorAccess<dim>
{
  public:
    /**
     * Constructor.
     */
    Inversion () = default;

    /**
     * @copydoc aspect::TimeStepping::Interface<dim>::execute()
     */
    double
    execute() override;

    /**
     * The main execute() function.
     */
    std::pair<TimeStepping::Reaction, double>
    determine_reaction(const TimeStepping::TimeStepInfo &info) override;

    static
    void
    declare_parameters (dealii::ParameterHandler &prm);

    void
    parse_parameters (dealii::ParameterHandler &prm) override;

  private:

};


template <int dim>
double
Inversion<dim>::execute()
{
  this->get_pcout() << "Inversion<dim>::execute()" << std::endl;
  Simulator<dim> &sim = const_cast<Simulator<dim>&>(this->get_simulator());
  sim.postprocess();

  const bool am_i_rank_0 = dealii::Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0;

  while (true)
    {
      std::string line;
      if (am_i_rank_0)
        {
          std ::cout << "?" << std::endl;
          std::getline(std::cin, line);
        }

      line = dealii::Utilities::MPI::broadcast(this->get_mpi_communicator(), line, 0);

      const std::vector<std::string> parts = aspect::Utilities::split_string_list(line," ");

      if (parts[0] == "wb" && parts.size()==2)
        {
          this->get_pcout() << "loading WB " << parts[1] << std::endl;
          sim.world_builder = std::make_shared<WorldBuilder::World>(parts[1]);
        }
      else if ((parts.size() == 2 && ((parts[0] == "cohesion_factor") || (parts[0] == "friction_angle_factor") || (parts[0] == "strain_factor"))) ||
      (parts[0] == "continue"))
        {
          const double new_value = dealii::Utilities::string_to_double(parts[1]);
          MaterialModel::ViscoPlastic<dim> *material = dynamic_cast<MaterialModel::ViscoPlastic<dim>*>(
                                                         const_cast<MaterialModel::Interface<dim>*>(&this->get_material_model())
                                                       );

          if (material != nullptr && parts[0] == "cohesion_factor")
          {
            this->get_pcout() << "cohesion value: " << new_value << std::endl;
            material->strain_dependent_rheology->friction_strain_weakening_factors[0] = new_value;
          }
          else if (material != nullptr && parts[0] == "friction_angle_factor")
          {
            this->get_pcout() << "friction angle value: " << new_value << std::endl;
            material->strain_dependent_rheology->friction_strain_weakening_factors[0] = new_value;
          }        
          else if (material != nullptr && parts[0] == "strain_factor")
          {
            this->get_pcout() << "strain weakening value: " << new_value << std::endl;
            material->strain_dependent_rheology->viscous_strain_weakening_factors[0]  = new_value;
          }
            
        }
    }



  // re-initialize the initial condition plugins to find the new world builder
  for (auto &x: this->get_initial_composition_manager().get_active_initial_composition_conditions())
    x->initialize();

  for (auto &x: this->get_initial_temperature_manager().get_active_initial_temperature_conditions())
    x->initialize();

  return std::numeric_limits<double>::max();
}



template <int dim>
std::pair<TimeStepping::Reaction, double>
Inversion<dim>::determine_reaction (const TimeStepping::TimeStepInfo &/*info*/)
{
  // always repeat
  this->get_pcout() << "Inversion<dim>::determine_reaction()" << std::endl;
  return std::make_pair<TimeStepping::Reaction, double>
         (TimeStepping::Reaction::restart, this->get_timestep());
}



template <int dim>
void
Inversion<dim>::declare_parameters (ParameterHandler &prm)
{
  prm.enter_subsection("Time stepping");
  {
  }
  prm.leave_subsection();
}



template <int dim>
void
Inversion<dim>::parse_parameters (ParameterHandler &prm)
{
  prm.enter_subsection("Time stepping");
  prm.leave_subsection();
}


ASPECT_REGISTER_TIME_STEPPING_MODEL(Inversion,
                                    "inversion",
                                    "")
