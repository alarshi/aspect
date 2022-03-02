/*
  Copyright (C) 2014 - 2018 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include "equilibrium_grain_size.h"
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/utilities.h>
#include <aspect/initial_temperature/interface.h>
#include <aspect/initial_temperature/adiabatic_boundary.h>
#include <aspect/simulator_signals.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/signaling_nan.h>

#include <iostream>

using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {
    /**
     * Additional output fields for the viscosity prior to scaling to be added to
     * the MaterialModel::MaterialModelOutputs structure and filled in the
     * MaterialModel::Interface::evaluate() function.
     */
    template <int dim>
    class UnscaledViscosityAdditionalOutputs : public NamedAdditionalMaterialOutputs<dim>
    {
      public:
        UnscaledViscosityAdditionalOutputs(const unsigned int n_points)
          : NamedAdditionalMaterialOutputs<dim>(std::vector<std::string>(1, "unscaled_viscosity"),
                                                n_points)
        {}
    };

    /**
     * Additional output fields for the the material type, decribing if we are in the
     * crust/lithosphere/asthenosphere/lower mantle.
     */
    template <int dim>
    class MaterialTypeAdditionalOutputs : public NamedAdditionalMaterialOutputs<dim>
    {
      public:
        MaterialTypeAdditionalOutputs(const unsigned int n_points)
          : NamedAdditionalMaterialOutputs<dim>(std::vector<std::string>(1, "material_type"),
                                                n_points)
        {}
    };
  }

  namespace internal
  {
    template <int dim>
    class FunctorDepthAverageUnscaledViscosity: public internal::FunctorBase<dim>
    {
      public:
        FunctorDepthAverageUnscaledViscosity()
        {}

        bool need_material_properties() const override
        {
          return true;
        }

        void
        create_additional_material_model_outputs (const unsigned int n_points,
                                                  MaterialModel::MaterialModelOutputs<dim> &outputs) const override
        {
          if (outputs.template get_additional_output<MaterialModel::UnscaledViscosityAdditionalOutputs<dim>>() == nullptr)
            {
              outputs.additional_outputs.push_back(
                std::make_unique<MaterialModel::UnscaledViscosityAdditionalOutputs<dim>> (n_points));
            }
        }

        void operator()(const MaterialModel::MaterialModelInputs<dim> &,
                        const MaterialModel::MaterialModelOutputs<dim> &out,
                        const FEValues<dim> &,
                        const LinearAlgebra::BlockVector &,
                        std::vector<double> &output) override
        {
          const MaterialModel::UnscaledViscosityAdditionalOutputs<dim> *unscaled_viscosity_outputs
            = out.template get_additional_output<const MaterialModel::UnscaledViscosityAdditionalOutputs<dim> >();

          Assert(unscaled_viscosity_outputs != nullptr,ExcInternalError());

          for (unsigned int q=0; q<output.size(); ++q)
            output[q] = unscaled_viscosity_outputs->output_values[0][q];
        }
    };
  }

  namespace MaterialModel
  {
    template <int dim>
    void
    EquilibriumGrainSize<dim>::initialize()
    {
      // Get reference viscosity profile from the ascii data
      if (use_depth_dependent_viscosity)
        reference_viscosity_coordinates = reference_viscosity_profile->get_interpolation_point_coordinates();

      // Get column index for density scaling
      if (use_depth_dependent_rho_vs)
        {
          rho_vs_depth_profile.initialize(this->get_mpi_communicator());
          density_scaling_index = rho_vs_depth_profile.get_column_index_from_name("density_scaling");
        }

      if (use_depth_dependent_thermal_expansivity)
        {
          thermal_expansivity_profile.initialize(this->get_mpi_communicator());
          thermal_expansivity_column_index = thermal_expansivity_profile.get_column_index_from_name("thermal_expansivity");
        }

      if (use_depth_dependent_dT_vs)
        {
          dT_vs_depth_profile.initialize(this->get_mpi_communicator());
          temperature_scaling_index = dT_vs_depth_profile.get_column_index_from_name("temperature_scaling");
        }

      // Get column for crustal depths
      std::set<types::boundary_id> surface_boundary_set;
      surface_boundary_set.insert(this->get_geometry_model().translate_symbolic_boundary_name_to_id("top"));
      crustal_boundary_depth.initialize(surface_boundary_set, 1);

      this->get_signals().post_stokes_solver.connect([&](const SimulatorAccess<dim> &,
                                                         const unsigned int ,
                                                         const unsigned int ,
                                                         const SolverControl &,
                                                         const SolverControl &)
      {
        this->update();
      });

      this->get_signals().post_advection_solver.connect([&](const SimulatorAccess<dim> &,
                                                            const unsigned int ,
                                                            const unsigned int ,
                                                            const SolverControl &)
      {
        this->update();
      });

      n_material_data = material_file_names.size();
      for (unsigned i = 0; i < n_material_data; i++)
        {
          if (material_file_format == perplex)
            material_lookup.push_back(std::shared_ptr<MaterialUtilities::Lookup::MaterialLookup>
                                      (new MaterialUtilities::Lookup::PerplexReader(datadirectory+material_file_names[i],
                                                                                    /*use_bilinear_interpolation*/ true,
                                                                                    this->get_mpi_communicator())));
          else if (material_file_format == hefesto)
            material_lookup.push_back(std::shared_ptr<MaterialUtilities::Lookup::MaterialLookup>
                                      (new MaterialUtilities::Lookup::HeFESToReader(datadirectory+material_file_names[i],
                                                                                    datadirectory+derivatives_file_names[i],
                                                                                    /*use_bilinear_interpolation*/ true,
                                                                                    this->get_mpi_communicator())));
          else
            AssertThrow (false, ExcNotImplemented());
        }
    }



    template <int dim>
    void
    EquilibriumGrainSize<dim>::update()
    {
      if (use_depth_dependent_viscosity)
        {
          std::vector<std::unique_ptr<internal::FunctorBase<dim> > > lateral_averaging_properties;
          lateral_averaging_properties.push_back(std::make_unique<internal::FunctorDepthAverageUnscaledViscosity<dim>>());

          std::vector<std::vector<double>> averages =
                                          this->get_lateral_averaging().compute_lateral_averages(reference_viscosity_coordinates,
                                              lateral_averaging_properties);

          average_viscosity_profile.swap(averages[0]);

          for (const auto &lateral_viscosity_average: average_viscosity_profile)
            AssertThrow(numbers::is_finite(lateral_viscosity_average),
                        ExcMessage("In computing depth averages, there is at"
                                   " least one depth band that does not have"
                                   " any quadrature points in it."
                                   " Consider reducing number of depth layers"
                                   " for averaging."));
        }

      initialized = true;
    }



    template <int dim>
    std::pair< double, unsigned int>
    EquilibriumGrainSize<dim>::get_reference_viscosity (const double depth) const
    {
      // Make maximal depth slightly larger to ensure depth < maximal_depth
      const double maximal_depth = this->get_geometry_model().maximal_depth() *
                                   (1.0+std::numeric_limits<double>::epsilon());
      (void) maximal_depth;

      Assert(depth < maximal_depth, ExcInternalError());

      unsigned int depth_index;
      if (depth < reference_viscosity_coordinates.front())
        {
          depth_index = 0;
        }
      else if (depth > reference_viscosity_coordinates.back())
        {
          depth_index = reference_viscosity_coordinates.size() - 1;
        }
      else
        {
          depth_index = std::distance(reference_viscosity_coordinates.begin(),
                                      std::lower_bound(reference_viscosity_coordinates.begin(),
                                                       reference_viscosity_coordinates.end(),
                                                       depth));
        }
      if (depth_index > 0)
        depth_index -= 1;

      // When evaluating reference viscosity, evaluate at the next lower depth that is stored
      // in the reference profile instead of the actual depth. This makes the profile piecewise
      // constant. This will be specific to the viscosity profile used (and ignore the entry with
      // the largest depth in the profile).
      double reference_viscosity = reference_viscosity_profile->compute_viscosity(reference_viscosity_coordinates.at(depth_index));

      // This parameter is only because we change the asthenosphere viscosity in our models.
      // By default, it is set to the value in reference steinberger profile.
      if (depth_index == 1)
        reference_viscosity = asthenosphere_viscosity;

      return std::make_pair (reference_viscosity, depth_index);

    }



    template <int dim>
    double
    EquilibriumGrainSize<dim>::get_uppermost_mantle_thickness () const
    {
      return uppermost_mantle_thickness;
    }



    template <int dim>
    double
    EquilibriumGrainSize<dim>::compute_viscosity_scaling (const double depth) const
    {
      Assert(average_viscosity_profile.size() != 0,
             ExcMessage("The average viscosity profile has not yet been computed. "
                        "Unable to scale viscosities"));

      const std::pair<double, unsigned int> reference_viscosity_and_depth_index = get_reference_viscosity (depth);

      const double average_viscosity = std::pow(10, average_viscosity_profile[reference_viscosity_and_depth_index.second]);

      return reference_viscosity_and_depth_index.first / average_viscosity;
    }



    template <int dim>
    void
    EquilibriumGrainSize<dim>::compute_equilibrium_grain_size(const typename Interface<dim>::MaterialModelInputs &in,
                                                              typename Interface<dim>::MaterialModelOutputs &out) const
    {
      PrescribedFieldOutputs<dim> *prescribed_field_out = out.template get_additional_output<PrescribedFieldOutputs<dim> >();
      DislocationViscosityOutputs<dim> *disl_viscosities_out = out.template get_additional_output<DislocationViscosityOutputs<dim> >();

      UnscaledViscosityAdditionalOutputs<dim> *unscaled_viscosity_out =
        out.template get_additional_output<MaterialModel::UnscaledViscosityAdditionalOutputs<dim> >();

      const unsigned int craton_composition = this->introspection().compositional_index_for_name("continents");

      const InitialTemperature::AdiabaticBoundary<dim> &adiabatic_boundary =
        this->get_initial_temperature_manager().template get_matching_initial_temperature_model<InitialTemperature::AdiabaticBoundary<dim> >();

      const unsigned int surface_boundary_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id("outer");

      const unsigned int grain_size_index = this->introspection().compositional_index_for_name("grain_size");
      const unsigned int fault_index = (use_faults)
                                       ?
                                       this->introspection().compositional_index_for_name("faults")
                                       :
                                       numbers::invalid_unsigned_int;

      for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
        {
          const double depth = this->get_geometry_model().depth(in.position[i]);

          // Use the adiabatic pressure instead of the real one, because of oscillations
          const double pressure = (this->get_adiabatic_conditions().is_initialized())
                                  ?
                                  this->get_adiabatic_conditions().pressure(in.position[i])
                                  :
                                  in.pressure[i];

          Assert(pressure >= 0.0,
                 ExcMessage("Pressure has to be non-negative for the viscosity computation. Instead it is: "
                            + std::to_string(pressure)));

          double lithosphere_thickness = 0.;
          // Get variable lithosphere using an adiabatic boundary ascii file
          if (this->get_adiabatic_conditions().is_initialized())
            lithosphere_thickness = adiabatic_boundary.get_data_component(surface_boundary_id, in.position[i], 0);

          const unsigned int phase_index = get_phase_index(in.position[i], in.temperature[i], pressure);

          // Computed according to equation (7) in Dannberg et al., 2016, using the paleowattmeter grain size.
          // Austin and Evans (2007): Paleowattmeters: A scaling relation for dynamically recrystallized grain size. Geology 35, 343-346.
          const double prefactor = geometric_constant[phase_index] * grain_boundary_energy[phase_index] * grain_growth_rate_constant[phase_index]
                                   / (boundary_area_change_work_fraction[phase_index] * grain_growth_exponent[phase_index]);
          const double exponential = std::exp(- (grain_growth_activation_energy[phase_index] + pressure * grain_growth_activation_volume[phase_index])
                                              / (constants::gas_constant * in.temperature[i]));

          // get the dislocation viscosity to compute the dislocation strain rate
          // If we do not have the strain rate yet, set it to a low value.
          // TODO: make minimum strain rate an input parameter
          const SymmetricTensor<2,dim> shear_strain_rate = in.strain_rate[i] - 1./dim * trace(in.strain_rate[i]) * unit_symmetric_tensor<dim>();
          const double second_strain_rate_invariant = std::sqrt(std::abs(second_invariant(shear_strain_rate)));

          // TODO: if we update the interface of the diffusion_viscosity and dislocation_viscosity functions,
          // we don't need this vector anymore
          std::vector<double> composition (in.composition[i]);
          double grain_size = in.composition[i][this->introspection().compositional_index_for_name("grain_size")];
          grain_size = std::max(min_grain_size, grain_size);

          // Only consider dislocation creep and equilibrium grain size if we have a sufficient strain rate
          if (std::abs(second_strain_rate_invariant > 1e-30))
            {
              unsigned int j = 0;
              double old_grain_size = 0.0;

              // because the diffusion viscosity depends on the grain size itself, and we need it to compute the dislocation strain rate,
              // we have to iterate in the computation of the equilibrium grain size
              while ((std::abs((grain_size-old_grain_size) / grain_size) > dislocation_viscosity_iteration_threshold)
                     && (j < dislocation_viscosity_iteration_number))
                {
                  composition[grain_size_index] = std::max(min_grain_size, grain_size);
                  const double diff_viscosity = diffusion_viscosity(in.temperature[i], pressure, composition, in.strain_rate[i], in.position[i]);
                  const double disl_viscosity = dislocation_viscosity(in.temperature[i], pressure, composition, in.strain_rate[i], in.position[i], diff_viscosity);

                  if (disl_viscosities_out != NULL)
                    {
                      disl_viscosities_out->dislocation_viscosities[i] = std::min(std::max(min_eta,disl_viscosity),1e30);
                      disl_viscosities_out->boundary_area_change_work_fractions[i] =
                        boundary_area_change_work_fraction[phase_index];
                    }

                  out.viscosities[i] = diff_viscosity * disl_viscosity / (disl_viscosity + diff_viscosity);

                  Assert(out.viscosities[i] > 0.0,
                         ExcMessage("Negative viscosity is not allowed. Current viscosity is: " + std::to_string(out.viscosities[i])));

                  // This follows from Equation (S25 - S30)
                  const double dislocation_strain_rate_invariant = second_strain_rate_invariant
                                                                   * out.viscosities[i] / disl_viscosity;
                  const double stress_term = 4.0 * out.viscosities[i] * second_strain_rate_invariant * dislocation_strain_rate_invariant;

                  old_grain_size = grain_size;

                  if (equilibrate_grain_size)
                    grain_size = 0.9 * old_grain_size + 0.1 * std::pow(prefactor/stress_term * exponential,1./(1+grain_growth_exponent[phase_index]));

                  ++j;
                }
            }
          else
            {
              out.viscosities[i] = diffusion_viscosity(in.temperature[i], pressure, composition, in.strain_rate[i], in.position[i]);

              if (disl_viscosities_out != NULL)
                {
                  disl_viscosities_out->dislocation_viscosities[i] = 1e30;
                  disl_viscosities_out->boundary_area_change_work_fractions[i] =
                    boundary_area_change_work_fraction[phase_index];
                }
            }

          if (unscaled_viscosity_out != nullptr)
            unscaled_viscosity_out->output_values[0][i] = std::log10(out.viscosities[i]);

          if (use_depth_dependent_viscosity)
            {
              const double viscosity_scaling_below_this_depth = 60e3;

              // Scale viscosity so that laterally averaged viscosity == reference viscosity profile
              // Only scale if average viscosity is already available and we are below a specified depth.
              if (average_viscosity_profile.size() != 0 && depth > viscosity_scaling_below_this_depth)
                out.viscosities[i] *= compute_viscosity_scaling(this->get_geometry_model().depth(in.position[i]));
            }

          // Ensure we respect viscosity bounds
          out.viscosities[i] = std::min(std::max(min_eta, out.viscosities[i]),max_eta);

          // If using faults, use the composition value to compute the viscosity instead
          // We extend the faults to depths 40 km more than the lithospheric depths because otherwise
          // faults do not extend through the lithosphere in our high-resolution models.
          if (use_faults && in.composition[i][fault_index] > 0. && depth <= lithosphere_thickness + 40e3)
            {
              const double background_viscosity_log = std::log10(out.viscosities[i]);
              out.viscosities[i] = std::pow(10,
                                            std::log10(fault_viscosity) * in.composition[i][fault_index]
                                            + background_viscosity_log * (1. - in.composition[i][fault_index]));
            }

          // If using cratons, use strong viscosity
          if (in.composition[i][craton_composition] > 0. && depth <= lithosphere_thickness + 40e3)
            out.viscosities[i] = 1e26;

          Assert(out.viscosities[i] > 0,
                 ExcMessage("Viscosity has to be positive. Instead it is: " + std::to_string(out.viscosities[i])));

          // Fill the prescribed outputs for grain size and assign faults to a prescribed field for diffusion.
          if (prescribed_field_out != NULL)
            for (unsigned int c=0; c<composition.size(); ++c)
              {
                if (c == grain_size_index)
                  prescribed_field_out->prescribed_field_outputs[i][c] = std::max(min_grain_size, grain_size);
                else if (c == fault_index)
                  prescribed_field_out->prescribed_field_outputs[i][c] = in.composition[i][fault_index];
                else
                  prescribed_field_out->prescribed_field_outputs[i][c] = 0.;
              }

          if (this->get_nonlinear_iteration() == 0 && use_depth_dependent_viscosity)
            out.viscosities[i] = get_reference_viscosity(depth).first;
        }
      return;
    }



    template <int dim>
    double
    EquilibriumGrainSize<dim>::
    phase_function (const Point<dim> &position,
                    const double temperature,
                    const double pressure,
                    const unsigned int phase) const
    {
      Assert(phase < transition_depths.size(),
             ExcMessage("Error: Phase index is too large. This phase index does not exist!"));

      // if we already have the adiabatic conditions, we can use them
      if (this->get_adiabatic_conditions().is_initialized())
        {
          // first, get the pressure at which the phase transition occurs normally
          const Point<dim,double> transition_point = this->get_geometry_model().representative_point(transition_depths[phase]);
          const double transition_pressure = this->get_adiabatic_conditions().pressure(transition_point);

          // then calculate the deviation from the transition point (both in temperature
          // and in pressure)
          const double pressure_deviation = pressure - transition_pressure
                                            - transition_slopes[phase] * (temperature - transition_temperatures[phase]);

          // last, calculate the percentage of material that has undergone the transition
          return (pressure_deviation > 0) ? 1 : 0;
        }

      // if we do not have the adiabatic conditions, we have to use the depth instead
      // this is less precise, because we do not have the exact pressure gradient, instead we use pressure/depth
      // (this is for calculating e.g. the density in the adiabatic profile)
      else
        {
          const double depth = this->get_geometry_model().depth(position);
          const double depth_deviation = (pressure > 0
                                          ?
                                          depth - transition_depths[phase]
                                          - transition_slopes[phase] * (depth / pressure) * (temperature - transition_temperatures[phase])
                                          :
                                          depth - transition_depths[phase]
                                          - transition_slopes[phase] / (this->get_gravity_model().gravity_vector(position).norm() * reference_rho)
                                          * (temperature - transition_temperatures[phase]));

          return (depth_deviation > 0) ? 1 : 0;
        }
    }



    template <int dim>
    unsigned int
    EquilibriumGrainSize<dim>::
    get_phase_index (const Point<dim> &position,
                     const double temperature,
                     const double pressure) const
    {
      Assert(grain_growth_activation_energy.size()>0,
             ExcMessage("Error: No grain evolution parameters are given!"));

      unsigned int phase_index = 0;
      if (transition_depths.size()>0)
        if (phase_function(position, temperature, pressure, transition_depths.size()-1) == 1)
          phase_index = transition_depths.size();

      for (unsigned int j=1; j<transition_depths.size(); ++j)
        if (phase_function(position, temperature, pressure, j) != phase_function(position, temperature, pressure, j-1))
          phase_index = j;

      return phase_index;
    }



    template <int dim>
    double
    EquilibriumGrainSize<dim>::
    diffusion_viscosity (const double                  temperature,
                         const double                  pressure,
                         const std::vector<double>    &composition,
                         const SymmetricTensor<2,dim> &strain_rate,
                         const Point<dim>             &position) const
    {
      const SymmetricTensor<2,dim> shear_strain_rate = strain_rate - 1./dim * trace(strain_rate) * unit_symmetric_tensor<dim>();
      const double second_strain_rate_invariant = std::sqrt(std::abs(second_invariant(shear_strain_rate)));

      const double grain_size = composition[this->introspection().compositional_index_for_name("grain_size")];

      // Currently this will never be called without adiabatic_conditions initialized, but just in case
      const double adiabatic_pressure = this->get_adiabatic_conditions().is_initialized()
                                        ?
                                        this->get_adiabatic_conditions().pressure(position)
                                        :
                                        pressure;

      // find out in which phase we are
      const unsigned int phase_index = get_phase_index(position, temperature, adiabatic_pressure);

      Assert(temperature > 0.0,
             ExcMessage("Temperature has to be positive for diffusion creep. Instead it is: " + std::to_string(temperature)));
      Assert(pressure >= 0.0,
             ExcMessage("Adiabatic pressure has to be non-negative for diffusion creep. Instead it is: " + std::to_string(adiabatic_pressure)));

      double energy_term = std::exp((diffusion_activation_energy[phase_index] + diffusion_activation_volume[phase_index] * adiabatic_pressure)
                                    / (diffusion_creep_exponent[phase_index] * constants::gas_constant * temperature));

      Assert(energy_term > 0.0,
             ExcMessage("Energy term has to be positive for diffusion creep. Instead it is: " + std::to_string(energy_term)));

      // If the adiabatic profile is already calculated we can use it to limit
      // variations in viscosity due to temperature.
      if (this->get_adiabatic_conditions().is_initialized())
        {
          const double adiabatic_temperature = this->get_adiabatic_conditions().temperature(position);
          Assert(adiabatic_temperature > 0.0,
                 ExcMessage("Adiabatic temperature has to be positive for diffusion creep. Instead it is: " + std::to_string(adiabatic_temperature)));

          const double adiabatic_energy_term
            = std::exp((diffusion_activation_energy[phase_index] + diffusion_activation_volume[phase_index] * adiabatic_pressure)
                       / (diffusion_creep_exponent[phase_index] * constants::gas_constant * adiabatic_temperature));

          Assert(adiabatic_energy_term > 0.0,
                 ExcMessage("Adiabatic energy term has to be positive for diffusion creep. Instead it is: " + std::to_string(adiabatic_energy_term)));

          const double temperature_dependence = energy_term / adiabatic_energy_term;
          if (temperature_dependence > max_temperature_dependence_of_eta)
            energy_term = adiabatic_energy_term * max_temperature_dependence_of_eta;
          if (temperature_dependence < 1.0 / max_temperature_dependence_of_eta)
            energy_term = adiabatic_energy_term / max_temperature_dependence_of_eta;
        }

      const double strain_rate_dependence = (1.0 - diffusion_creep_exponent[phase_index]) / diffusion_creep_exponent[phase_index];

      const double diffusion_viscosity = std::pow(diffusion_creep_prefactor[phase_index],-1.0/diffusion_creep_exponent[phase_index])
                                         * std::pow(second_strain_rate_invariant,strain_rate_dependence)
                                         * std::pow(grain_size, diffusion_creep_grain_size_exponent[phase_index]/diffusion_creep_exponent[phase_index])
                                         * energy_term;

      Assert(diffusion_viscosity > 0.0,
             ExcMessage("Diffusion viscosity has to be positive. Instead it is: " + std::to_string(diffusion_viscosity)));

      return diffusion_viscosity;
    }



    template <int dim>
    double
    EquilibriumGrainSize<dim>::
    dislocation_viscosity (const double      temperature,
                           const double      pressure,
                           const std::vector<double> &composition,
                           const SymmetricTensor<2,dim> &strain_rate,
                           const Point<dim> &position,
                           const double viscosity_guess) const
    {
      const double diff_viscosity = diffusion_viscosity(temperature,pressure,composition,strain_rate,position);

      // If we do not have the strain rate yet, set it to a low value.
      // TODO: make minimum strain rate an input parameter
      const SymmetricTensor<2,dim> shear_strain_rate = strain_rate - 1./dim * trace(strain_rate) * unit_symmetric_tensor<dim>();
      const double second_strain_rate_invariant = std::max(std::sqrt(std::abs(second_invariant(shear_strain_rate))),1e-30);

      // Start the iteration with the full strain rate
      double dis_viscosity = viscosity_guess;
      if (viscosity_guess == 0)
        dis_viscosity = dislocation_viscosity_fixed_strain_rate(temperature,pressure,second_strain_rate_invariant,position);

      double dis_viscosity_old = 0;
      unsigned int i = 0;
      while ((std::abs((dis_viscosity-dis_viscosity_old) / dis_viscosity) > dislocation_viscosity_iteration_threshold)
             && (i < dislocation_viscosity_iteration_number))
        {
          const double dislocation_strain_rate = diff_viscosity / (diff_viscosity + dis_viscosity) * second_strain_rate_invariant;
          dis_viscosity_old = dis_viscosity;
          dis_viscosity = dislocation_viscosity_fixed_strain_rate(temperature,
                                                                  pressure,
                                                                  dislocation_strain_rate,
                                                                  position);
          ++i;

          Assert(dis_viscosity > 0.0,
                 ExcMessage("Encountered negative dislocation viscosity in iteration " + std::to_string(i) +
                            ". Dislocation viscosity is: " + std::to_string(dis_viscosity)));
        }

      Assert(i<dislocation_viscosity_iteration_number,ExcInternalError());

      return dis_viscosity;
    }



    template <int dim>
    double
    EquilibriumGrainSize<dim>::
    dislocation_viscosity_fixed_strain_rate (const double      temperature,
                                             const double      pressure,
                                             const double      second_strain_rate_invariant,
                                             const Point<dim> &position) const
    {
      Assert(temperature > 0.0,
             ExcMessage("Temperature has to be positive for dislocation creep. Instead it is: " + std::to_string(temperature)));
      Assert(pressure >= 0.0,
             ExcMessage("Pressure has to be non-negative for dislocation creep. Instead it is: " + std::to_string(pressure)));

      // find out in which phase we are
      const unsigned int phase_index = get_phase_index(position, temperature, pressure);

      double energy_term = std::exp((dislocation_activation_energy[phase_index] + dislocation_activation_volume[phase_index] * pressure)
                                    / (dislocation_creep_exponent[phase_index] * constants::gas_constant * temperature));

      Assert(energy_term > 0.0,
             ExcMessage("Energy term has to be positive for dislocation creep. Instead it is: " + std::to_string(energy_term)));

      // If we are past the initialization of the adiabatic profile, use it to
      // limit viscosity variations due to temperature.
      if (this->get_adiabatic_conditions().is_initialized())
        {
          const double adiabatic_energy_term
            = exp((dislocation_activation_energy[phase_index] + dislocation_activation_volume[phase_index] * pressure)
                  / (dislocation_creep_exponent[phase_index] * constants::gas_constant * this->get_adiabatic_conditions().temperature(position)));

          Assert(adiabatic_energy_term > 0.0,
                 ExcMessage("Adiabatic energy term has to be positive for dislocation creep. Instead it is: " + std::to_string(adiabatic_energy_term)));

          const double temperature_dependence = energy_term / adiabatic_energy_term;
          if (temperature_dependence > max_temperature_dependence_of_eta)
            energy_term = adiabatic_energy_term * max_temperature_dependence_of_eta;
          if (temperature_dependence < 1.0 / max_temperature_dependence_of_eta)
            energy_term = adiabatic_energy_term / max_temperature_dependence_of_eta;
        }

      const double strain_rate_dependence = (1.0 - dislocation_creep_exponent[phase_index]) / dislocation_creep_exponent[phase_index];

      const double dislocation_viscosity = std::pow(dislocation_creep_prefactor[phase_index],-1.0/dislocation_creep_exponent[phase_index])
                                           * std::pow(second_strain_rate_invariant,strain_rate_dependence)
                                           * energy_term;

      Assert(dislocation_viscosity > 0.0,
             ExcMessage("Dislocation viscosity has to be positive. Instead it is: " + std::to_string(dislocation_viscosity)));

      return dislocation_viscosity;
    }



    template <int dim>
    double
    EquilibriumGrainSize<dim>::
    seismic_Vp (const double      temperature,
                const double      pressure,
                const std::vector<double> &compositional_fields,
                const Point<dim> &/*position*/) const
    {
      AssertThrow ((reference_compressibility != 0.0) || use_table_properties,
                   ExcMessage("Currently only compressible models are supported."));

      double vp = 0.0;
      if (n_material_data == 1)
        vp = material_lookup[0]->seismic_Vp(temperature,pressure);
      else
        {
          for (unsigned i = 0; i < n_material_data; i++)
            vp += compositional_fields[i] * material_lookup[i]->seismic_Vp(temperature,pressure);
        }
      return vp;
    }



    template <int dim>
    double
    EquilibriumGrainSize<dim>::
    seismic_Vs (const double      temperature,
                const double      pressure,
                const std::vector<double> &compositional_fields,
                const Point<dim> &/*position*/) const
    {
      AssertThrow ((reference_compressibility != 0.0) || use_table_properties,
                   ExcMessage("Currently only compressible models are supported."));

      double vs = 0.0;
      if (n_material_data == 1)
        vs = material_lookup[0]->seismic_Vs(temperature,pressure);
      else
        {
          for (unsigned i = 0; i < n_material_data; i++)
            vs += compositional_fields[i] * material_lookup[i]->seismic_Vs(temperature,pressure);
        }
      return vs;
    }



    template <int dim>
    double
    EquilibriumGrainSize<dim>::
    reference_viscosity () const
    {
      return eta;
    }



    template <int dim>
    double
    EquilibriumGrainSize<dim>::
    density (const double temperature,
             const double pressure,
             const std::vector<double> &compositional_fields, /*composition*/
             const Point<dim> &) const
    {
      if (!use_table_properties)
        {
          return reference_rho * std::exp(reference_compressibility * (pressure - this->get_surface_pressure()))
                 * (1 - thermal_alpha * (temperature - reference_T));
        }
      else
        {
          double rho = 0.0;
          if (n_material_data == 1)
            {
              rho = material_lookup[0]->density(temperature,pressure);
            }
          else
            {
              for (unsigned i = 0; i < n_material_data; i++)
                rho += compositional_fields[i] * material_lookup[i]->density(temperature,pressure);
            }

          return rho;
        }
    }



    template <int dim>
    bool
    EquilibriumGrainSize<dim>::
    is_compressible () const
    {
      return (reference_compressibility != 0)
             || use_table_properties;
    }



    template <int dim>
    double
    EquilibriumGrainSize<dim>::
    compressibility (const double temperature,
                     const double pressure,
                     const std::vector<double> &compositional_fields,
                     const Point<dim> &position) const
    {
      if (!use_table_properties)
        return reference_compressibility;

      double dRhodp = 0.0;
      if (n_material_data == 1)
        dRhodp = material_lookup[0]->dRhodp(temperature,pressure);
      else
        {
          for (unsigned i = 0; i < n_material_data; i++)
            dRhodp += compositional_fields[i] * material_lookup[i]->dRhodp(temperature,pressure);
        }
      const double rho = density(temperature,pressure,compositional_fields,position);
      return (1/rho)*dRhodp;
    }



    template <int dim>
    double
    EquilibriumGrainSize<dim>::
    thermal_expansivity (const double temperature,
                         const double pressure,
                         const std::vector<double> &compositional_fields,
                         const Point<dim> &) const
    {
      if (!use_table_properties)
        return thermal_alpha;

      double alpha = 0.0;
      if (n_material_data == 1)
        alpha = material_lookup[0]->thermal_expansivity(temperature,pressure);
      else
        {
          for (unsigned i = 0; i < n_material_data; i++)
            alpha += compositional_fields[i] * material_lookup[i]->thermal_expansivity(temperature,pressure);
        }
      return alpha;
    }



    template <int dim>
    void
    EquilibriumGrainSize<dim>::
    evaluate(const typename Interface<dim>::MaterialModelInputs &in, typename Interface<dim>::MaterialModelOutputs &out) const
    {
      // Determine some properties that are constant for all points
      const unsigned int vs_anomaly_index = this->introspection().compositional_index_for_name("vs_anomaly");
      const unsigned int surface_boundary_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id("outer");
      const unsigned int craton_composition = this->introspection().compositional_index_for_name("continents");

      const InitialTemperature::AdiabaticBoundary<dim> &adiabatic_boundary =
        this->get_initial_temperature_manager().template get_matching_initial_temperature_model<InitialTemperature::AdiabaticBoundary<dim> >();

      // This function will fill the outputs for grain size, viscosity, and dislocation viscosity
      if (in.requests_property(MaterialProperties::viscosity))
        compute_equilibrium_grain_size(in, out);

      for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
        {
          // Use the adiabatic pressure instead of the real one, because of oscillations
          const double pressure = (this->get_adiabatic_conditions().is_initialized())
                                  ?
                                  this->get_adiabatic_conditions().pressure(in.position[i])
                                  :
                                  in.pressure[i];

          out.thermal_conductivities[i] = k_value;
          out.specific_heat[i] = reference_specific_heat;

          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c] = 0.0;

          const double depth = this->get_geometry_model().depth(in.position[i]);
          const double delta_log_vs = in.composition[i][vs_anomaly_index];

          // For computing the adiabatic conditions, we want to use the temperature it hands
          // over as input. Once the adiabatic conditions are computed, we want to use the
          // temperature based on seismic tomography and the input temperature model.
          double new_temperature = in.temperature[i];

          // For the adiabatic conditions, assume a characteristic crustal and lithosphere thickness.
          double crustal_thickness = 40000.;
          double lithosphere_thickness = 100000.;

          if (this->get_adiabatic_conditions().is_initialized())
            {
              // Get variable lithosphere and crustal depths using an adiabatic boundary ascii file
              lithosphere_thickness = adiabatic_boundary.get_data_component(surface_boundary_id, in.position[i], 0);
              crustal_thickness = crustal_boundary_depth.get_data_component(surface_boundary_id, in.position[i], 0);

              // This variable stores the seismic tomography based temperatures
              double mantle_temperature;

              // This does not work when it is called before the adiabatic conditions are initialized, because
              // the AsciiDataBoundary plugin needs the time, which is not yet initialized at that point.
              const double initial_temperature = this->get_initial_temperature_manager().initial_temperature(in.position[i]);
              const double reference_temperature = this->get_adiabatic_conditions().temperature(in.position[i]);

              if (use_depth_dependent_dT_vs)
                {
                  // We use the dlnvs/dT profile from Steinberger and Calderwood (2006). The values in profile are in units 1/K .
                  mantle_temperature = reference_temperature +
                                       delta_log_vs * dT_vs_depth_profile.get_data_component(Point<1>(depth), temperature_scaling_index);
                }

              else
                {
                  // compute the temperature below the asthenosphere using the parameters given in the table by Becker (2006).
                  mantle_temperature = reference_temperature + delta_log_vs * -4.2 * 1785.;
                }

              const double sigmoid_width = 2.e4;
              const double sigmoid = 1.0 / (1.0 + std::exp( (uppermost_mantle_thickness - depth)/sigmoid_width));

              new_temperature = initial_temperature + (mantle_temperature - initial_temperature) * sigmoid;
            }

          if (PrescribedTemperatureOutputs<dim> *prescribed_temperature_out = out.template get_additional_output<PrescribedTemperatureOutputs<dim> >())
            prescribed_temperature_out->prescribed_temperature_outputs[i] = new_temperature;

          if (use_table_properties)
            {
              // fill seismic velocities outputs if they exist
              if (SeismicAdditionalOutputs<dim> *seismic_out = out.template get_additional_output<SeismicAdditionalOutputs<dim> >())
                {
                  seismic_out->vp[i] = seismic_Vp(in.temperature[i], in.pressure[i], in.composition[i], in.position[i]);
                  seismic_out->vs[i] = seismic_Vs(in.temperature[i], in.pressure[i], in.composition[i], in.position[i]);
                }

              out.densities[i] = density(new_temperature, pressure, in.composition[i], in.position[i]);
              out.thermal_expansion_coefficients[i] = thermal_expansivity(new_temperature, pressure, in.composition[i], in.position[i]);
              out.compressibilities[i] = compressibility(new_temperature, pressure, in.composition[i], in.position[i]);
            }
          else
            {
              // Temperature and density of the upper part of the mantle are computed separately,
              // based on Tutu et al. (2018).
              // Reference temperature is 20 C in Tutu et al. (2018).
              double deltaT  = new_temperature - 293.;

              unsigned int material_type = 0;

              // Density computation
              if (depth <= crustal_thickness)
                {
                  out.thermal_expansion_coefficients[i] = 2.7e-5;
                  out.compressibilities[i] = 1./6.3e10;
                  out.densities[i] = 2.85e3 * (1. - out.thermal_expansion_coefficients[i] * deltaT
                                               + pressure * out.compressibilities[i]);
                  material_type = 1;
                }
              else if (depth > crustal_thickness && depth <= lithosphere_thickness)
                {
                  out.thermal_expansion_coefficients[i] = 3.e-5;
                  out.compressibilities[i] = 1./12.2e10;
                  out.densities[i] = 3.27e3 * (1. - out.thermal_expansion_coefficients[i] * deltaT
                                               + pressure * out.compressibilities[i]);
                  material_type = 2;
                }
              else if (depth > lithosphere_thickness && depth <= uppermost_mantle_thickness)
                {
                  out.thermal_expansion_coefficients[i] = 3.e-5;
                  out.compressibilities[i] = 1./12.2e10;
                  out.densities[i] = 3.3e3 * (1. - out.thermal_expansion_coefficients[i] * deltaT
                                              + pressure * out.compressibilities[i]);
                  material_type = 3;
                }
              else
                {
                  // Densities below 300 km are computed using the scaling relationship from the velocity anomalies
                  double density_vs_scaling;
                  if (use_depth_dependent_rho_vs)
                    density_vs_scaling = rho_vs_depth_profile.get_data_component(Point<1>(depth), density_scaling_index);
                  else
                    // Values from Becker [2006], GJI
                    density_vs_scaling = 0.15;

                  const double density_anomaly = delta_log_vs * density_vs_scaling;
                  const double reference_density = this->get_adiabatic_conditions().is_initialized()
                                                   ?
                                                   this->get_adiabatic_conditions().density(in.position[i])
                                                   :
                                                   reference_rho;

                  out.densities[i] = reference_density * (1. + density_anomaly);

                  if (use_depth_dependent_thermal_expansivity)
                    out.thermal_expansion_coefficients[i] = thermal_expansivity_profile.get_data_component(Point<1>(depth), thermal_expansivity_column_index);
                  else
                    out.thermal_expansion_coefficients[i] = thermal_alpha;

                  out.compressibilities[i] = reference_compressibility;

                  material_type = 4;
                }

              if (MaterialTypeAdditionalOutputs<dim> *material_type_out = out.template get_additional_output<MaterialTypeAdditionalOutputs<dim> >())
                material_type_out->output_values[0][i] = material_type;
            }
        }
    }



    template <int dim>
    void
    EquilibriumGrainSize<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Equilibrium grain size model");
        {
          prm.declare_entry ("Reference density", "3300",
                             Patterns::Double (0),
                             "The reference density $\\rho_0$. Units: $kg/m^3$.");
          prm.declare_entry ("Reference temperature", "293",
                             Patterns::Double (0),
                             "The reference temperature $T_0$. Units: $K$.");
          prm.declare_entry ("Viscosity", "5e24",
                             Patterns::Double (0),
                             "The value of the constant viscosity. Units: $kg/m/s$.");
          prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0),
                             "The value of the thermal conductivity $k$. "
                             "Units: $W/m/K$.");
          prm.declare_entry ("Reference specific heat", "1250",
                             Patterns::Double (0),
                             "The value of the specific heat $cp$. "
                             "Units: $J/kg/K$.");
          prm.declare_entry ("Thermal expansion coefficient", "2e-5",
                             Patterns::Double (0),
                             "The value of the thermal expansion coefficient $\\alpha$. "
                             "Units: $1/K$.");
          prm.declare_entry ("Reference compressibility", "4e-12",
                             Patterns::Double (0),
                             "The value of the reference compressibility. "
                             "Units: $1/Pa$.");
          prm.declare_entry ("Phase transition depths", "",
                             Patterns::List (Patterns::Double(0)),
                             "A list of depths where phase transitions occur. Values must "
                             "monotonically increase. "
                             "Units: $m$.");
          prm.declare_entry ("Phase transition temperatures", "",
                             Patterns::List (Patterns::Double(0)),
                             "A list of temperatures where phase transitions occur. Higher or lower "
                             "temperatures lead to phase transition ocurring in smaller or greater "
                             "depths than given in Phase transition depths, depending on the "
                             "Clapeyron slope given in Phase transition Clapeyron slopes. "
                             "List must have the same number of entries as Phase transition depths. "
                             "Units: $K$.");
          prm.declare_entry ("Phase transition widths", "",
                             Patterns::List (Patterns::Double(0)),
                             "A list of widths for each phase transition. This is only use to specify "
                             "the region where the recrystallized grain size is assigned after material "
                             "has crossed a phase transition and should accordingly be chosen similar "
                             "to the maximum cell width expected at the phase transition."
                             "List must have the same number of entries as Phase transition depths. "
                             "Units: $m$.");
          prm.declare_entry ("Phase transition Clapeyron slopes", "",
                             Patterns::List (Patterns::Double()),
                             "A list of Clapeyron slopes for each phase transition. A positive "
                             "Clapeyron slope indicates that the phase transition will occur in "
                             "a greater depth, if the temperature is higher than the one given in "
                             "Phase transition temperatures and in a smaller depth, if the "
                             "temperature is smaller than the one given in Phase transition temperatures. "
                             "For negative slopes the other way round. "
                             "List must have the same number of entries as Phase transition depths. "
                             "Units: $Pa/K$.");
          prm.declare_entry ("Grain growth activation energy", "3.5e5",
                             Patterns::List (Patterns::Double(0)),
                             "The activation energy for grain growth $E_g$. "
                             "Units: $J/mol$.");
          prm.declare_entry ("Grain growth activation volume", "8e-6",
                             Patterns::List (Patterns::Double(0)),
                             "The activation volume for grain growth $V_g$. "
                             "Units: $m^3/mol$.");
          prm.declare_entry ("Grain growth exponent", "3",
                             Patterns::List (Patterns::Double(0)),
                             "The exponent of the grain growth law $p_g$. This is an experimentally determined "
                             "grain growth constant. "
                             "Units: none.");
          prm.declare_entry ("Grain growth rate constant", "1.5e-5",
                             Patterns::List (Patterns::Double(0)),
                             "The prefactor for the Ostwald ripening grain growth law $G_0$. "
                             "This is dependent on water content, which is assumed to be "
                             "50 H/$10^6$ Si for the default value. "
                             "Units: $m^{p_g}/s$.");
          prm.declare_entry ("Minimum grain size", "5e-6",
                             Patterns::Double(0),
                             "The minimum allowable grain size. The grain size will be limited to be "
                             "larger than this value. This can be used to damp out oscillations, or "
                             "to limit the viscosity variation due to grain size. "
                             "Units: $m$.");
          prm.declare_entry ("Reciprocal required strain", "10",
                             Patterns::List (Patterns::Double(0)),
                             "This parameter ($\\lambda$) gives an estimate of the strain necessary "
                             "to achieve a new grain size. ");
          prm.declare_entry ("Recrystallized grain size", "",
                             Patterns::List (Patterns::Double(0)),
                             "The grain size $d_{ph}$ to that a phase will be reduced to when crossing a phase transition. "
                             "When set to zero, grain size will not be reduced. "
                             "Units: m.");
          prm.declare_entry ("Use equilibrium grain size", "true",
                             Patterns::Bool (),
                             "A flag indicating whether the computation should use the equilibrium grain size "
                             "when computing the viscosity (if true), or use a constant grain size instead (if "
                             "false).");
          prm.declare_entry ("Use paleowattmeter", "true",
                             Patterns::Bool (),
                             "A flag indicating whether the computation should use the "
                             "paleowattmeter approach of Austin and Evans (2007) for grain size reduction "
                             "in the dislocation creep regime (if true) or the paleopiezometer approach "
                             "from Hall and Parmetier (2003) (if false).");
          prm.declare_entry ("Average specific grain boundary energy", "1.0",
                             Patterns::List (Patterns::Double(0)),
                             "The average specific grain boundary energy $\\gamma$. "
                             "Units: $J/m^2$.");
          prm.declare_entry ("Work fraction for boundary area change", "0.1",
                             Patterns::List (Patterns::Double(0)),
                             "The fraction $\\chi$ of work done by dislocation creep to change the grain boundary area. "
                             "Units: $J/m^2$.");
          prm.declare_entry ("Geometric constant", "3",
                             Patterns::List (Patterns::Double(0)),
                             "The geometric constant $c$ used in the paleowattmeter grain size reduction law. "
                             "Units: none.");
          prm.declare_entry ("Dislocation viscosity iteration threshold", "1e-3",
                             Patterns::Double(0),
                             "We need to perform an iteration inside the computation "
                             "of the dislocation viscosity, because it depends on the "
                             "dislocation strain rate, which depends on the dislocation "
                             "viscosity itself. This number determines the termination "
                             "accuracy, i.e. if the dislocation viscosity changes by less "
                             "than this factor we terminate the iteration.");
          prm.declare_entry ("Dislocation viscosity iteration number", "100",
                             Patterns::Integer(0),
                             "We need to perform an iteration inside the computation "
                             "of the dislocation viscosity, because it depends on the "
                             "dislocation strain rate, which depends on the dislocation "
                             "viscosity itself. This number determines the maximum "
                             "number of iterations that are performed. ");
          prm.declare_entry ("Dislocation creep exponent", "3.5",
                             Patterns::List (Patterns::Double(0)),
                             "The power-law exponent $n_{dis}$ for dislocation creep. "
                             "Units: none.");
          prm.declare_entry ("Dislocation activation energy", "4.8e5",
                             Patterns::List (Patterns::Double(0)),
                             "The activation energy for dislocation creep $E_{dis}$. "
                             "Units: $J/mol$.");
          prm.declare_entry ("Dislocation activation volume", "1.1e-5",
                             Patterns::List (Patterns::Double(0)),
                             "The activation volume for dislocation creep $V_{dis}$. "
                             "Units: $m^3/mol$.");
          prm.declare_entry ("Dislocation creep prefactor", "4.5e-15",
                             Patterns::List (Patterns::Double(0)),
                             "The prefactor for the dislocation creep law $A_{dis}$. "
                             "Units: $Pa^{-n_{dis}}/s$.");
          prm.declare_entry ("Diffusion creep exponent", "1",
                             Patterns::List (Patterns::Double(0)),
                             "The power-law exponent $n_{diff}$ for diffusion creep. "
                             "Units: none.");
          prm.declare_entry ("Diffusion activation energy", "3.35e5",
                             Patterns::List (Patterns::Double(0)),
                             "The activation energy for diffusion creep $E_{diff}$. "
                             "Units: $J/mol$.");
          prm.declare_entry ("Diffusion activation volume", "4e-6",
                             Patterns::List (Patterns::Double(0)),
                             "The activation volume for diffusion creep $V_{diff}$. "
                             "Units: $m^3/mol$.");
          prm.declare_entry ("Diffusion creep prefactor", "7.4e-15",
                             Patterns::List (Patterns::Double(0)),
                             "The prefactor for the diffusion creep law $A_{diff}$. "
                             "Units: $m^{p_{diff}} Pa^{-n_{diff}}/s$.");
          prm.declare_entry ("Diffusion creep grain size exponent", "3",
                             Patterns::List (Patterns::Double(0)),
                             "The diffusion creep grain size exponent $p_{diff}$ that determines the "
                             "dependence of viscosity on grain size. "
                             "Units: none.");
          prm.declare_entry ("Maximum temperature dependence of viscosity", "100",
                             Patterns::Double (1.0),
                             "The factor by which viscosity at adiabatic temperature and ambient temperature "
                             "are allowed to differ (a value of x means that the viscosity can be x times higher "
                             "or x times lower compared to the value at adiabatic temperature. This parameter "
                             "is introduced to limit local viscosity contrasts, but still allow for a widely "
                             "varying viscosity over the whole mantle range. "
                             "Units: none.");
          prm.declare_entry ("Minimum viscosity", "1e18",
                             Patterns::Double (0),
                             "The minimum viscosity that is allowed in the whole model domain. "
                             "Units: Pa \\, s.");
          prm.declare_entry ("Maximum viscosity", "1e26",
                             Patterns::Double (0),
                             "The maximum viscosity that is allowed in the whole model domain. "
                             "Units: Pa \\, s.");
          prm.declare_entry ("Minimum specific heat", "500",
                             Patterns::Double (0),
                             "The minimum specific heat that is allowed in the whole model domain. "
                             "Units: J/kg/K.");
          prm.declare_entry ("Maximum specific heat", "6000",
                             Patterns::Double (0),
                             "The maximum specific heat that is allowed in the whole model domain. "
                             "Units: J/kg/K.");
          prm.declare_entry ("Minimum thermal expansivity", "1e-5",
                             Patterns::Double (),
                             "The minimum thermal expansivity that is allowed in the whole model domain. "
                             "Units: 1/K.");
          prm.declare_entry ("Maximum thermal expansivity", "1e-3",
                             Patterns::Double (),
                             "The maximum thermal expansivity that is allowed in the whole model domain. "
                             "Units: 1/K.");
          prm.declare_entry ("Minimum grain size", "1e-5",
                             Patterns::Double (0),
                             "The minimum grain size that is used for the material model. This parameter "
                             "is introduced to limit local viscosity contrasts, but still allows for a widely "
                             "varying viscosity over the whole mantle range. "
                             "Units: m.");
          prm.declare_entry ("Data directory", "$ASPECT_SOURCE_DIR/data/material-model/steinberger/",
                             Patterns::DirectoryName (),
                             "The path to the model data. The path may also include the special "
                             "text '$ASPECT_SOURCE_DIR' which will be interpreted as the path "
                             "in which the ASPECT source files were located when ASPECT was "
                             "compiled. This interpretation allows, for example, to reference "
                             "files located in the 'data/' subdirectory of ASPECT. ");
          prm.declare_entry ("Material file names", "pyr-ringwood88.txt",
                             Patterns::List (Patterns::Anything()),
                             "The file names of the material data. "
                             "List with as many components as active "
                             "compositional fields (material data is assumed to "
                             "be in order with the ordering of the fields). ");
          prm.declare_entry ("Derivatives file names", "",
                             Patterns::List (Patterns::Anything()),
                             "The file names of the enthalpy derivatives data. "
                             "List with as many components as active "
                             "compositional fields (material data is assumed to "
                             "be in order with the ordering of the fields). ");
          prm.declare_entry ("Use table properties", "false",
                             Patterns::Bool(),
                             "This parameter determines whether to use the table properties "
                             "also for density, thermal expansivity and specific heat. "
                             "If false the properties are generated as in the "
                             "simple compressible plugin.");
          prm.declare_entry ("Material file format", "perplex",
                             Patterns::Selection ("perplex|hefesto"),
                             "The material file format to be read in the property "
                             "tables.");
          prm.declare_entry ("Use depth dependent viscosity", "false",
                             Patterns::Bool (),
                             "This parameter value determines if we want to use the layered depth dependent "
                             "rheology, which is input as an ascii data file.");
          prm.declare_entry ("Use faults", "false",
                             Patterns::Bool (),
                             "This parameter value determines if we want to use the faults/plate boundaries as "
                             "a composition field, currently input as a world builder file.");
          prm.declare_entry ("Fault viscosity", "1e20",
                             Patterns::Double(0),
                             "This parameter value determines the viscosity of faults or plate boundaries. "
                             "We would want to have weak faults/plate boundaries relative to the surrounding "
                             "lithosphere. "
                             "Units: Pa.s");
          prm.declare_entry ("Asthenosphere viscosity", "2.4e20",
                             Patterns::Double(0),
                             "This parameter value determines the asthenosphere layer in the reference file. "
                             "Units: Pa.s");
          prm.declare_entry ("Use depth dependent density scaling", "false",
                             Patterns::Bool (),
                             "This parameter value determines if we want to use depth-dependent scaling files "
                             "for computing the density from seismic velocities (if true) or use a constant "
                             "value (if false).");
          prm.declare_entry ("Use thermal expansivity profile", "true",
                             Patterns::Bool (),
                             "This parameter determines if we use a depth-dependent thermal expansivity read "
                             "from a data file (if true), or if we use the constant reference value given by "
                             "the 'Thermal expansion coefficient' (if false). In the case that material "
                             "properties are read from a look-up table, as determined by the input parameter "
                             "'Use table properties', this value is irrelevant and will be ignored.");
          prm.declare_entry ("Uppermost mantle thickness", "300000",
                             Patterns::Double (0),
                             "The depth of the base of the uppoermost mantle, which marks the transition between "
                             "using the temperature model of Tutu et al. (above) and derived from seismic "
                             "tomography (below). "
                             "Units: m.");
          prm.declare_entry ("Use depth dependent temperature scaling", "false",
                             Patterns::Bool (),
                             "This parameter value determines if we want to use a depth-dependent scaling read "
                             "in from a data file for computing the temperature from seismic velocities (if true) "
                             "or use a constant value (if false).");

          // Depth-dependent viscosity parameters
          Rheology::AsciiDepthProfile<dim>::declare_parameters(prm);

          // Depth-dependent density scaling parameters
          Utilities::AsciiDataBase<dim>::declare_parameters(prm, "../../input_data/", "rho_vs_scaling.txt", "Density velocity scaling");

          // Depth-dependent density scaling parameters
          Utilities::AsciiDataBase<dim>::declare_parameters(prm, "../../input_data/", "dT_vs_scaling.txt", "Temperature velocity scaling");

          // Depth-dependent thermal expansivity parameters
          Utilities::AsciiDataBase<dim>::declare_parameters(prm, "../../input_data/", "thermal_expansivity_steinberger_calderwood.txt", "Thermal expansivity profile");

          // Crustal boundary depths parameters
          Utilities::AsciiDataBoundary<dim>::declare_parameters(prm,  "../../input_data/crust1.0/", "crustal_structure.txt", "Crustal depths");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    EquilibriumGrainSize<dim>::parse_parameters (ParameterHandler &prm)
    {
      AssertThrow (this->introspection().compositional_name_exists("grain_size"),
                   ExcMessage("The 'grain size' material model only works if a compositional "
                              "field with name 'grain_size' is present. Please use another material "
                              "model or add such a field."));

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Equilibrium grain size model");
        {
          reference_rho              = prm.get_double ("Reference density");
          reference_T                = prm.get_double ("Reference temperature");
          eta                        = prm.get_double ("Viscosity");
          k_value                    = prm.get_double ("Thermal conductivity");
          reference_specific_heat    = prm.get_double ("Reference specific heat");
          thermal_alpha              = prm.get_double ("Thermal expansion coefficient");
          reference_compressibility  = prm.get_double ("Reference compressibility");


          transition_depths         = Utilities::string_to_double
                                      (Utilities::split_string_list(prm.get ("Phase transition depths")));
          transition_temperatures   = Utilities::string_to_double
                                      (Utilities::split_string_list(prm.get ("Phase transition temperatures")));
          transition_slopes         = Utilities::string_to_double
                                      (Utilities::split_string_list(prm.get ("Phase transition Clapeyron slopes")));
          recrystallized_grain_size = Utilities::string_to_double
                                      (Utilities::split_string_list(prm.get ("Recrystallized grain size")));
          transition_widths         = Utilities::string_to_double
                                      (Utilities::split_string_list(prm.get ("Phase transition widths")));

          if (transition_temperatures.size() != transition_depths.size() ||
              transition_slopes.size() != transition_depths.size() ||
              transition_widths.size() != transition_depths.size() ||
              recrystallized_grain_size.size() != transition_depths.size() )
            AssertThrow(false,
                        ExcMessage("Error: At least one list that gives input parameters for the phase transitions has the wrong size."));

          if (transition_depths.size()>1)
            for (unsigned int i=0; i<transition_depths.size()-2; ++i)
              AssertThrow(transition_depths[i]<transition_depths[i+1],
                          ExcMessage("Error: Phase transition depths have to be sorted in ascending order!"));

          // grain evolution parameters
          grain_growth_activation_energy        = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Grain growth activation energy")));
          grain_growth_activation_volume        = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Grain growth activation volume")));
          grain_growth_rate_constant            = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Grain growth rate constant")));
          grain_growth_exponent                 = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Grain growth exponent")));
          minimum_grain_size                    = prm.get_double("Minimum grain size");
          reciprocal_required_strain            = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Reciprocal required strain")));
          equilibrate_grain_size                = prm.get_bool ("Use equilibrium grain size");

          use_paleowattmeter                    = prm.get_bool ("Use paleowattmeter");
          grain_boundary_energy                 = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Average specific grain boundary energy")));
          boundary_area_change_work_fraction    = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Work fraction for boundary area change")));
          geometric_constant                    = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Geometric constant")));

          // rheology parameters
          dislocation_viscosity_iteration_threshold = prm.get_double("Dislocation viscosity iteration threshold");
          dislocation_viscosity_iteration_number = prm.get_integer("Dislocation viscosity iteration number");
          dislocation_creep_exponent            = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Dislocation creep exponent")));
          dislocation_activation_energy         = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Dislocation activation energy")));
          dislocation_activation_volume         = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Dislocation activation volume")));
          dislocation_creep_prefactor           = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Dislocation creep prefactor")));
          diffusion_creep_exponent              = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Diffusion creep exponent")));
          diffusion_activation_energy           = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Diffusion activation energy")));
          diffusion_activation_volume           = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Diffusion activation volume")));
          diffusion_creep_prefactor             = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Diffusion creep prefactor")));
          diffusion_creep_grain_size_exponent   = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Diffusion creep grain size exponent")));
          max_temperature_dependence_of_eta     = prm.get_double ("Maximum temperature dependence of viscosity");
          min_eta                               = prm.get_double ("Minimum viscosity");
          max_eta                               = prm.get_double ("Maximum viscosity");
          min_specific_heat                     = prm.get_double ("Minimum specific heat");
          max_specific_heat                     = prm.get_double ("Maximum specific heat");
          min_thermal_expansivity               = prm.get_double ("Minimum thermal expansivity");
          max_thermal_expansivity               = prm.get_double ("Maximum thermal expansivity");
          min_grain_size                        = prm.get_double ("Minimum grain size");

          if (grain_growth_activation_energy.size() != grain_growth_activation_volume.size() ||
              grain_growth_activation_energy.size() != grain_growth_rate_constant.size() ||
              grain_growth_activation_energy.size() != grain_growth_exponent.size() ||
              grain_growth_activation_energy.size() != dislocation_creep_exponent.size() ||
              grain_growth_activation_energy.size() != dislocation_activation_energy.size() ||
              grain_growth_activation_energy.size() != dislocation_activation_volume.size() ||
              grain_growth_activation_energy.size() != dislocation_creep_prefactor.size() ||
              grain_growth_activation_energy.size() != diffusion_creep_exponent.size() ||
              grain_growth_activation_energy.size() != diffusion_activation_energy.size() ||
              grain_growth_activation_energy.size() != diffusion_activation_volume.size() ||
              grain_growth_activation_energy.size() != diffusion_creep_prefactor.size() ||
              grain_growth_activation_energy.size() != diffusion_creep_grain_size_exponent.size() )
            AssertThrow(false,
                        ExcMessage("Error: The lists of grain size evolution and flow law parameters "
                                   "need to have the same length!"));

          if (use_paleowattmeter)
            {
              if (grain_growth_activation_energy.size() != grain_boundary_energy.size() ||
                  grain_growth_activation_energy.size() != boundary_area_change_work_fraction.size() ||
                  grain_growth_activation_energy.size() != geometric_constant.size() )
                AssertThrow(false,
                            ExcMessage("Error: One of the lists of grain size evolution parameters "
                                       "given for the paleowattmeter does not have the correct length!"));
            }
          else
            AssertThrow(grain_growth_activation_energy.size() == reciprocal_required_strain.size(),
                        ExcMessage("Error: The list of grain size evolution parameters in the "
                                   "paleopiezometer does not have the correct length!"));

          AssertThrow(grain_growth_activation_energy.size() == transition_depths.size()+1,
                      ExcMessage("Error: The lists of grain size evolution and flow law parameters need to "
                                 "have exactly one more entry than the number of phase transitions "
                                 "(which is defined by the length of the lists of phase transition depths, ...)!"));

          // parameters for reading in tables with material properties
          datadirectory        = prm.get ("Data directory");
          {
            const std::string subst_text = "$ASPECT_SOURCE_DIR";
            std::string::size_type position;
            while (position = datadirectory.find (subst_text),  position!=std::string::npos)
              datadirectory.replace (datadirectory.begin()+position,
                                     datadirectory.begin()+position+subst_text.size(),
                                     ASPECT_SOURCE_DIR);
          }
          material_file_names  = Utilities::split_string_list
                                 (prm.get ("Material file names"));
          derivatives_file_names = Utilities::split_string_list
                                   (prm.get ("Derivatives file names"));

          use_table_properties                    = prm.get_bool ("Use table properties");
          use_depth_dependent_viscosity           = prm.get_bool ("Use depth dependent viscosity");
          use_faults                              = prm.get_bool ("Use faults");
          use_depth_dependent_rho_vs              = prm.get_bool ("Use depth dependent density scaling");
          use_depth_dependent_dT_vs               = prm.get_bool ("Use depth dependent temperature scaling");
          use_depth_dependent_thermal_expansivity = prm.get_bool ("Use thermal expansivity profile");
          uppermost_mantle_thickness              = prm.get_double ("Uppermost mantle thickness");
          fault_viscosity                        = prm.get_double ("Fault viscosity");
          asthenosphere_viscosity                 = prm.get_double ("Asthenosphere viscosity");

          // Parse all depth-dependent parameters
          if (use_depth_dependent_viscosity)
            {
              reference_viscosity_profile = std::make_unique<Rheology::AsciiDepthProfile<dim>>();
              reference_viscosity_profile->initialize_simulator (this->get_simulator());
              reference_viscosity_profile->parse_parameters(prm);
              reference_viscosity_profile->initialize();
            }

          if (use_depth_dependent_rho_vs)
            rho_vs_depth_profile.parse_parameters(prm, "Density velocity scaling");

          if (use_depth_dependent_thermal_expansivity)
            thermal_expansivity_profile.parse_parameters(prm, "Thermal expansivity profile");

          if (use_depth_dependent_dT_vs)
            dT_vs_depth_profile.parse_parameters(prm, "Temperature velocity scaling");

          crustal_boundary_depth.initialize_simulator (this->get_simulator());
          crustal_boundary_depth.parse_parameters(prm, "Crustal depths");

          // Make sure the grain size field comes after all potential material
          // data fields. Otherwise our material model calculation uses the
          // wrong compositional fields.
          if (use_table_properties && material_file_names.size() > 1)
            {
              AssertThrow(this->introspection().compositional_index_for_name("grain_size") >= material_file_names.size(),
                          ExcMessage("The compositional fields indicating the major element composition need to be first in the "
                                     "list of compositional fields, but the grain size field seems to have a lower index than the number "
                                     "of provided data files. This is likely inconsistent. Please check the number of provided data "
                                     "files and the order of compositional fields."));
            }

          if (prm.get ("Material file format") == "perplex")
            material_file_format = perplex;
          else if (prm.get ("Material file format") == "hefesto")
            material_file_format = hefesto;
          else
            AssertThrow (false, ExcNotImplemented());
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();


      // Declare dependencies on solution variables
      this->model_dependence.thermal_conductivity = NonlinearDependence::none;

      this->model_dependence.viscosity = NonlinearDependence::temperature
                                         | NonlinearDependence::pressure
                                         | NonlinearDependence::strain_rate
                                         | NonlinearDependence::compositional_fields;

      this->model_dependence.density = NonlinearDependence::none;
      this->model_dependence.compressibility = NonlinearDependence::none;
      this->model_dependence.specific_heat = NonlinearDependence::none;

      if (use_table_properties)
        {
          this->model_dependence.density |= NonlinearDependence::temperature
                                            | NonlinearDependence::pressure
                                            | NonlinearDependence::compositional_fields;
          this->model_dependence.compressibility = NonlinearDependence::temperature
                                                   | NonlinearDependence::pressure
                                                   | NonlinearDependence::compositional_fields;
          this->model_dependence.specific_heat = NonlinearDependence::temperature
                                                 | NonlinearDependence::pressure
                                                 | NonlinearDependence::compositional_fields;
        }
      else
        {
          if (thermal_alpha != 0)
            this->model_dependence.density |=NonlinearDependence::temperature;
          if (reference_compressibility != 0)
            this->model_dependence.density |=NonlinearDependence::pressure;
        }
    }



    template <int dim>
    void
    EquilibriumGrainSize<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      // These properties are useful as output, but will also be used by the
      // heating model to reduce shear heating by the amount of work done to
      // reduce grain size.
      if (out.template get_additional_output<DislocationViscosityOutputs<dim> >() == nullptr)
        {
          const unsigned int n_points = out.viscosities.size();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::DislocationViscosityOutputs<dim>> (n_points));
        }

      // We need the prescribed field outputs to interpolate the grain size onto a compositional field.
      if (out.template get_additional_output<PrescribedFieldOutputs<dim> >() == nullptr)
        {
          const unsigned int n_points = out.viscosities.size();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::PrescribedFieldOutputs<dim>> (n_points, this->n_compositional_fields()));
        }

      // These properties are only output properties. But we should only create them if they are filled.
      if (use_table_properties
          && out.template get_additional_output<SeismicAdditionalOutputs<dim> >() == nullptr)
        {
          const unsigned int n_points = out.viscosities.size();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::SeismicAdditionalOutputs<dim>> (n_points));
        }

      if (this->get_parameters().temperature_method == Parameters<dim>::AdvectionFieldMethod::prescribed_field &&
          out.template get_additional_output<PrescribedTemperatureOutputs<dim> >() == nullptr)
        {
          const unsigned int n_points = out.viscosities.size();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::PrescribedTemperatureOutputs<dim>> (n_points));
        }

      // We need additional field outputs for the unscaled viscosity
      if (out.template get_additional_output<UnscaledViscosityAdditionalOutputs<dim>>() == nullptr)
        {
          const unsigned int n_points = out.viscosities.size();
          out.additional_outputs.push_back(
            std::make_unique<UnscaledViscosityAdditionalOutputs<dim>> (n_points));
        }

      if (out.template get_additional_output<MaterialTypeAdditionalOutputs<dim>>() == nullptr)
        {
          const unsigned int n_points = out.viscosities.size();
          out.additional_outputs.push_back(
            std::make_unique<MaterialTypeAdditionalOutputs<dim>> (n_points));
        }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(EquilibriumGrainSize,
                                   "equilibrium grain size",
                                   "A material model that relies on compositional "
                                   "fields that correspond to the average grain sizes of a "
                                   "mineral phase and source terms that determine the grain "
                                   "size evolution in terms of the strain rate, "
                                   "temperature, phase transitions, and the creep regime. "
                                   "This material model only works if a compositional field "
                                   "named 'grain_size' is present. "
                                   "In the diffusion creep regime, the viscosity depends "
                                   "on this grain size field. "
                                   "We use the grain size evolution laws described in Behn "
                                   "et al., 2009. Implications of grain size evolution on the "
                                   "seismic structure of the oceanic upper mantle, "
                                   "Earth Planet. Sci. Letters, 282, 178–189. "
                                   "Other material parameters are either prescribed similar "
                                   "to the 'simple' material model, or read from data files "
                                   "that were generated by the Perplex or Hefesto software. "
                                   "This material model "
                                   "is described in more detail in Dannberg, J., Z. Eilon, "
                                   "U. Faul, R. Gassmoeller, P. Moulik, and R. Myhill (2017), "
                                   "The importance of grain size to mantle dynamics and "
                                   "seismological observations, Geochem. Geophys. Geosyst., "
                                   "18, 3034–3061, doi:10.1002/2017GC006944.")
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {

#define INSTANTIATE(dim) \
  template class UnscaledViscosityAdditionalOutputs<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)
#undef INSTANTIATE
  }
}
