/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

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

#include <aspect/initial_composition/lithosphere_porosity.h>
#include <aspect/initial_temperature/interface.h>
#include <aspect/postprocess/interface.h>
#include <aspect/utilities.h>
#include <deal.II/base/signaling_nan.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/material_model/interface.h>
#include <aspect/melt.h>

namespace aspect
{

  namespace InitialComposition
  {
    template <int dim>
    double
    LithospherePorosity<dim>::
    initial_composition (const Point<dim> &position,
                         const unsigned int n_comp) const
    {

      AssertThrow(this->introspection().compositional_name_exists("porosity"),
                  ExcMessage("The initial composition plugin `lithosphere porosity' did not find a "
                             "compositional field called `porosity' to initialize. Please add a "
                             "compositional field with this name."));

      const MaterialModel::MeltFractionModel<dim> *material_model =
        dynamic_cast<const MaterialModel::MeltFractionModel<dim>* > (&this->get_material_model());
      AssertThrow(material_model != nullptr,
                  ExcMessage("The used material model is not derived from the 'MeltFractionModel' class, "
                             "and therefore does not support computing equilibrium melt fractions. "
                             "This is incompatible with the `porosity' "
                             "initial composition plugin, which needs to compute these melt fractions."));

      MaterialModel::MaterialModelInputs<dim> in(1, this->n_compositional_fields());

      in.position[0] = position;
      in.temperature[0] = this->get_initial_temperature_manager().initial_temperature(position);
      in.pressure[0] = this->get_adiabatic_conditions().pressure(position);
      in.pressure_gradient[0] = 0.0;
      in.velocity[0] = 0.0;

      const double depth = this->get_geometry_model().depth(position);
      const double temperature = this->get_initial_temperature_manager().initial_temperature(position);
      const unsigned int porosity_index = this->introspection().compositional_index_for_name("porosity");
      const unsigned int peridotite_index = this->introspection().compositional_index_for_name("peridotite");


    // Crustal composition
      if (depth < moho && n_comp == 0)
        return 1.;
      else if (depth >= moho && temperature < LAB_isotherm && n_comp == 1)
        return 1.;

      else if (temperature > LAB_isotherm && n_comp == porosity_index)
        {
          in.composition[0][porosity_index] = 0.0;

          in.strain_rate[0] = SymmetricTensor<2,dim>();

          std::vector<double> melt_fraction(1);
          material_model->melt_fractions(in,melt_fraction);
          return melt_fraction[0];
        }

      else if (temperature > LAB_isotherm && n_comp == peridotite_index)
        {
          in.composition[0][peridotite_index] = 1.0;

          in.strain_rate[0] = SymmetricTensor<2,dim>();

          std::vector<double> melt_fraction(1);
          material_model->melt_fractions(in,melt_fraction);
          return 1 - melt_fraction[0];
        }


      else
        return 0.;
    }


    template <int dim>
    void
    LithospherePorosity<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("Lithosphere porosity");
        {
          prm.declare_entry ("Moho", "30000.0",
                             Patterns::Double (0),
                             "Moho depth. Units: $m.");
          prm.declare_entry ("LAB isotherm", "1673",
                             Patterns::Double (0),
                             "Temperature at the base of the lithosphere. Units: $K.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    LithospherePorosity<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("Lithosphere porosity");
        {
          moho                            = prm.get_double ("Moho");
          LAB_isotherm                    = prm.get_double ("LAB isotherm");
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
    ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(LithospherePorosity,
                                              "lithosphere porosity",
                                              "Specify the composition in terms of an explicit formula. The format of these "
                                              "functions follows the syntax understood by the "
                                              "muparser library, see Section~\\ref{sec:muparser-format}.")
  }
}
