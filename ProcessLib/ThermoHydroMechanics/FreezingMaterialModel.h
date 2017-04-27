#ifndef PROCESS_LIB_FREEZINGMATERIALMODEL_H_
#define PROCESS_LIB_FREEZINGMATERIALMODEL_H_

#include <math.h>

namespace ProcessLib
{

namespace ThermoHydroMechanics {


static double CalcIceVolFrac(double T_in_dC, double freezing_sigmoid_coeff, double porosity)
{
   double phi_i = 0.0;

   phi_i = porosity* (1.0 - 1.0 / (1.0 + exp(-1.0 * freezing_sigmoid_coeff * T_in_dC)));

   return phi_i;
}

// dPhi_i/dT
static double Calcsigmoidderive(double freezing_sigmoid_coeff, double porosity, double T_in_dC)
{
  double sigmoid_derive = 0.0;
  double logistic = 1.0 / (1.0 + exp(-1.0 * freezing_sigmoid_coeff * T_in_dC));

  sigmoid_derive = -porosity* freezing_sigmoid_coeff * (1 - logistic) * logistic;

  return sigmoid_derive;
}

// second derivative the negative symbol is cancelled out
static double Calcsigmoidsecondderive(double freezing_sigmoid_coeff, double porosity, double T_in_dC)
{
  double sigmoid_derive = 0.0;
  double logistic = 1.0 / (1.0 + exp(-1.0 * freezing_sigmoid_coeff * T_in_dC));

  sigmoid_derive = porosity* freezing_sigmoid_coeff * freezing_sigmoid_coeff*
          (1 - logistic) * logistic *(1 - 2*logistic);

  return sigmoid_derive;
}

static double EquaHeatCapacity(double phi_i, double density_w, double density_s, double density_i,
double specific_heat_capacity_soil, double specific_heat_capacity_ice, double specific_heat_capacity_water,
double porosity, double sigmoid_derive, double latent_heat)
{
    double heat_capacity;

    heat_capacity = (porosity - phi_i) *specific_heat_capacity_water*density_w + (1.0 - porosity) *specific_heat_capacity_soil* density_s
    + phi_i * specific_heat_capacity_ice * density_i - density_i * sigmoid_derive * latent_heat ;

    return heat_capacity;
}

static double TotalThermalConductivity(double porosity, double phi_i, double thermal_conductivity_ice,
double thermal_conductivity_soil, double thermal_conductivity_water)
{
   double thermal_conductivity;

   thermal_conductivity = (porosity - phi_i)*thermal_conductivity_water + (1 - porosity)*thermal_conductivity_soil
           +  phi_i*thermal_conductivity_ice ;

   return thermal_conductivity;
}

/*static double KozenyKarman(double hydraulic_conductiviy, double porosity, double phi_i)
{
    double real_hydraulic_conductivity;
    real_hydraulic_conductivity = hydraulic_conductiviy* pow((porosity-phi_i)/porosity, 3)* pow((1-porosity)/(1+phi_i-porosity), 3);
    return real_hydraulic_conductivity;
}*/


static double Relative_permeability(double phi_i, double porosity)
{
    double Sw = (porosity - phi_i)/porosity;
    double Omega = 50 ;
    double Kr;
    Kr = pow(10, -Omega*porosity*(1-Sw));
    if (Kr < pow(10, -4))
            Kr = pow(10, -4);
    return Kr;
}

static double Hydraulic_conductivity(double Kint, double Kr, double mu)
{
    Kr = 1.0 ;
    return Kint*Kr/mu;
}

static double DensityWater_T(double density0, double temperature, double temperature0, double beta)
{
   return density0*(1 - beta*(temperature - temperature0));
}

static double Viscosity(double viscosity0, double temperature, double temperature_con, double temperature0)
{
    return viscosity0*exp(-(temperature - temperature0)/temperature_con); 
}

}  // FREEZING

}  // ProcessLib


#endif  // PROCESS_LIB_FREEZINGMATERIALMODEL_H_
