#ifndef PROCESS_LIB_FREEZINGMATERIALMODEL_H_
#define PROCESS_LIB_FREEZINGMATERIALMODEL_H_

namespace ProcessLib
{

namespace Freezing {


static double CalcIceVolFrac(double T_in_dC, double freezing_sigmoid_coeff, double porosity)
{
   double phi_i = 0.0;

   phi_i = porosity* (1.0 - 1.0 / (1.0 + std::exp(-1.0 * freezing_sigmoid_coeff * T_in_dC)));

   return phi_i;
}

static double Calcsigmoidderive(double phi_i, double freezing_sigmoid_coeff, double porosity)
{
  double sigmoid_derive = 0.0;

  sigmoid_derive = -porosity* freezing_sigmoid_coeff * (1 - phi_i) * phi_i;

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

static double KozenyKarman(double hydraulic_conductiviy, double porosity, double phi_i)
{
    double real_hydraulic_conductivity;
    real_hydraulic_conductivity = hydraulic_conductiviy* pow((porosity-phi_i)/porosity, 3)* pow((1-porosity)/(1+phi_i-porosity), 3);
    return real_hydraulic_conductivity;
}


}  // FREEZING

}  // ProcessLib


#endif  // PROCESS_LIB_FREEZINGMATERIALMODEL_H_
