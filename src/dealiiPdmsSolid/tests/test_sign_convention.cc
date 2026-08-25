#include <array>
#include <cassert>
#include <iostream>

int main()
{
  // Flat roof: the fluid outward normal points from fluid into solid.
  const std::array<double,3> normal{{0.0,0.0,1.0}};
  const std::array<double,3> viscous{{0.0,0.0,0.0}};
  const double pressure=100.0;
  std::array<double,3> fluid_traction;
  std::array<double,3> solid_load;
  for (unsigned int d=0; d<3; ++d)
    {
      fluid_traction[d]=viscous[d]-pressure*normal[d];
      solid_load[d]=-fluid_traction[d];
    }
  assert(fluid_traction[2] == -pressure);
  assert(solid_load[2] == pressure);

  // A positive virtual displacement into the solid has positive external work.
  const std::array<double,3> virtual_displacement{{0.0,0.0,1.0}};
  double external_work=0.0;
  for (unsigned int d=0; d<3; ++d)
    external_work += solid_load[d]*virtual_displacement[d];
  assert(external_work > 0.0);

  std::cout << "PASS: t_f=sigma_f*n_f and f_s=-t_f give compressive +normal load"
            << std::endl;
}
