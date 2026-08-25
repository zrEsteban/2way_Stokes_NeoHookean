#include "TimeIntegration.H"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

namespace ti = dealii_pdms::time_integration;

bool close(const double a, const double b, const double tolerance=1e-12)
{
  return std::abs(a-b) <= tolerance*std::max({1.0,std::abs(a),std::abs(b)});
}

int main()
{
  const double dt=0.2;

  // d=t^2 at t={0,dt,2dt}; BDF2 reproduces v=2t and a=2 exactly.
  const double d_older=0.0;
  const double d_old=dt*dt;
  const double d_next=4.0*dt*dt;
  const double v_older=0.0;
  const double v_old=2.0*dt;
  const double v_next=ti::velocity(d_next,d_old,d_older,dt,true);
  const double a_next=ti::acceleration(v_next,v_old,v_older,dt,true);
  assert(close(v_next,4.0*dt));
  assert(close(a_next,2.0));

  // Backward Euler startup on linear displacement gives exact constant speed.
  const double be_velocity=ti::velocity(0.8,0.6,0.0,dt,false);
  const double be_acceleration=ti::acceleration(be_velocity,1.0,0.0,dt,false);
  assert(close(be_velocity,1.0));
  assert(close(be_acceleration,0.0));

  // Finite-difference the inertial acceleration with respect to d^{n+1}.
  const double epsilon=1e-7;
  const auto discrete_acceleration = [&](const double displacement,
                                         const bool bdf2)
  {
    const double velocity=ti::velocity(displacement,d_old,d_older,dt,bdf2);
    return ti::acceleration(velocity,v_old,v_older,dt,bdf2);
  };
  const double bdf2_jacobian=
    (discrete_acceleration(d_next+epsilon,true)
     -discrete_acceleration(d_next-epsilon,true))/(2.0*epsilon);
  const double be_jacobian=
    (discrete_acceleration(d_next+epsilon,false)
     -discrete_acceleration(d_next-epsilon,false))/(2.0*epsilon);
  assert(close(bdf2_jacobian,ti::acceleration_coefficient(dt,true),1e-9));
  assert(close(be_jacobian,ti::acceleration_coefficient(dt,false),1e-9));

  std::cout << "PASS: BE/BDF2 velocity, acceleration, startup and inertial Jacobian"
            << std::endl;
}
