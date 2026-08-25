#ifndef dealiiPdmsSolid_TimeIntegration_H
#define dealiiPdmsSolid_TimeIntegration_H

namespace dealii_pdms
{
namespace time_integration
{

inline double velocity_coefficient(const double dt, const bool bdf2)
{
  return bdf2 ? 3.0/(2.0*dt) : 1.0/dt;
}

inline double acceleration_coefficient(const double dt, const bool bdf2)
{
  return bdf2 ? 9.0/(4.0*dt*dt) : 1.0/(dt*dt);
}

template<class Value>
Value velocity(const Value &next_displacement,
               const Value &old_displacement,
               const Value &older_displacement,
               const double dt,
               const bool bdf2)
{
  return bdf2
    ? (3.0*next_displacement-4.0*old_displacement+older_displacement)/(2.0*dt)
    : (next_displacement-old_displacement)/dt;
}

template<class Value>
Value acceleration(const Value &next_velocity,
                   const Value &old_velocity,
                   const Value &older_velocity,
                   const double dt,
                   const bool bdf2)
{
  return bdf2
    ? (3.0*next_velocity-4.0*old_velocity+older_velocity)/(2.0*dt)
    : (next_velocity-old_velocity)/dt;
}

} // namespace time_integration
} // namespace dealii_pdms

#endif
