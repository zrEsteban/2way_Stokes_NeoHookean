#ifndef dealii_pdms_neo_hookean_material_h
#define dealii_pdms_neo_hookean_material_h

#include <deal.II/base/exceptions.h>
#include <deal.II/base/tensor.h>

#include <algorithm>
#include <cmath>

namespace dealii_pdms
{
using namespace dealii;

template <int dim>
Tensor<2,dim> first_piola(const Tensor<2,dim> &F,
                          const double mu,
                          const double bulk)
{
  const double J = determinant(F);
  AssertThrow(J > 0.0, ExcMessage("neo-Hookean trial has det(F) <= 0"));
  const Tensor<2,dim> FinvT = transpose(invert(F));
  const double I1 = scalar_product(F,F);
  return mu*std::pow(J,-2.0/3.0)*(F-(I1/3.0)*FinvT)
         + bulk*J*(J-1.0)*FinvT;
}

template <int dim>
Tensor<4,dim> consistent_tangent(const Tensor<2,dim> &F,
                                 const double mu,
                                 const double bulk)
{
  Tensor<4,dim> tangent;
  for (unsigned int k=0; k<dim; ++k)
    for (unsigned int l=0; l<dim; ++l)
      {
        const double h = 1e-7*std::max(1.0,std::abs(F[k][l]));
        Tensor<2,dim> plus=F, minus=F;
        plus[k][l] += h;
        minus[k][l] -= h;
        const Tensor<2,dim> derivative =
          (first_piola<dim>(plus,mu,bulk)-first_piola<dim>(minus,mu,bulk))/(2*h);
        for (unsigned int i=0; i<dim; ++i)
          for (unsigned int j=0; j<dim; ++j)
            tangent[i][j][k][l] = derivative[i][j];
      }
  return tangent;
}
}

#endif

