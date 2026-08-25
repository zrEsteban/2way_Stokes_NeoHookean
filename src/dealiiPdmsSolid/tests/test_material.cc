#include "NeoHookeanMaterial.H"

#include <deal.II/base/symmetric_tensor.h>

#include <cmath>
#include <iostream>

using namespace dealii;

int main()
{
  constexpr unsigned int dim=3;
  constexpr double mu=778200.0;
  constexpr double bulk=3631600.0;
  Tensor<2,dim> identity=unit_symmetric_tensor<dim>();
  const auto zero=dealii_pdms::first_piola<dim>(identity,mu,bulk);
  AssertThrow(zero.norm()<1e-8, ExcMessage("P(I) is not zero"));

  Tensor<2,dim> F=identity;
  F[0][0]=1.08;
  F[1][1]=0.97;
  F[2][2]=1.02;
  F[0][2]=0.04;
  const auto tangent=dealii_pdms::consistent_tangent<dim>(F,mu,bulk);
  double major_defect=0;
  for (unsigned int i=0;i<dim;++i)
    for (unsigned int j=0;j<dim;++j)
      for (unsigned int k=0;k<dim;++k)
        for (unsigned int l=0;l<dim;++l)
          major_defect=std::max(major_defect,
            std::abs(tangent[i][j][k][l]-tangent[k][l][i][j]));
  AssertThrow(major_defect<1e-2, ExcMessage("material tangent lost major symmetry"));

  bool rejected=false;
  F=identity;
  F[2][2]=-0.1;
  try { (void)dealii_pdms::first_piola<dim>(F,mu,bulk); }
  catch (const std::exception &) { rejected=true; }
  AssertThrow(rejected, ExcMessage("det(F)<=0 was not rejected"));
  std::cout << "PASS: neo-Hookean identity, tangent symmetry, and J gate; defect="
            << major_defect << std::endl;
}

