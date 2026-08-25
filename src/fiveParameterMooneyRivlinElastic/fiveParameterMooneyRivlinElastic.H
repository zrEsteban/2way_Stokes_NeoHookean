/*---------------------------------------------------------------------------*\
  Five-parameter Mooney-Rivlin law for solids4foam v2.3 / OpenFOAM v2512.
\*---------------------------------------------------------------------------*/

#ifndef fiveParameterMooneyRivlinElastic_H
#define fiveParameterMooneyRivlinElastic_H

#include "mechanicalLaw.H"

namespace Foam
{

class fiveParameterMooneyRivlinElastic
:
    public mechanicalLaw
{
    volScalarField C10_;
    volScalarField C01_;
    volScalarField C20_;
    volScalarField C11_;
    volScalarField C02_;

    surfaceScalarField C10f_;
    surfaceScalarField C01f_;
    surfaceScalarField C20f_;
    surfaceScalarField C11f_;
    surfaceScalarField C02f_;

    scalar coefficientScale_;

    fiveParameterMooneyRivlinElastic
    (
        const fiveParameterMooneyRivlinElastic&
    );
    void operator=(const fiveParameterMooneyRivlinElastic&);

public:

    TypeName("fiveParameterMooneyRivlinElastic");

    fiveParameterMooneyRivlinElastic
    (
        const word& name,
        const fvMesh& mesh,
        const dictionary& dict,
        const nonLinearGeometry::nonLinearType& nonLinGeom
    );

    virtual ~fiveParameterMooneyRivlinElastic() {}

    virtual tmp<volScalarField> impK() const;
    virtual void correct(volSymmTensorField& sigma);
    virtual void correct(surfaceSymmTensorField& sigma);
};

} // End namespace Foam

#endif

