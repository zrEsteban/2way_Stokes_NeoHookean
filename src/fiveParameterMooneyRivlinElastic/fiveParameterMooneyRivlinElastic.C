#include "fiveParameterMooneyRivlinElastic.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include <cmath>

namespace Foam
{
    defineTypeNameAndDebug(fiveParameterMooneyRivlinElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw,
        fiveParameterMooneyRivlinElastic,
        nonLinGeomMechLaw
    );
}

namespace
{
Foam::dimensionedScalar scaledCoefficient
(
    const Foam::dictionary& dict,
    const Foam::word& key,
    const Foam::scalar scale
)
{
    if (!dict.found(key))
    {
        FatalIOErrorInFunction(dict)
            << "Required five-parameter Mooney-Rivlin coefficient '" << key
            << "' is missing. Supply C10, C01, C20, C11 and C02."
            << Foam::exit(Foam::FatalIOError);
    }

    Foam::dimensionedScalar value(dict.lookup(key));
    const Foam::dimensionSet pressureDims(1, -1, -2, 0, 0, 0, 0);
    if (value.dimensions() != pressureDims)
    {
        FatalIOErrorInFunction(dict)
            << "Coefficient '" << key << "' must have pressure dimensions "
            << pressureDims << ", but has " << value.dimensions()
            << Foam::exit(Foam::FatalIOError);
    }
    return Foam::dimensionedScalar
    (
        value.name(), value.dimensions(), scale*value.value()
    );
}
}

Foam::fiveParameterMooneyRivlinElastic::fiveParameterMooneyRivlinElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    C10_(IOobject("C10", mesh.time().timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE), mesh, scaledCoefficient(dict, "C10", dict.lookupOrDefault<scalar>("coefficientScale", 1.0))),
    C01_(IOobject("C01", mesh.time().timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE), mesh, scaledCoefficient(dict, "C01", dict.lookupOrDefault<scalar>("coefficientScale", 1.0))),
    C20_(IOobject("C20", mesh.time().timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE), mesh, scaledCoefficient(dict, "C20", dict.lookupOrDefault<scalar>("coefficientScale", 1.0))),
    C11_(IOobject("C11", mesh.time().timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE), mesh, scaledCoefficient(dict, "C11", dict.lookupOrDefault<scalar>("coefficientScale", 1.0))),
    C02_(IOobject("C02", mesh.time().timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE), mesh, scaledCoefficient(dict, "C02", dict.lookupOrDefault<scalar>("coefficientScale", 1.0))),
    C10f_(fvc::interpolate(C10_)),
    C01f_(fvc::interpolate(C01_)),
    C20f_(fvc::interpolate(C20_)),
    C11f_(fvc::interpolate(C11_)),
    C02f_(fvc::interpolate(C02_)),
    coefficientScale_(dict.lookupOrDefault<scalar>("coefficientScale", 1.0))
{
    if (!std::isfinite(coefficientScale_) || coefficientScale_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "coefficientScale must be finite and strictly positive; got "
            << coefficientScale_ << exit(FatalIOError);
    }

    this->mu(2.0*(C10_ + C01_));
    if (gMin(mu()()) <= 0)
    {
        WarningInFunction
            << "The initial shear modulus 2*(C10 + C01) is non-positive "
            << "somewhere. Coefficients are retained unchanged, but the "
            << "undeformed state may be unstable." << endl;
    }

    if (dict.found("bulkModulus") && dict.found("K"))
    {
        FatalIOErrorInFunction(dict)
            << "Specify only one of bulkModulus and K" << exit(FatalIOError);
    }
    if (dict.found("bulkModulus"))
    {
        this->K(dimensionedScalar(dict.lookup("bulkModulus")));
    }
    else if (dict.found("K"))
    {
        this->K(dimensionedScalar(dict.lookup("K")));
    }
    else if (dict.found("nu"))
    {
        const dimensionedScalar nu(dict.lookup("nu"));
        this->K(2.0*mu()*(1.0 + nu)/(3.0*(1.0 - 2.0*nu)));
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Specify bulkModulus (preferred), K, or nu"
            << exit(FatalIOError);
    }

    if (gMin(K()()) <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "The bulk modulus must be strictly positive"
            << exit(FatalIOError);
    }

    Info<< "Five-parameter Mooney-Rivlin properties" << nl
        << "    coefficientScale = " << coefficientScale_ << nl
        << "    min/max(mu0) = " << gMin(mu()()) << ' '
        << gMax(mu()()) << nl
        << "    min/max(K) = " << gMin(K()()) << ' '
        << gMax(K()()) << endl;

    this->muf(fvc::interpolate(mu()));
    this->Kf(fvc::interpolate(K()));
}

Foam::tmp<Foam::volScalarField>
Foam::fiveParameterMooneyRivlinElastic::impK() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject("impK", mesh().time().timeName(), mesh(), IOobject::READ_IF_PRESENT, IOobject::NO_WRITE),
            (4.0/3.0)*mu() + K()
        )
    );
}

void Foam::fiveParameterMooneyRivlinElastic::correct(volSymmTensorField& sigma)
{
    if (updateF(sigma, mu(), K())) return;

    const volScalarField J(det(F()));
    scalar minJ = gMin(J);
    forAll(J.boundaryField(), patchI)
    {
        if (J.boundaryField()[patchI].size())
        {
            minJ = min(minJ, gMin(J.boundaryField()[patchI]));
        }
    }
    if (!std::isfinite(minJ) || minJ <= SMALL)
    {
        FatalErrorInFunction
            << "Non-positive/invalid deformation Jacobian: min(J) = " << minJ
            << ". The constitutive response is undefined for J <= 0."
            << abort(FatalError);
    }

    const volSymmTensorField isoB(pow(J, -2.0/3.0)*symm(F() & F().T()));
    const volScalarField I1(tr(isoB));
    const volScalarField I2(0.5*(sqr(I1) - tr(symm(isoB & isoB))));
    const volScalarField x(I1 - 3.0);
    const volScalarField y(I2 - 3.0);
    const volScalarField W1(C10_ + 2.0*C20_*x + C11_*y);
    const volScalarField W2(C01_ + C11_*x + 2.0*C02_*y);
    const volSymmTensorField s(2.0*W1*isoB - 2.0*W2*inv(isoB));

    updateSigmaHyd(0.5*K()*(sqr(J) - 1.0), (4.0/3.0)*mu() + K());
    sigma = (dev(s) + sigmaHyd()*I + symm(F() & sigma0() & F().T()))/J;
}

void Foam::fiveParameterMooneyRivlinElastic::correct(surfaceSymmTensorField& sigma)
{
    if (updateF(sigma, muf(), Kf())) return;

    const surfaceScalarField J(det(Ff()));
    scalar minJ = gMin(J);
    forAll(J.boundaryField(), patchI)
    {
        if (J.boundaryField()[patchI].size())
        {
            minJ = min(minJ, gMin(J.boundaryField()[patchI]));
        }
    }
    if (!std::isfinite(minJ) || minJ <= SMALL)
    {
        FatalErrorInFunction
            << "Non-positive/invalid face deformation Jacobian: min(J) = "
            << minJ << abort(FatalError);
    }

    const surfaceSymmTensorField isoB(pow(J, -2.0/3.0)*symm(Ff() & Ff().T()));
    const surfaceScalarField I1(tr(isoB));
    const surfaceScalarField I2(0.5*(sqr(I1) - tr(symm(isoB & isoB))));
    const surfaceScalarField x(I1 - 3.0);
    const surfaceScalarField y(I2 - 3.0);
    const surfaceScalarField W1(C10f_ + 2.0*C20f_*x + C11f_*y);
    const surfaceScalarField W2(C01f_ + C11f_*x + 2.0*C02f_*y);
    const surfaceSymmTensorField s(2.0*W1*isoB - 2.0*W2*inv(isoB));
    const surfaceScalarField sigmaHydf(0.5*Kf()*(sqr(J) - 1.0));

    sigma =
        (dev(s) + sigmaHydf*I
       + symm(Ff() & linearInterpolate(sigma0()) & Ff().T()))/J;
}
