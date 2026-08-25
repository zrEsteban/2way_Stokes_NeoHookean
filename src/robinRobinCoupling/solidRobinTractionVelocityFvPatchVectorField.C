#include "solidRobinTractionVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

namespace Foam
{

makePatchTypeField
(
    fvPatchVectorField,
    solidRobinTractionVelocityFvPatchVectorField
);

solidRobinTractionVelocityFvPatchVectorField::
solidRobinTractionVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(p, iF),
    robinTraction_(p.size(), vector::zero),
    targetVelocity_(p.size(), vector::zero),
    impedance_("solidImpedance", dimDensity*dimVelocity, 0)
{}

solidRobinTractionVelocityFvPatchVectorField::
solidRobinTractionVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    solidTractionFvPatchVectorField(p, iF, dict),
    robinTraction_
    (
        dict.found("robinTraction")
      ? vectorField("robinTraction", dict, p.size())
      : vectorField(p.size(), vector::zero)
    ),
    targetVelocity_
    (
        dict.found("targetVelocity")
      ? vectorField("targetVelocity", dict, p.size())
      : vectorField(p.size(), vector::zero)
    ),
    impedance_
    (
        "solidImpedance",
        dimDensity*dimVelocity,
        dict
    )
{
    if (impedance_.dimensions() != dimDensity*dimVelocity)
    {
        FatalIOErrorInFunction(dict)
            << "solidImpedance must have dimensions "
            << dimDensity*dimVelocity << exit(FatalIOError);
    }
    if (impedance_.value() <= SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "solidImpedance must be strictly positive" << exit(FatalIOError);
    }
}

solidRobinTractionVelocityFvPatchVectorField::
solidRobinTractionVelocityFvPatchVectorField
(
    const solidRobinTractionVelocityFvPatchVectorField& rhs,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidTractionFvPatchVectorField(rhs, p, iF, mapper),
    robinTraction_(rhs.robinTraction_, mapper),
    targetVelocity_(rhs.targetVelocity_, mapper),
    impedance_(rhs.impedance_)
{}

#ifndef OPENFOAM_ORG
solidRobinTractionVelocityFvPatchVectorField::
solidRobinTractionVelocityFvPatchVectorField
(
    const solidRobinTractionVelocityFvPatchVectorField& rhs
)
:
    solidTractionFvPatchVectorField(rhs),
    robinTraction_(rhs.robinTraction_),
    targetVelocity_(rhs.targetVelocity_),
    impedance_(rhs.impedance_)
{}
#endif

solidRobinTractionVelocityFvPatchVectorField::
solidRobinTractionVelocityFvPatchVectorField
(
    const solidRobinTractionVelocityFvPatchVectorField& rhs,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(rhs, iF),
    robinTraction_(rhs.robinTraction_),
    targetVelocity_(rhs.targetVelocity_),
    impedance_(rhs.impedance_)
{}

void solidRobinTractionVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& mapper
)
{
    solidTractionFvPatchVectorField::autoMap(mapper);
    robinTraction_.autoMap(mapper);
    targetVelocity_.autoMap(mapper);
}

void solidRobinTractionVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    solidTractionFvPatchVectorField::rmap(ptf, addr);
    const solidRobinTractionVelocityFvPatchVectorField& rhs =
        refCast<const solidRobinTractionVelocityFvPatchVectorField>(ptf);
    robinTraction_.rmap(rhs.robinTraction_, addr);
    targetVelocity_.rmap(rhs.targetVelocity_, addr);
}

void solidRobinTractionVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar dt = db().time().deltaTValue();
    if (dt <= SMALL)
    {
        FatalErrorInFunction << "deltaT must be strictly positive"
            << abort(FatalError);
    }

    // DD is an incremental displacement in the active updated-Lagrangian
    // solid.  For total-displacement fields, use the old-time increment.
    vectorField solidVelocity(size(), vector::zero);
    if (internalField().name() == "DD")
    {
        solidVelocity = vectorField(*this)/dt;
    }
    else
    {
        FatalErrorInFunction
            << typeName << " currently supports the incremental displacement "
            << "field DD; received " << internalField().name()
            << abort(FatalError);
    }

    traction() =
        robinTraction_
      + impedance_.value()*(targetVelocity_ - solidVelocity);

    solidTractionFvPatchVectorField::updateCoeffs();
}

void solidRobinTractionVelocityFvPatchVectorField::write(Ostream& os) const
{
    solidTractionFvPatchVectorField::write(os);
    robinTraction_.writeEntry("robinTraction", os);
    targetVelocity_.writeEntry("targetVelocity", os);
    impedance_.writeEntry("solidImpedance", os);
}

} // End namespace Foam
