#ifndef solidRobinTractionVelocityFvPatchVectorField_H
#define solidRobinTractionVelocityFvPatchVectorField_H

#include "solidTractionFvPatchVectorField.H"
#include "dimensionedScalar.H"

namespace Foam
{

class solidRobinTractionVelocityFvPatchVectorField
:
    public solidTractionFvPatchVectorField
{
    vectorField robinTraction_;
    vectorField targetVelocity_;
    dimensionedScalar impedance_;

public:

    TypeName("solidRobinTractionVelocity");

    solidRobinTractionVelocityFvPatchVectorField
    (
        const fvPatch&,
        const DimensionedField<vector, volMesh>&
    );

    solidRobinTractionVelocityFvPatchVectorField
    (
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const dictionary&
    );

    solidRobinTractionVelocityFvPatchVectorField
    (
        const solidRobinTractionVelocityFvPatchVectorField&,
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const fvPatchFieldMapper&
    );

#ifndef OPENFOAM_ORG
    solidRobinTractionVelocityFvPatchVectorField
    (
        const solidRobinTractionVelocityFvPatchVectorField&
    );
#endif

    solidRobinTractionVelocityFvPatchVectorField
    (
        const solidRobinTractionVelocityFvPatchVectorField&,
        const DimensionedField<vector, volMesh>&
    );

    virtual tmp<fvPatchVectorField> clone() const
    {
        return tmp<fvPatchVectorField>
        (
            new solidRobinTractionVelocityFvPatchVectorField(*this)
        );
    }

    virtual tmp<fvPatchVectorField> clone
    (
        const DimensionedField<vector, volMesh>& iF
    ) const
    {
        return tmp<fvPatchVectorField>
        (
            new solidRobinTractionVelocityFvPatchVectorField(*this, iF)
        );
    }

    vectorField& robinTraction() { return robinTraction_; }
    const vectorField& robinTraction() const { return robinTraction_; }
    vectorField& targetVelocity() { return targetVelocity_; }
    const vectorField& targetVelocity() const { return targetVelocity_; }
    scalar impedance() const { return impedance_.value(); }

    virtual void autoMap(const fvPatchFieldMapper&);
    virtual void rmap(const fvPatchVectorField&, const labelList&);
    virtual void updateCoeffs();
    virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif
