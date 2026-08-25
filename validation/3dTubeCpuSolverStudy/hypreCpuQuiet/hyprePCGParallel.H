#ifndef Foam_hyprePCGParallel_H
#define Foam_hyprePCGParallel_H

#include "lduMatrix.H"

namespace Foam
{

class hyprePCGParallel
:
    public lduMatrix::solver
{
    hyprePCGParallel(const hyprePCGParallel&) = delete;
    void operator=(const hyprePCGParallel&) = delete;

protected:

    solverPerformance solveKrylov
    (
        scalarField& psi,
        const scalarField& source,
        const direction cmpt,
        const bool flexGmres
    ) const;

public:

    TypeName("hyprePCGParallel");

    hyprePCGParallel
    (
        const word& fieldName,
        const lduMatrix& matrix,
        const FieldField<Field, scalar>& interfaceBouCoeffs,
        const FieldField<Field, scalar>& interfaceIntCoeffs,
        const lduInterfaceFieldPtrsList& interfaces,
        const dictionary& solverControls
    );

    virtual ~hyprePCGParallel();

    virtual solverPerformance solve
    (
        scalarField& psi,
        const scalarField& source,
        const direction cmpt=0
    ) const;
};


class hypreFlexGMRESParallel
:
    public hyprePCGParallel
{
    hypreFlexGMRESParallel(const hypreFlexGMRESParallel&) = delete;
    void operator=(const hypreFlexGMRESParallel&) = delete;

public:

    TypeName("hypreFlexGMRESParallel");

    hypreFlexGMRESParallel
    (
        const word& fieldName,
        const lduMatrix& matrix,
        const FieldField<Field, scalar>& interfaceBouCoeffs,
        const FieldField<Field, scalar>& interfaceIntCoeffs,
        const lduInterfaceFieldPtrsList& interfaces,
        const dictionary& solverControls
    );

    virtual ~hypreFlexGMRESParallel() = default;

    virtual solverPerformance solve
    (
        scalarField& psi,
        const scalarField& source,
        const direction cmpt=0
    ) const;
};

} // End namespace Foam

#endif
