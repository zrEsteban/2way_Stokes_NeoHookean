#include "robinRobinCouplingInterface.H"
#include "RectangularMatrix.H"
#include "solidRobinTractionVelocityFvPatchVectorField.H"
#include "pdmsElasticWallPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"

#include <cstdlib>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

namespace Foam
{
namespace fluidSolidInterfaces
{

defineTypeNameAndDebug(robinRobinCouplingInterface, 0);
addToRunTimeSelectionTable
(
    fluidSolidInterface,
    robinRobinCouplingInterface,
    dictionary
);

robinRobinCouplingInterface::robinRobinCouplingInterface
(
    Time& runTime,
    const word& region
)
:
    fluidSolidInterface(typeName, runTime, region),
    relaxationFactor_
    (
        fsiProperties().lookupOrAddDefault<scalar>("relaxationFactor", 1.0)
    ),
    couplingScheme_
    (
        fsiProperties().lookupOrAddDefault<word>
        (
            "couplingScheme", "robinRobin"
        )
    ),
    semiImplicitBeta_
    (
        fsiProperties().lookupOrAddDefault<scalar>("semiImplicitBeta", 1.0)
    ),
    accelerationMethod_
    (
        fsiProperties().lookupOrAddDefault<word>
        (
            "accelerationMethod", "fixed"
        )
    ),
    relaxationFactorMin_
    (
        fsiProperties().lookupOrAddDefault<scalar>
        (
            "relaxationFactorMin", 0.01
        )
    ),
    relaxationFactorMax_
    (
        fsiProperties().lookupOrAddDefault<scalar>
        (
            "relaxationFactorMax", 1.0
        )
    ),
    aitkenRelaxationFactors_(nGlobalPatches(), relaxationFactor_),
    iqnMaxModes_
    (
        fsiProperties().lookupOrAddDefault<label>("iqnMaxModes", 20)
    ),
    iqnMaxCorrectionFactor_
    (
        fsiProperties().lookupOrAddDefault<scalar>
        (
            "iqnMaxCorrectionFactor", 2.0
        )
    ),
    iqnReuseAcrossTimeSteps_
    (
        fsiProperties().lookupOrAddDefault<bool>
        (
            "iqnReuseAcrossTimeSteps", false
        )
    ),
    iqnModeFilterTolerance_
    (
        fsiProperties().lookupOrAddDefault<scalar>
        (
            "iqnModeFilterTolerance", 1e-10
        )
    ),
    iqnV_(nGlobalPatches()),
    iqnW_(nGlobalPatches()),
    iqnCurrentStepModes_(nGlobalPatches(), 0),
    iqnDataV_(nGlobalPatches()),
    iqnDataW_(nGlobalPatches()),
    iqnDataResidualPrev_(nGlobalPatches()),
    iqnDataOutputPrev_(nGlobalPatches()),
    robinDataRelaxationFactor_
    (
        fsiProperties().lookupOrAddDefault<scalar>
        (
            "robinDataRelaxationFactor", 1.0
        )
    ),
    predictSolid_(fsiProperties().lookupOrAddDefault<bool>("predictSolid", false)),
    velocityTolerance_
    (
        fsiProperties().lookupOrAddDefault<scalar>
        (
            "robinVelocityTolerance", 1e-4
        )
    ),
    forceTolerance_
    (
        fsiProperties().lookupOrAddDefault<scalar>("robinForceTolerance", 1e-4)
    ),
    powerTolerance_
    (
        fsiProperties().lookupOrAddDefault<scalar>("robinPowerTolerance", 1e-3)
    ),
    velocityScaleFloor_
    (
        fsiProperties().lookupOrAddDefault<scalar>("robinVelocityScale", 0)
    ),
    tractionScaleFloor_
    (
        fsiProperties().lookupOrAddDefault<scalar>("robinTractionScale", 0)
    ),
    transferredTraction_(),
    transferredVelocity_(),
    useDealII_(fsiProperties().lookupOrAddDefault<bool>("useDealII", false)),
    dealIICommand_
    (
        fsiProperties().lookupOrAddDefault<string>("dealIICommand", "")
    ),
    dealIIAcceptCommand_
    (
        fsiProperties().lookupOrAddDefault<string>("dealIIAcceptCommand", "")
    ),
    dealIIInput_
    (
        fsiProperties().lookupOrAddDefault<fileName>("dealIIInput", "robin-in.csv")
    ),
    dealIIQueryInput_
    (
        fsiProperties().lookupOrAddDefault<fileName>
        (
            "dealIIQueryInput", "robin-query.csv"
        )
    ),
    dealIIOutput_
    (
        fsiProperties().lookupOrAddDefault<fileName>("dealIIOutput", "robin-out.csv")
    ),
    dealIIImpedance_
    (
        fsiProperties().lookupOrAddDefault<scalar>("dealIIImpedance", 0)
    ),
    stateAudit_
    (
        fsiProperties().lookupOrAddDefault<bool>("stateAudit", false)
    ),
    acceptedStepCount_(0),
    dealIIPointIncrement_(),
    dealIIFaceVelocity_(),
    dealIIFaceTraction_(),
    dealIIFaceAcceleration_()
{
    if
    (
        couplingScheme_ != "robinRobin"
     && couplingScheme_ != "dirichletNeumann"
     && couplingScheme_ != "semiImplicitBeta"
    )
    {
        FatalIOErrorInFunction(fsiProperties())
            << "couplingScheme must be robinRobin, dirichletNeumann or "
            << "semiImplicitBeta"
            << exit(FatalIOError);
    }
    if (relaxationFactor_ <= 0 || relaxationFactor_ > 1)
    {
        FatalIOErrorInFunction(fsiProperties())
            << "relaxationFactor must be in (0,1]" << exit(FatalIOError);
    }
    if
    (
        couplingScheme_ == "semiImplicitBeta"
     && mag(semiImplicitBeta_ - 1.0) > SMALL
    )
    {
        FatalIOErrorInFunction(fsiProperties())
            << "The progressive semiImplicitBeta gate currently supports "
            << "only beta=1" << exit(FatalIOError);
    }
    if
    (
        accelerationMethod_ != "fixed"
     && accelerationMethod_ != "aitken"
     && accelerationMethod_ != "iqnils"
     && accelerationMethod_ != "iqnilsCombined"
    )
    {
        FatalIOErrorInFunction(fsiProperties())
            << "accelerationMethod must be fixed, aitken, iqnils or "
            << "iqnilsCombined"
            << exit(FatalIOError);
    }
    if
    (
        iqnMaxModes_ < 1
     || iqnMaxCorrectionFactor_ <= 0
     || iqnModeFilterTolerance_ <= 0
    )
    {
        FatalIOErrorInFunction(fsiProperties())
            << "IQN-ILS mode count, correction factor and filter tolerance "
            << "must be positive"
            << exit(FatalIOError);
    }
    if
    (
        relaxationFactorMin_ <= 0
     || relaxationFactorMax_ < relaxationFactorMin_
    )
    {
        FatalIOErrorInFunction(fsiProperties())
            << "Require 0 < relaxationFactorMin <= relaxationFactorMax"
            << exit(FatalIOError);
    }

    transferredTraction_.setSize(nGlobalPatches());
    transferredVelocity_.setSize(nGlobalPatches());
    dealIIPointIncrement_.setSize(nGlobalPatches());
    dealIIFaceVelocity_.setSize(nGlobalPatches());
    dealIIFaceTraction_.setSize(nGlobalPatches());
    dealIIFaceAcceleration_.setSize(nGlobalPatches());
    if
    (
        useDealII_
     &&
        (
            dealIICommand_.empty()
         || dealIIImpedance_ < 0
         || (couplingScheme_ == "robinRobin" && dealIIImpedance_ <= SMALL)
         ||
            (
                couplingScheme_ != "robinRobin"
             && dealIIImpedance_ > SMALL
            )
        )
    )
    {
        FatalIOErrorInFunction(fsiProperties())
            << "useDealII requires a command, positive impedance for "
            << "robinRobin, and zero impedance for the Neumann structural "
            << "participant in dirichletNeumann/semiImplicitBeta"
            << exit(FatalIOError);
    }
    if
    (
        robinDataRelaxationFactor_ <= 0
     || robinDataRelaxationFactor_ > 1
    )
    {
        FatalIOErrorInFunction(fsiProperties())
            << "robinDataRelaxationFactor must be in (0,1]"
            << exit(FatalIOError);
    }
    if
    (
        velocityTolerance_ <= 0
     || forceTolerance_ <= 0
     || powerTolerance_ <= 0
    )
    {
        FatalIOErrorInFunction(fsiProperties())
            << "Robin physical tolerances must be strictly positive"
            << exit(FatalIOError);
    }
}

void robinRobinCouplingInterface::solveDealIISolid
(
    const label interfaceI,
    const standAlonePatch& fluidZone,
    const vectorField& fluidTraction,
    const vectorField& fluidVelocity
)
{
    if (Pstream::master())
    {
        std::ofstream output(dealIIInput_.c_str());
        if (!output.good())
        {
            FatalErrorInFunction << "Cannot write " << dealIIInput_
                << abort(FatalError);
        }
        output << "x,y,z,tx,ty,tz,vx,vy,vz\n" << std::setprecision(17);
        const vectorField centres(fluidZone.faceCentres());
        forAll(centres, faceI)
        {
            const vector solidT = -fluidTraction[faceI];
            output << centres[faceI].x() << ',' << centres[faceI].y() << ','
                << centres[faceI].z() << ',' << solidT.x() << ','
                << solidT.y() << ',' << solidT.z() << ','
                << fluidVelocity[faceI].x() << ',' << fluidVelocity[faceI].y()
                << ',' << fluidVelocity[faceI].z() << '\n';
        }
        output.close();
        std::ofstream query(dealIIQueryInput_.c_str());
        if (!query.good())
        {
            FatalErrorInFunction << "Cannot write " << dealIIQueryInput_
                << abort(FatalError);
        }
        query << "x,y,z\n" << std::setprecision(17);
        const pointField& points=fluidZone.localPoints();
        forAll(points,pointI)
        {
            query << points[pointI].x() << ',' << points[pointI].y() << ','
                << points[pointI].z() << '\n';
        }
        query.close();
        const int status = std::system(dealIICommand_.c_str());
        if (status != 0)
        {
            FatalErrorInFunction << "deal.II command failed with status "
                << status << ": " << dealIICommand_ << abort(FatalError);
        }
    }
    UPstream::barrier(UPstream::worldComm);

    std::ifstream input(dealIIOutput_.c_str());
    if (!input.good())
    {
        FatalErrorInFunction << "Cannot read " << dealIIOutput_
            << abort(FatalError);
    }
    std::string line;
    std::getline(input,line);
    std::vector<point> locations;
    std::vector<vector> velocities;
    std::vector<vector> tractions;
    std::vector<vector> accelerations;
    while (std::getline(input,line))
    {
        std::replace(line.begin(),line.end(),',',' ');
        std::istringstream row(line);
        scalar x,y,z,ux,uy,uz,vx,vy,vz,tx,ty,tz,ax,ay,az;
        if (row >> x >> y >> z >> ux >> uy >> uz >> vx >> vy >> vz
                >> tx >> ty >> tz >> ax >> ay >> az)
        {
            locations.push_back(point(x,y,z));
            velocities.push_back(vector(vx,vy,vz));
            tractions.push_back(vector(tx,ty,tz));
            accelerations.push_back(vector(ax,ay,az));
        }
    }
    if (locations.empty())
    {
        FatalErrorInFunction << dealIIOutput_ << " has no interface samples"
            << abort(FatalError);
    }

    const pointField& fluidPoints = fluidZone.localPoints();
    dealIIPointIncrement_[interfaceI].setSize(fluidPoints.size());
    forAll(fluidPoints, pointI)
    {
        label nearestI=0;
        scalar nearestDistance=GREAT;
        for (label sampleI=0; sampleI<label(locations.size()); ++sampleI)
        {
            const scalar distance=magSqr(fluidPoints[pointI]-locations[sampleI]);
            if (distance < nearestDistance)
            {
                nearestDistance=distance;
                nearestI=sampleI;
            }
        }
        dealIIPointIncrement_[interfaceI][pointI] =
            runTime().deltaTValue()*velocities[nearestI];
    }

    const vectorField fluidCentres(fluidZone.faceCentres());
    dealIIFaceVelocity_[interfaceI].setSize(fluidCentres.size());
    dealIIFaceTraction_[interfaceI].setSize(fluidCentres.size());
    dealIIFaceAcceleration_[interfaceI].setSize(fluidCentres.size());
    forAll(fluidCentres, faceI)
    {
        const face& interfaceFace=fluidZone[faceI];
        vector solidVelocity=vector::zero;
        vector solidAcceleration=vector::zero;
        forAll(interfaceFace,facePointI)
        {
            const label pointI=interfaceFace[facePointI];
            solidVelocity += dealIIPointIncrement_[interfaceI][pointI]
                /runTime().deltaTValue();
            label nearestI=0;
            scalar nearestDistance=GREAT;
            for (label sampleI=0; sampleI<label(locations.size()); ++sampleI)
            {
                const scalar distance=
                    magSqr(fluidPoints[pointI]-locations[sampleI]);
                if (distance < nearestDistance)
                {
                    nearestDistance=distance;
                    nearestI=sampleI;
                }
            }
            solidAcceleration += accelerations[nearestI];
        }
        solidVelocity /= interfaceFace.size();
        solidAcceleration /= interfaceFace.size();
        dealIIFaceVelocity_[interfaceI][faceI]=solidVelocity;
        dealIIFaceAcceleration_[interfaceI][faceI]=solidAcceleration;
        dealIIFaceTraction_[interfaceI][faceI]=
            -fluidTraction[faceI]
          + dealIIImpedance_*(fluidVelocity[faceI]-solidVelocity);
    }
}

void robinRobinCouplingInterface::updateDisplacement()
{
    Info<< nl << "Time = " << fluid().runTime().timeName()
        << ", " << couplingScheme_ << " iteration: " << outerCorr() << endl;

    forAll(fluid().globalPatches(), interfaceI)
    {
        if (accelerationMethod_ == "iqnils")
        {
            if (outerCorr() == 1)
            {
                iqnCurrentStepModes_[interfaceI]=0;
                if (!iqnReuseAcrossTimeSteps_)
                {
                    iqnV_[interfaceI].clear();
                    iqnW_[interfaceI].clear();
                }
            }
            else if (outerCorr() == 2)
            {
                solidZonesPointsDisplsRef()[interfaceI]=
                    solidZonesPointsDispls()[interfaceI];
                fluidZonesPointsDisplsRef()[interfaceI]=
                    fluidZonesPointsDispls()[interfaceI];
            }
            else
            {
                const vectorField newV
                (
                    residuals()[interfaceI]
                  - (
                        solidZonesPointsDisplsRef()[interfaceI]
                      - fluidZonesPointsDisplsRef()[interfaceI]
                    )
                );
                const vectorField newW
                (
                    solidZonesPointsDispls()[interfaceI]
                  - solidZonesPointsDisplsRef()[interfaceI]
                );
                vectorField orthogonalPart(newV);
                const scalar originalNorm=Foam::sqrt(gSum(newV & newV));
                forAll(iqnV_[interfaceI], modeI)
                {
                    const scalar modeNormSqr=
                        gSum(iqnV_[interfaceI][modeI]
                           & iqnV_[interfaceI][modeI]);
                    if (modeNormSqr > VSMALL)
                    {
                        orthogonalPart -= iqnV_[interfaceI][modeI]
                           *(gSum(iqnV_[interfaceI][modeI]
                               & orthogonalPart)/modeNormSqr);
                    }
                }
                const scalar independentNorm=
                    Foam::sqrt(gSum(orthogonalPart & orthogonalPart));
                if
                (
                    originalNorm > VSMALL
                 && independentNorm
                    > iqnModeFilterTolerance_*originalNorm
                )
                {
                    iqnV_[interfaceI].append(newV);
                    iqnW_[interfaceI].append(newW);
                    ++iqnCurrentStepModes_[interfaceI];
                }
                else
                {
                    Info<< couplingScheme_ << " interface " << interfaceI
                        << " rejected nearly dependent IQN mode" << endl;
                }
                if (iqnV_[interfaceI].size() > iqnMaxModes_)
                {
                    for (label i=0; i<iqnV_[interfaceI].size()-1; ++i)
                    {
                        iqnV_[interfaceI][i]=iqnV_[interfaceI][i+1];
                        iqnW_[interfaceI][i]=iqnW_[interfaceI][i+1];
                    }
                    iqnV_[interfaceI].remove();
                    iqnW_[interfaceI].remove();
                }
            }

            fluidZonesPointsDisplsPrev()[interfaceI]=
                fluidZonesPointsDispls()[interfaceI];
            if
            (
                iqnV_[interfaceI].size() > 1
             && iqnCurrentStepModes_[interfaceI] > 1
            )
            {
                const label columns=iqnV_[interfaceI].size();
                RectangularMatrix<scalar> R(columns,columns,0);
                RectangularMatrix<scalar> coefficients(columns,1,0);
                DynamicList<vectorField> Q;
                for (label i=0; i<columns; ++i)
                {
                    Q.append(iqnV_[interfaceI][columns-1-i]);
                }
                scalar maximumDiagonal=0;
                for (label i=0; i<columns; ++i)
                {
                    R[i][i]=Foam::sqrt(gSum(Q[i] & Q[i]));
                    maximumDiagonal=max(maximumDiagonal,mag(R[i][i]));
                    if (R[i][i] > VSMALL) Q[i] /= R[i][i];
                    for (label j=i+1; j<columns; ++j)
                    {
                        R[i][j]=gSum(Q[i] & Q[j]);
                        Q[j] -= R[i][j]*Q[i];
                    }
                    coefficients[i][0]=gSum
                    (
                        Q[i]
                      & (
                            fluidZonesPointsDispls()[interfaceI]
                          - solidZonesPointsDispls()[interfaceI]
                        )
                    );
                }
                const scalar cutoff=max(1e-12*maximumDiagonal,VSMALL);
                for (label i=0; i<columns; ++i)
                {
                    if (mag(R[i][i]) > cutoff)
                    {
                        for (label j=i+1; j<columns; ++j)
                            R[i][j] /= R[i][i];
                        coefficients[i][0] /= R[i][i];
                        R[i][i]=1;
                    }
                }
                for (label j=columns-1; j>=0; --j)
                {
                    if (mag(R[j][j]) > cutoff)
                    {
                        for (label i=0; i<j; ++i)
                            coefficients[i][0] -=
                                coefficients[j][0]*R[i][j];
                    }
                    else coefficients[j][0]=0;
                }
                vectorField candidate
                (
                    solidZonesPointsDispls()[interfaceI]
                );
                for (label i=0; i<columns; ++i)
                {
                    candidate +=
                        iqnW_[interfaceI][i]
                       *coefficients[columns-1-i][0];
                }
                vectorField correction
                (
                    candidate-fluidZonesPointsDispls()[interfaceI]
                );
                const scalar correctionNorm=gMax(mag(correction));
                const scalar correctionLimit=iqnMaxCorrectionFactor_
                   *gMax(mag(residuals()[interfaceI]));
                if
                (
                    correctionNorm > correctionLimit
                 && correctionNorm > SMALL
                )
                {
                    correction *= correctionLimit/correctionNorm;
                }
                fluidZonesPointsDispls()[interfaceI] += correction;
                Info<< couplingScheme_ << " interface " << interfaceI
                    << " IQN-ILS modes = " << columns << endl;
            }
            else
            {
                fluidZonesPointsDispls()[interfaceI] +=
                    relaxationFactor_*residuals()[interfaceI];
                Info<< couplingScheme_ << " interface " << interfaceI
                    << " IQN-ILS startup factor = " << relaxationFactor_
                    << endl;
            }
            continue;
        }
        scalar factor=relaxationFactor_;
        if (accelerationMethod_ == "aitken")
        {
            if (outerCorr() <= 2)
            {
                aitkenRelaxationFactors_[interfaceI]=relaxationFactor_;
            }
            else
            {
                const vectorField deltaResidual
                (
                    residuals()[interfaceI]-residualsPrev()[interfaceI]
                );
                const scalar denominator=gSum
                (
                    deltaResidual & deltaResidual
                );
                if (denominator > VSMALL)
                {
                    aitkenRelaxationFactors_[interfaceI] *= -gSum
                    (
                        residualsPrev()[interfaceI] & deltaResidual
                    )/denominator;
                }
                aitkenRelaxationFactors_[interfaceI]=min
                (
                    relaxationFactorMax_,
                    max
                    (
                        relaxationFactorMin_,
                        mag(aitkenRelaxationFactors_[interfaceI])
                    )
                );
            }
            factor=aitkenRelaxationFactors_[interfaceI];
        }
        Info<< couplingScheme_ << " interface " << interfaceI
            << " relaxation factor (" << accelerationMethod_ << ") = "
            << factor << endl;
        fluidZonesPointsDisplsPrev()[interfaceI] =
            fluidZonesPointsDispls()[interfaceI];
        fluidZonesPointsDispls()[interfaceI] +=
            factor*residuals()[interfaceI];
    }

    updateMovingWallPressureAcceleration();
    updateElasticWallPressureAcceleration();
    updatePdmsElasticWallPressure();
    syncFluidZonePointsDispl(fluidZonesPointsDispls());
}

void robinRobinCouplingInterface::updatePdmsElasticWallPressure()
{
    forAll(fluid().globalPatches(), interfaceI)
    {
        const label fluidPatchID = fluidPatchIndices()[interfaceI];
        fvPatchScalarField& patchP =
            fluid().solutionP().boundaryFieldRef()[fluidPatchID];

        if
        (
            patchP.type()
         != pdmsElasticWallPressureFvPatchScalarField::typeName
        )
        {
            continue;
        }

        const standAlonePatch& fluidZone =
            fluid().globalPatches()[interfaceI].globalPatch();
        vectorField fluidAcceleration(fluidZone.size(), vector::zero);
        if (useDealII_ && dealIIFaceAcceleration_[interfaceI].size())
        {
            fluidAcceleration=dealIIFaceAcceleration_[interfaceI];
        }
        else if (!useDealII_)
        {
            const standAlonePatch& solidZone =
                solid().globalPatches()[interfaceI].globalPatch();
            const vectorField solidAcceleration
            (
                solid().faceZoneAcceleration(interfaceI)
            );
            interfaceToInterfaceList()[interfaceI].transferFacesZoneToZone
            (
                solidZone, fluidZone, solidAcceleration, fluidAcceleration
            );
        }

        pdmsElasticWallPressureFvPatchScalarField& robinP =
            refCast<pdmsElasticWallPressureFvPatchScalarField>(patchP);
        robinP.prevAcceleration() =
            fluid().globalPatches()[interfaceI].globalFaceToPatch
            (
                fluidAcceleration
            );
        robinP.prevPressure() =
            fluid().patchSolutionPressureForce(fluidPatchID);
    }
}

void robinRobinCouplingInterface::updateForce()
{
    Info<< "Setting Robin traction and target velocity on solid interfaces"
        << endl;

    for (label interfaceI = 0; interfaceI < nGlobalPatches(); ++interfaceI)
    {
        const standAlonePatch& fluidZone =
            fluid().globalPatches()[interfaceI].globalPatch();
        const standAlonePatch& solidZone =
            solid().globalPatches()[interfaceI].globalPatch();

        const vectorField fluidTraction
        (
            fluid().faceZoneViscousForce(interfaceI)
          - fluid().faceZonePressureForce(interfaceI)*fluidZone.faceNormals()
        );
        vectorField solidTraction(solidZone.size(), vector::zero);
        interfaceToInterfaceList()[interfaceI].transferFacesZoneToZone
        (
            fluidZone, solidZone, fluidTraction, solidTraction
        );
        solidTraction = -solidTraction;

        const label fluidPatchID = fluidPatchIndices()[interfaceI];
        const vectorField fluidPatchVelocity
        (
            fluid().U().boundaryField()[fluidPatchID]
        );
        const vectorField fluidZoneVelocity
        (
            fluid().globalPatches()[interfaceI].patchFaceToGlobal
            (
                fluidPatchVelocity
            )
        );
        vectorField solidZoneVelocity(solidZone.size(), vector::zero);
        interfaceToInterfaceList()[interfaceI].transferFacesZoneToZone
        (
            fluidZone, solidZone, fluidZoneVelocity, solidZoneVelocity
        );

        if (useDealII_)
        {
            if (accelerationMethod_ == "iqnilsCombined")
            {
                const label nFaces=fluidTraction.size();
                vectorField output(2*nFaces, vector::zero);
                for (label faceI=0; faceI<nFaces; ++faceI)
                {
                    output[faceI]=fluidZoneVelocity[faceI]
                       /(velocityScaleFloor_+SMALL);
                    output[nFaces+faceI]=fluidTraction[faceI]
                       /(tractionScaleFloor_+SMALL);
                }

                if
                (
                    outerCorr() == 1
                 || transferredTraction_[interfaceI].size() != nFaces
                )
                {
                    iqnDataV_[interfaceI].clear();
                    iqnDataW_[interfaceI].clear();
                    transferredTraction_[interfaceI]=fluidTraction;
                    transferredVelocity_[interfaceI]=fluidZoneVelocity;
                    iqnDataResidualPrev_[interfaceI].clear();
                    iqnDataOutputPrev_[interfaceI].clear();
                }
                else
                {
                    vectorField input(2*nFaces, vector::zero);
                    for (label faceI=0; faceI<nFaces; ++faceI)
                    {
                        input[faceI]=transferredVelocity_[interfaceI][faceI]
                           /(velocityScaleFloor_+SMALL);
                        input[nFaces+faceI]=
                            transferredTraction_[interfaceI][faceI]
                           /(tractionScaleFloor_+SMALL);
                    }
                    const vectorField dataResidual(output-input);

                    if (iqnDataResidualPrev_[interfaceI].size())
                    {
                        iqnDataV_[interfaceI].append
                        (
                            vectorField
                            (
                                dataResidual
                              - iqnDataResidualPrev_[interfaceI]
                            )
                        );
                        iqnDataW_[interfaceI].append
                        (
                            vectorField
                            (
                                output-iqnDataOutputPrev_[interfaceI]
                            )
                        );
                        if (iqnDataV_[interfaceI].size() > iqnMaxModes_)
                        {
                            for
                            (
                                label i=0;
                                i<iqnDataV_[interfaceI].size()-1;
                                ++i
                            )
                            {
                                iqnDataV_[interfaceI][i]=
                                    iqnDataV_[interfaceI][i+1];
                                iqnDataW_[interfaceI][i]=
                                    iqnDataW_[interfaceI][i+1];
                            }
                            iqnDataV_[interfaceI].remove();
                            iqnDataW_[interfaceI].remove();
                        }
                    }

                    iqnDataResidualPrev_[interfaceI]=dataResidual;
                    iqnDataOutputPrev_[interfaceI]=output;
                    vectorField candidate(input);
                    const label columns=iqnDataV_[interfaceI].size();
                    if (columns > 1)
                    {
                        RectangularMatrix<scalar> R(columns,columns,0);
                        RectangularMatrix<scalar> C(columns,1,0);
                        DynamicList<vectorField> Q;
                        for (label i=0; i<columns; ++i)
                            Q.append(iqnDataV_[interfaceI][columns-1-i]);
                        scalar maximumDiagonal=0;
                        for (label i=0; i<columns; ++i)
                        {
                            R[i][i]=Foam::sqrt(gSum(Q[i] & Q[i]));
                            maximumDiagonal=max(maximumDiagonal,mag(R[i][i]));
                            if (R[i][i] > VSMALL) Q[i] /= R[i][i];
                            for (label j=i+1; j<columns; ++j)
                            {
                                R[i][j]=gSum(Q[i] & Q[j]);
                                Q[j] -= R[i][j]*Q[i];
                            }
                            C[i][0]=gSum(Q[i] & (-dataResidual));
                        }
                        const scalar cutoff=
                            max(1e-12*maximumDiagonal,VSMALL);
                        for (label i=0; i<columns; ++i)
                        {
                            if (mag(R[i][i]) > cutoff)
                            {
                                for (label j=i+1; j<columns; ++j)
                                    R[i][j] /= R[i][i];
                                C[i][0] /= R[i][i];
                                R[i][i]=1;
                            }
                        }
                        for (label j=columns-1; j>=0; --j)
                        {
                            if (mag(R[j][j]) > cutoff)
                            {
                                for (label i=0; i<j; ++i)
                                    C[i][0] -= C[j][0]*R[i][j];
                            }
                            else C[j][0]=0;
                        }
                        candidate=output;
                        for (label i=0; i<columns; ++i)
                            candidate += iqnDataW_[interfaceI][i]
                               *C[columns-1-i][0];
                    }
                    else candidate=input+relaxationFactor_*dataResidual;

                    vectorField correction(candidate-input);
                    const scalar correctionNorm=gMax(mag(correction));
                    const scalar correctionLimit=iqnMaxCorrectionFactor_
                       *gMax(mag(dataResidual));
                    if
                    (
                        correctionNorm > correctionLimit
                     && correctionNorm > SMALL
                    )
                        candidate=input
                           + correction*(correctionLimit/correctionNorm);
                    for (label faceI=0; faceI<nFaces; ++faceI)
                    {
                        transferredVelocity_[interfaceI][faceI]=
                            candidate[faceI]*velocityScaleFloor_;
                        transferredTraction_[interfaceI][faceI]=
                            candidate[nFaces+faceI]*tractionScaleFloor_;
                    }
                    Info<< "deal.II combined IQN-ILS modes = "
                        << columns << endl;
                }
                solveDealIISolid
                (
                    interfaceI, fluidZone,
                    transferredTraction_[interfaceI],
                    transferredVelocity_[interfaceI]
                );
                continue;
            }
            const scalar dataFactor=robinDataRelaxationFactor_;
            if (transferredTraction_[interfaceI].size() != fluidTraction.size())
            {
                transferredTraction_[interfaceI]=fluidTraction;
                transferredVelocity_[interfaceI]=fluidZoneVelocity;
            }
            else
            {
                transferredTraction_[interfaceI]=
                    (1-dataFactor)
                   *transferredTraction_[interfaceI]
                  + dataFactor*fluidTraction;
                transferredVelocity_[interfaceI]=
                    (1-dataFactor)
                   *transferredVelocity_[interfaceI]
                  + dataFactor*fluidZoneVelocity;
            }
            solveDealIISolid
            (
                interfaceI, fluidZone, transferredTraction_[interfaceI],
                transferredVelocity_[interfaceI]
            );
            Info<< "deal.II Robin data relaxation factor = "
                << dataFactor << endl;
            continue;
        }

        const label solidPatchID = solidPatchIndices()[interfaceI];
        fvPatchVectorField& patchD =
            solid().solutionD().boundaryFieldRef()[solidPatchID];

        if
        (
            patchD.type()
         != solidRobinTractionVelocityFvPatchVectorField::typeName
        )
        {
            FatalErrorInFunction
                << "Solid interface patch " << patchD.patch().name()
                << " must use "
                << solidRobinTractionVelocityFvPatchVectorField::typeName
                << ", found " << patchD.type() << abort(FatalError);
        }

        solidRobinTractionVelocityFvPatchVectorField& robinPatch =
            refCast<solidRobinTractionVelocityFvPatchVectorField>(patchD);
        const vectorField solidPatchVelocity
        (
            vectorField(robinPatch)/runTime().deltaTValue()
        );
        const vectorField transferredTraction
        (
            solid().globalPatches()[interfaceI].globalFaceToPatch
            (
                solidTraction
            )
        );
        const vectorField transferredVelocity
        (
            solid().globalPatches()[interfaceI].globalFaceToPatch
            (
                solidZoneVelocity
            )
        );

        transferredTraction_[interfaceI] = transferredTraction;
        transferredVelocity_[interfaceI] = transferredVelocity;

        robinPatch.robinTraction() =
            (1 - robinDataRelaxationFactor_)*robinPatch.robinTraction()
          + robinDataRelaxationFactor_*transferredTraction;
        robinPatch.targetVelocity() =
            (1 - robinDataRelaxationFactor_)*robinPatch.targetVelocity()
          + robinDataRelaxationFactor_*transferredVelocity;

        const scalar velocityJump = gMax
        (
            mag(robinPatch.targetVelocity() - solidPatchVelocity)
        );

        Info<< couplingScheme_ << " interface " << interfaceI
            << ": Zs=" << robinPatch.impedance()
            << ", fluid force="
            << totalForceOnInterface(fluidZone, fluidTraction)
            << ", solid transmitted force="
            << totalForceOnInterface(solidZone, solidTraction)
            << ", max|u_f-v_s|=" << velocityJump << endl;
    }
}

scalar robinRobinCouplingInterface::updateDealIIResidual()
{
    scalar maximum=0;
    for (label interfaceI=0; interfaceI<nGlobalPatches(); ++interfaceI)
    {
        const vectorField& mapped=dealIIPointIncrement_[interfaceI];
        solidZonesPointsDispls()[interfaceI]=mapped;
        residualsPrev()[interfaceI]=residuals()[interfaceI];
        residuals()[interfaceI]=mapped-fluidZonesPointsDispls()[interfaceI];
        interfacesPointsDisplsPrev()[interfaceI]=interfacesPointsDispls()[interfaceI];
        interfacesPointsDispls()[interfaceI]=mapped;
        const scalar residualNorm=Foam::sqrt(gSum(magSqr(residuals()[interfaceI])));
        const scalar displacementNorm=Foam::sqrt(gSum(magSqr(mapped)));
        maxResidualsNorm()[interfaceI]=max(maxResidualsNorm()[interfaceI],residualNorm);
        maxIntsDisplsNorm()[interfaceI]=max(maxIntsDisplsNorm()[interfaceI],displacementNorm);
        const scalar relative = velocityScaleFloor_ > SMALL
          ? gMax(mag(residuals()[interfaceI]))
           /(velocityScaleFloor_*runTime().deltaTValue()+SMALL)
          : residualNorm/(maxIntsDisplsNorm()[interfaceI]+SMALL);
        maximum=max(maximum,relative);
        Info<< "deal.II FSI displacement residual " << interfaceI
            << ": absolute=" << residualNorm << ", relative=" << relative
            << endl;
    }
    return maximum;
}

scalar robinRobinCouplingInterface::reportDealIIResiduals()
{
    scalar maximumRatio=0;
    for (label interfaceI=0; interfaceI<nGlobalPatches(); ++interfaceI)
    {
        const standAlonePatch& fluidZone=
            fluid().globalPatches()[interfaceI].globalPatch();
        const vectorField fluidTraction
        (
            fluid().faceZoneViscousForce(interfaceI)
          - fluid().faceZonePressureForce(interfaceI)*fluidZone.faceNormals()
        );
        const label fluidPatchID=fluidPatchIndices()[interfaceI];
        const vectorField fluidVelocity
        (
            fluid().globalPatches()[interfaceI].patchFaceToGlobal
            (
                vectorField(fluid().U().boundaryField()[fluidPatchID])
            )
        );
        const vectorField velocityResidual
        (
            fluidVelocity-dealIIFaceVelocity_[interfaceI]
        );
        const vectorField tractionResidual
        (
            dealIIFaceTraction_[interfaceI]+fluidTraction
        );
        const vectorField robinResidual
        (
            tractionResidual
          - dealIIImpedance_*velocityResidual
        );
        const scalar velocityScale=max
        (
            velocityScaleFloor_,
            max
            (
                gMax(mag(fluidVelocity)),
                gMax(mag(dealIIFaceVelocity_[interfaceI]))
            )
        );
        const scalar tractionScale=max
        (
            tractionScaleFloor_,
            max
            (
                gMax(mag(fluidTraction)),
                gMax(mag(dealIIFaceTraction_[interfaceI]))
            )
        );
        const scalar relativeVelocity=
            gMax(mag(velocityResidual))/(velocityScale+SMALL);
        const scalar relativeTraction=
            gMax(mag(tractionResidual))/(tractionScale+SMALL);
        const scalar robinScale=
            tractionScale+dealIIImpedance_*velocityScale+SMALL;
        const scalar relativeRobin=gMax(mag(robinResidual))/robinScale;
        const scalarField powerResidual
        (
            (fluidTraction & fluidVelocity)
          + (dealIIFaceTraction_[interfaceI]
             & dealIIFaceVelocity_[interfaceI])
        );
        const scalar relativePower=gMax(mag(powerResidual))
            /(tractionScale*velocityScale+SMALL);
        const scalarField faceArea(mag(fluidZone.faceAreas()));
        const scalar integratedPowerDefect=gSum(faceArea*powerResidual);
        const scalar integratedPowerScale=max
        (
            gSum(faceArea)*tractionScale*velocityScale,
            gSum
            (
                faceArea
               *(
                    mag(fluidTraction & fluidVelocity)
                  + mag
                    (
                        dealIIFaceTraction_[interfaceI]
                      & dealIIFaceVelocity_[interfaceI]
                    )
                )
            )
        );
        const vector integratedTractionDefect=gSum
        (
            faceArea*tractionResidual
        );
        const scalar relativeIntegratedPower=
            mag(integratedPowerDefect)/(integratedPowerScale+SMALL);
        const scalar ratio=max
        (
            relativeVelocity/velocityTolerance_,
            max
            (
                relativeTraction/forceTolerance_,
                max
                (
                    relativeRobin/forceTolerance_,
                    max(relativePower,relativeIntegratedPower)/powerTolerance_
                )
            )
        );
        maximumRatio=max(maximumRatio,ratio);
        Info<< "deal.II " << couplingScheme_ << " physical residuals "
            << interfaceI
            << ": max|r_u|=" << gMax(mag(velocityResidual))
            << ", max|r_t|=" << gMax(mag(tractionResidual))
            << ", max|r_R|=" << gMax(mag(robinResidual))
            << ", relative power defect=" << relativePower
            << ", integrated traction defect="
            << integratedTractionDefect << " N"
            << ", integrated power defect=" << integratedPowerDefect << " W"
            << ", relative integrated power defect="
            << relativeIntegratedPower
            << ", acceptance ratio=" << ratio << endl;
    }
    return maximumRatio;
}

scalar robinRobinCouplingInterface::reportInterfaceResiduals()
{
    const scalar dt = runTime().deltaTValue();
    scalar maxAcceptanceRatio = 0;

    for (label interfaceI = 0; interfaceI < nGlobalPatches(); ++interfaceI)
    {
        const label solidPatchID = solidPatchIndices()[interfaceI];
        const fvPatchVectorField& patchD =
            solid().solutionD().boundaryField()[solidPatchID];

        if
        (
            patchD.type()
         != solidRobinTractionVelocityFvPatchVectorField::typeName
        )
        {
            FatalErrorInFunction
                << "Solid interface patch " << patchD.patch().name()
                << " changed type during coupling; expected "
                << solidRobinTractionVelocityFvPatchVectorField::typeName
                << ", found " << patchD.type() << abort(FatalError);
        }

        const solidRobinTractionVelocityFvPatchVectorField& robinPatch =
            refCast
            <
                const solidRobinTractionVelocityFvPatchVectorField
            >(patchD);

        const vectorField incrementalVelocity(vectorField(robinPatch)/dt);
        const vectorField solidVelocity
        (
            solid().U().boundaryField()[solidPatchID]
        );
        const vectorField& physicalFluidVelocity =
            transferredVelocity_[interfaceI];
        const vectorField& physicalSolidTractionTarget =
            transferredTraction_[interfaceI];
        const vectorField velocityResidual
        (
            physicalFluidVelocity - solidVelocity
        );
        const vectorField tractionResidual
        (
            robinPatch.traction() - physicalSolidTractionTarget
        );
        const vectorField solidRobinResidual
        (
            robinPatch.traction() - robinPatch.robinTraction()
          - robinPatch.impedance()
           *(robinPatch.targetVelocity() - solidVelocity)
        );

        // In solid-patch orientation robinTraction=-t_f.  Therefore the
        // physical interface power defect is t_s.v_s+t_f.u_f.
        const scalar powerResidual = gSum
        (
            patchD.patch().magSf()
           *(
                (robinPatch.traction() & solidVelocity)
              - (physicalSolidTractionTarget & physicalFluidVelocity)
            )
        );
        const vector forceResidual = gSum
        (
            patchD.patch().magSf()*tractionResidual
        );
        const scalar forceScale = gSum
        (
            patchD.patch().magSf()*mag(physicalSolidTractionTarget)
        );
        const scalar powerScale = gSum
        (
            patchD.patch().magSf()
           *(
                mag(robinPatch.traction() & solidVelocity)
              + mag
                (
                    physicalSolidTractionTarget
                  & physicalFluidVelocity
                )
            )
        );
        const scalar velocityScale = max
        (
            gMax(mag(physicalFluidVelocity)),
            gMax(mag(solidVelocity))
        );
        const scalar relativeVelocityDefect =
            gMax(mag(velocityResidual))/(velocityScale + SMALL);
        const scalar relativeForceDefect =
            mag(forceResidual)/(forceScale + SMALL);
        const scalar relativePowerDefect =
            mag(powerResidual)/(powerScale + SMALL);

        maxAcceptanceRatio = max
        (
            maxAcceptanceRatio,
            max
            (
                relativeVelocityDefect/velocityTolerance_,
                max
                (
                    relativeForceDefect/forceTolerance_,
                    relativePowerDefect/powerTolerance_
                )
            )
        );

        Info<< "Robin physical residuals " << interfaceI
            << ": max|r_u|=" << gMax(mag(velocityResidual))
            << ", max|U_s-DD/dt|="
            << gMax(mag(solidVelocity - incrementalVelocity))
            << ", max|r_t|=" << gMax(mag(tractionResidual))
            << ", max|r_Rs|=" << gMax(mag(solidRobinResidual))
            << ", relative velocity defect=" << relativeVelocityDefect
            << ", relative force defect="
            << relativeForceDefect
            << ", power defect=" << powerResidual << " W"
            << ", relative power defect="
            << relativePowerDefect
            << ", acceptance ratio=" << maxAcceptanceRatio << endl;
    }

    return maxAcceptanceRatio;
}

bool robinRobinCouplingInterface::evolve()
{
    initializeFields();
    updateInterpolatorAndGlobalPatches();
    scalar residualNorm = 0;
    scalar physicalAcceptanceRatio = GREAT;

    if (!coupled()) updateCoupled();

    if (predictSolid_ && coupled())
    {
        updateForce();
        if (useDealII_)
        {
            physicalAcceptanceRatio=reportDealIIResiduals();
            residualNorm=updateDealIIResidual();
        }
        else
        {
            solid().evolve();
            physicalAcceptanceRatio = reportInterfaceResiduals();
            residualNorm = updateResidual();
        }
    }

    do
    {
        outerCorr()++;
        updateDisplacement();
        moveFluidMesh();
        fluid().evolve();

        if (coupled())
        {
            updateForce();
            if (useDealII_)
            {
                physicalAcceptanceRatio=reportDealIIResiduals();
                residualNorm=updateDealIIResidual();
            }
            else
            {
                solid().evolve();
                physicalAcceptanceRatio = reportInterfaceResiduals();
                residualNorm = updateResidual();
            }
        }
        else
        {
            residualNorm = 0;
            physicalAcceptanceRatio = 0;
        }

        if (writeResidualsToFile() && Pstream::master())
        {
            residualFile() << runTime().value() << ' ' << outerCorr()
                << ' ' << residualNorm << endl;
        }
    }
    while
    (
        (residualNorm > outerCorrTolerance() || physicalAcceptanceRatio > 1)
     && outerCorr() < nOuterCorr()
    );

    if
    (
        (residualNorm > outerCorrTolerance() || physicalAcceptanceRatio > 1)
     && outerCorr() >= nOuterCorr()
    )
    {
        FatalErrorInFunction
            << couplingScheme_ << " coupling did not converge: residual="
            << residualNorm << ", tolerance=" << outerCorrTolerance()
            << ", physical acceptance ratio=" << physicalAcceptanceRatio
            << ", iterations=" << outerCorr() << abort(FatalError);
    }

    if (useDealII_)
    {
        // Persist only the state associated with the converged interface.
        // During the outer loop these fields are provisional and are never
        // promoted to accepted-state.bin.
        updatePdmsElasticWallPressure();
        if (Pstream::master() && !dealIIAcceptCommand_.empty())
        {
            const int status=std::system(dealIIAcceptCommand_.c_str());
            if (status != 0)
            {
                FatalErrorInFunction << "deal.II state promotion failed: "
                    << dealIIAcceptCommand_ << abort(FatalError);
            }
        }
        UPstream::barrier(UPstream::worldComm);
        ++acceptedStepCount_;
        if (stateAudit_ && Pstream::master())
        {
            Info<< "G1 state audit: accepted time step "
                << acceptedStepCount_ << " exactly once at t="
                << runTime().value() << endl;
        }
    }
    else
    {
        solid().updateTotalFields();
    }
    return false;
}

} // End namespace fluidSolidInterfaces
} // End namespace Foam
