#ifndef robinRobinCouplingInterface_H
#define robinRobinCouplingInterface_H

#include "fluidSolidInterface.H"

namespace Foam
{
namespace fluidSolidInterfaces
{

class robinRobinCouplingInterface : public fluidSolidInterface
{
    const scalar relaxationFactor_;
    const word couplingScheme_;
    const scalar semiImplicitBeta_;
    const word accelerationMethod_;
    const scalar relaxationFactorMin_;
    const scalar relaxationFactorMax_;
    List<scalar> aitkenRelaxationFactors_;
    const label iqnMaxModes_;
    const scalar iqnMaxCorrectionFactor_;
    const bool iqnReuseAcrossTimeSteps_;
    const scalar iqnModeFilterTolerance_;
    List<DynamicList<vectorField>> iqnV_;
    List<DynamicList<vectorField>> iqnW_;
    List<label> iqnCurrentStepModes_;
    List<DynamicList<vectorField>> iqnDataV_;
    List<DynamicList<vectorField>> iqnDataW_;
    List<vectorField> iqnDataResidualPrev_;
    List<vectorField> iqnDataOutputPrev_;
    const scalar robinDataRelaxationFactor_;
    const bool predictSolid_;
    const scalar velocityTolerance_;
    const scalar forceTolerance_;
    const scalar powerTolerance_;
    const scalar velocityScaleFloor_;
    const scalar tractionScaleFloor_;
    List<vectorField> transferredTraction_;
    List<vectorField> transferredVelocity_;
    const bool useDealII_;
    const string dealIICommand_;
    const string dealIIAcceptCommand_;
    const fileName dealIIInput_;
    const fileName dealIIQueryInput_;
    const fileName dealIIOutput_;
    const scalar dealIIImpedance_;
    const bool stateAudit_;
    label acceptedStepCount_;
    List<vectorField> dealIIPointIncrement_;
    List<vectorField> dealIIFaceVelocity_;
    List<vectorField> dealIIFaceTraction_;
    List<vectorField> dealIIFaceAcceleration_;

    robinRobinCouplingInterface(const robinRobinCouplingInterface&);
    void operator=(const robinRobinCouplingInterface&);

public:
    TypeName("robinRobin");

    robinRobinCouplingInterface
    (
        Time& runTime,
        const word& region = dynamicFvMesh::defaultRegion
    );

    virtual ~robinRobinCouplingInterface() {}
    virtual scalar relaxationFactor() const { return relaxationFactor_; }
    virtual bool evolve();
    virtual void updateDisplacement();
    virtual void updateForce();
    void updatePdmsElasticWallPressure();
    scalar reportInterfaceResiduals();
    scalar updateDealIIResidual();
    scalar reportDealIIResiduals();
    void solveDealIISolid
    (
        const label interfaceI,
        const standAlonePatch& fluidZone,
        const vectorField& fluidTraction,
        const vectorField& fluidVelocity
    );
};

} // End namespace fluidSolidInterfaces
} // End namespace Foam

#endif
