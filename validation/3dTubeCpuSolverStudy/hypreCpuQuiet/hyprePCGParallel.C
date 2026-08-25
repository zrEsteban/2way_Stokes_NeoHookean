#include "hyprePCGParallel.H"

#include "HYPRE.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_config.h"
#include "Pstream.H"
#include "processorLduInterface.H"
#include "cyclicLduInterface.H"
#include "clockTime.H"
#include "regIOobject.H"
#include "objectRegistry.H"
#include "IOdictionary.H"
#include "Time.H"
#include "mpi.h"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <map>
#include <memory>
#include <sstream>
#include <vector>

namespace Foam
{
    defineTypeNameAndDebug(hyprePCGParallel, 0);
    lduMatrix::solver::addsymMatrixConstructorToTable<hyprePCGParallel>
        addHyprePCGSymMatrixConstructorToTable_;
    defineTypeNameAndDebug(hypreFlexGMRESParallel, 0);
    lduMatrix::solver::addsymMatrixConstructorToTable<hypreFlexGMRESParallel>
        addHypreFlexGMRESSymMatrixConstructorToTable_;
    lduMatrix::solver::addasymMatrixConstructorToTable<hypreFlexGMRESParallel>
        addHypreFlexGMRESAsymMatrixConstructorToTable_;
}

namespace
{

using Foam::Info;
using Foam::FatalError;
using Foam::label;
using Foam::nl;
using Foam::scalar;
using Foam::word;

void checkHypre(const HYPRE_Int code, const char* operation)
{
    if (code)
    {
        FatalErrorInFunction
            << "HYPRE operation failed: " << operation
            << " (error code " << code << ')'
            << Foam::exit(FatalError);
    }
}

class HypreRuntime
{
    bool ownsMpi_;
    bool active_;
    label references_;

public:
    HypreRuntime()
    :
        ownsMpi_(false),
        active_(false),
        references_(0)
    {}

    void acquire()
    {
        ++references_;
        if (active_) return;

        int initialized = 0;
        if (MPI_Initialized(&initialized) != MPI_SUCCESS)
        {
            FatalErrorInFunction << "MPI_Initialized failed"
                << Foam::exit(FatalError);
        }
        if (!initialized)
        {
            if (MPI_Init(nullptr, nullptr) != MPI_SUCCESS)
            {
                FatalErrorInFunction << "MPI_Init failed"
                    << Foam::exit(FatalError);
            }
            ownsMpi_ = true;
        }

        if (!HYPRE_Initialized())
        {
            checkHypre(HYPRE_Initialize(), "HYPRE_Initialize");
        }
#if defined(HYPRE_USING_CUDA)
        // Phase 1 targets distributed-memory CPU execution. Sharing a single
        // CUDA device between many MPI ranks is neither scalable nor a valid
        // substitute for rank-to-device binding on a multi-GPU system.
        checkHypre(HYPRE_SetMemoryLocation(HYPRE_MEMORY_HOST),
            "HYPRE_SetMemoryLocation(host)");
        checkHypre(HYPRE_SetExecutionPolicy(HYPRE_EXEC_HOST),
            "HYPRE_SetExecutionPolicy(host)");
#if defined(HYPRE_USING_CURAND)
        checkHypre(HYPRE_SetUseGpuRand(1), "HYPRE_SetUseGpuRand");
#endif
#endif
        active_ = true;
        Info<< "HYPRE runtime initialized" << Foam::endl;
    }

    void release()
    {
        if (references_ > 0) --references_;
    }

    // OpenFOAM can finalize MPI explicitly before static library objects are
    // destroyed. At process shutdown, avoid HYPRE cleanup routines that call
    // MPI_Comm_free after MPI_Finalize; the OS will reclaim these resources.
    void abandon()
    {
        active_ = false;
        ownsMpi_ = false;
        references_ = 0;
    }

    void finalize()
    {
        if (active_ && !HYPRE_Finalized())
        {
            checkHypre(HYPRE_Finalize(), "HYPRE_Finalize");
            Info<< "HYPRE runtime finalized" << Foam::endl;
        }
        if (ownsMpi_)
        {
            int finalized = 0;
            if (MPI_Finalized(&finalized) == MPI_SUCCESS && !finalized)
            {
                MPI_Finalize();
            }
        }
        active_ = false;
    }

    ~HypreRuntime()
    {
        finalize();
    }
};

enum class ReuseMode { none, matrix, preconditioner };

const char* modeName(const ReuseMode mode)
{
    switch (mode)
    {
        case ReuseMode::none: return "none";
        case ReuseMode::matrix: return "matrix";
        default: return "preconditioner";
    }
}

struct Controls
{
    ReuseMode mode;
    scalar matrixTolerance;
    bool verifyCoefficients;
    bool strictReuse;
    bool reportTimings;
    bool reportLifecycle;
    bool reportPerSolve;
    bool reportSummary;
    label printLevel;
    label logging;
    label coarsenType;
    label interpType;
    label relaxType;
    label relaxOrder;
    label numSweeps;
    label maxLevels;
    label kDim;
    scalar strongThreshold;
};

std::uint64_t fnvBytes
(
    std::uint64_t hash,
    const void* data,
    const std::size_t size
)
{
    const unsigned char* bytes = static_cast<const unsigned char*>(data);
    for (std::size_t i = 0; i < size; ++i)
    {
        hash ^= std::uint64_t(bytes[i]);
        hash *= UINT64_C(1099511628211);
    }
    return hash;
}

template<class List>
std::uint64_t hashList(std::uint64_t hash, const List& list)
{
    for (const auto& value : list)
    {
        hash = fnvBytes(hash, &value, sizeof(value));
    }
    return hash;
}

struct MatrixSignature
{
    std::uint64_t structure;
    std::uint64_t coefficients;
};

MatrixSignature signature
(
    const Foam::lduMatrix& matrix,
    const Foam::FieldField<Foam::Field, scalar>& interfaceBouCoeffs,
    const Foam::lduInterfaceFieldPtrsList& interfaces
)
{
    std::uint64_t structure = UINT64_C(14695981039346656037);
    std::uint64_t coefficients = structure;
    const label nRows = matrix.diag().size();
    const label nFaces = matrix.upper().size();
    const bool symmetric = matrix.symmetric();
    structure = fnvBytes(structure, &nRows, sizeof(nRows));
    structure = fnvBytes(structure, &nFaces, sizeof(nFaces));
    structure = fnvBytes(structure, &symmetric, sizeof(symmetric));
    structure = hashList(structure, matrix.lduAddr().lowerAddr());
    structure = hashList(structure, matrix.lduAddr().upperAddr());
    coefficients = hashList(coefficients, matrix.diag());
    coefficients = hashList(coefficients, matrix.lower());
    coefficients = hashList(coefficients, matrix.upper());
    forAll(interfaces, patchi)
    {
        if (interfaces.set(patchi))
        {
            coefficients =
                hashList(coefficients, interfaceBouCoeffs[patchi]);
        }
    }
    return {structure, coefficients};
}

struct Timings
{
    scalar matrixInspection = 0;
    scalar contextCreation = 0;
    scalar matrixValueUpdate = 0;
    scalar matrixAssembly = 0;
    scalar rhsUpdate = 0;
    scalar solutionUpdate = 0;
    scalar amgSetup = 0;
    scalar pcgSolve = 0;
    scalar solutionCopy = 0;
    scalar resourceCleanup = 0;
    scalar totalHypreSolve = 0;
};

struct HypreSolverContext
{
    HYPRE_IJMatrix ijMatrix = nullptr;
    HYPRE_ParCSRMatrix parMatrix = nullptr;
    HYPRE_IJVector ijRhs = nullptr;
    HYPRE_IJVector ijSolution = nullptr;
    HYPRE_ParVector parRhs = nullptr;
    HYPRE_ParVector parSolution = nullptr;
    HYPRE_Solver pcg = nullptr;
    HYPRE_Solver amg = nullptr;
    bool flexGmres = false;

    label nRows = 0;
    HYPRE_BigInt rowStart = 0;
    HYPRE_BigInt rowEnd = -1;
    HYPRE_BigInt globalRows = 0;
    label nInternalFaces = 0;
    label matrixRevision = 0;
    bool matrixValid = false;
    bool vectorsValid = false;
    bool amgSetupValid = false;
    scalar matrixScale = 1;
    std::uint64_t structureHash = 0;
    std::uint64_t coefficientHash = 0;

    label nSolveCalls = 0;
    label nMatrixAssemblies = 0;
    label nAMGSetups = 0;
    label nAMGReuses = 0;
    label nContextRebuilds = 0;
    label nTopologyRebuilds = 0;
    label nStructuralRebuilds = 0;
    label lastMeshTopologyEvent = 0;
    scalar cumulativeSetupTime = 0;
    scalar cumulativeSolveTime = 0;

    std::vector<HYPRE_Int> rowSizes;
    std::vector<HYPRE_BigInt> rows;
    std::vector<HYPRE_BigInt> columns;
    std::vector<HYPRE_BigInt> indices;
    std::vector<HYPRE_Complex> values;
    std::vector<HYPRE_Complex> rhsValues;
    std::vector<HYPRE_Complex> solutionValues;
    std::vector<scalar> coefficientSnapshot;
    std::vector<std::vector<HYPRE_BigInt>> interfaceColumns;

    void destroy()
    {
        if (pcg)
        {
            if (flexGmres)
            {
                checkHypre
                (
                    HYPRE_ParCSRFlexGMRESDestroy(pcg),
                    "FlexGMRESDestroy"
                );
            }
            else
            {
                checkHypre(HYPRE_ParCSRPCGDestroy(pcg), "PCGDestroy");
            }
        }
        if (amg) checkHypre(HYPRE_BoomerAMGDestroy(amg), "AMGDestroy");
        if (ijMatrix) checkHypre(HYPRE_IJMatrixDestroy(ijMatrix), "IJMatrixDestroy");
        if (ijRhs) checkHypre(HYPRE_IJVectorDestroy(ijRhs), "IJVectorDestroy(rhs)");
        if (ijSolution) checkHypre(HYPRE_IJVectorDestroy(ijSolution), "IJVectorDestroy(solution)");
        pcg = nullptr; amg = nullptr; ijMatrix = nullptr;
        flexGmres = false;
        ijRhs = nullptr; ijSolution = nullptr;
        parMatrix = nullptr; parRhs = nullptr; parSolution = nullptr;
        matrixValid = false; vectorsValid = false; amgSetupValid = false;
    }

    ~HypreSolverContext() { destroy(); }
};

struct Totals
{
    label solves = 0;
    label rebuilds = 0;
    label assemblies = 0;
    label setups = 0;
    label reuses = 0;
    label topologyRebuilds = 0;
    label structuralRebuilds = 0;
    label krylovIterations = 0;
    scalar setupTime = 0;
    scalar solveTime = 0;
    scalar inspectionTime = 0;
    scalar contextTime = 0;
    scalar matrixUpdateTime = 0;
    scalar matrixAssemblyTime = 0;
    scalar vectorUpdateTime = 0;
    scalar solutionCopyTime = 0;
    scalar totalTime = 0;
};

class ContextRegistry
{
    HypreRuntime runtime_;
    std::map<std::string, std::unique_ptr<HypreSolverContext>> contexts_;
    Totals totals_;
    bool report_;
    bool shutDown_;

public:
    ContextRegistry() : report_(false), shutDown_(false) { runtime_.acquire(); }

    HypreSolverContext& get(const std::string& key)
    {
        auto& ptr = contexts_[key];
        if (!ptr) ptr.reset(new HypreSolverContext);
        return *ptr;
    }

    void account
    (
        const label rebuilds,
        const label assemblies,
        const label setups,
        const label reuses,
        const label topologyRebuilds,
        const label structuralRebuilds,
        const label krylovIterations,
        const scalar setupTime,
        const scalar solveTime,
        const Timings& timing,
        const bool report
    )
    {
        ++totals_.solves;
        totals_.rebuilds += rebuilds;
        totals_.assemblies += assemblies;
        totals_.setups += setups;
        totals_.reuses += reuses;
        totals_.topologyRebuilds += topologyRebuilds;
        totals_.structuralRebuilds += structuralRebuilds;
        totals_.krylovIterations += krylovIterations;
        totals_.setupTime += setupTime;
        totals_.solveTime += solveTime;
        totals_.inspectionTime += timing.matrixInspection;
        totals_.contextTime += timing.contextCreation;
        totals_.matrixUpdateTime += timing.matrixValueUpdate;
        totals_.matrixAssemblyTime += timing.matrixAssembly;
        totals_.vectorUpdateTime += timing.rhsUpdate + timing.solutionUpdate;
        totals_.solutionCopyTime += timing.solutionCopy;
        totals_.totalTime += timing.totalHypreSolve;
        report_ = report_ || report;
    }

    void shutdown()
    {
        if (shutDown_) return;
        shutDown_ = true;
        if (report_)
        {
            Info<< nl << "HYPRE lifecycle summary" << nl
                << "    solve calls:              " << totals_.solves << nl
                << "    context rebuilds:         " << totals_.rebuilds << nl
                << "    matrix assemblies:        " << totals_.assemblies << nl
                << "    AMG setups:               " << totals_.setups << nl
                << "    AMG reuses:               " << totals_.reuses << nl
                << "    topology rebuilds:        "
                << totals_.topologyRebuilds << nl
                << "    structural rebuilds:      "
                << totals_.structuralRebuilds << nl
                << "    Krylov iterations:        "
                << totals_.krylovIterations << nl
                << "    cumulative setup time:    " << totals_.setupTime << nl
                << "    cumulative solve time:    " << totals_.solveTime << nl
                << "    matrix inspection time:   " << totals_.inspectionTime << nl
                << "    context creation time:    " << totals_.contextTime << nl
                << "    matrix value update time: " << totals_.matrixUpdateTime << nl
                << "    matrix assembly time:     " << totals_.matrixAssemblyTime << nl
                << "    vector update time:       " << totals_.vectorUpdateTime << nl
                << "    solution copy time:       " << totals_.solutionCopyTime << nl
                << "    total wrapper time:       " << totals_.totalTime << nl
                << "    average solve time:       "
                << (totals_.solves ? totals_.solveTime/totals_.solves : 0)
                << Foam::endl;
        }

        int mpiFinalized = 0;
        if
        (
            MPI_Finalized(&mpiFinalized) == MPI_SUCCESS
         && mpiFinalized
        )
        {
            // HYPRE object destructors call MPI_Comm_free. Since OpenFOAM has
            // already finalized MPI, intentionally release ownership at final
            // process teardown instead of making an invalid MPI call.
            for (auto& item : contexts_)
            {
                item.second.release();
            }
            contexts_.clear();
            runtime_.abandon();
            return;
        }

        // HYPRE objects must die before HYPRE_Finalize.
        contexts_.clear();
        runtime_.release();
        runtime_.finalize();
    }

    ~ContextRegistry() { shutdown(); }
};

ContextRegistry& registry()
{
    static ContextRegistry instance;
    return instance;
}

class HypreMeshLifecycle : public Foam::regIOobject
{
public:
    TypeName("hypreMeshLifecycle");

    explicit HypreMeshLifecycle(const Foam::objectRegistry& db)
    :
        regIOobject
        (
            Foam::IOobject
            (
                "hypreMeshLifecycle",
                db.time().timeName(),
                db,
                Foam::IOobject::NO_READ,
                Foam::IOobject::NO_WRITE,
                true
            )
        )
    {}

    ~HypreMeshLifecycle() override
    {
        registry().shutdown();
    }

    bool writeData(Foam::Ostream&) const override { return true; }
};

defineTypeNameAndDebug(HypreMeshLifecycle, 0);

void ensureMeshLifecycle(const Foam::lduMesh& mesh)
{
    const Foam::objectRegistry& db = mesh.thisDb();
    if (!db.foundObject<HypreMeshLifecycle>("hypreMeshLifecycle"))
    {
        Foam::regIOobject::store(new HypreMeshLifecycle(db));
    }
}

Controls parseControls(const Foam::dictionary& solverDict)
{
    const Foam::dictionary& dict = solverDict.optionalSubDict("hypre");
    const Foam::dictionary& reuse = dict.optionalSubDict("reuse");
    const word mode(reuse.lookupOrDefault<word>("mode", "none"));
    ReuseMode parsed;
    if (mode == "none") parsed = ReuseMode::none;
    else if (mode == "matrix") parsed = ReuseMode::matrix;
    else if (mode == "preconditioner") parsed = ReuseMode::preconditioner;
    else
    {
        FatalIOErrorInFunction(reuse)
            << "Unknown hypre reuse mode '" << mode
            << "'. Valid modes are: none, matrix, preconditioner"
            << Foam::exit(Foam::FatalIOError);
        parsed = ReuseMode::none;
    }
    Controls c
    {
        parsed,
        reuse.lookupOrDefault<scalar>("matrixTolerance", 0),
        reuse.lookupOrDefault<bool>("verifyCoefficients", true),
        reuse.lookupOrDefault<bool>("strictReuse", false),
        dict.lookupOrDefault<bool>("reportTimings", true),
        dict.lookupOrDefault<bool>("reportLifecycle", true),
        dict.lookupOrDefault<bool>("reportPerSolve", false),
        dict.lookupOrDefault<bool>("reportSummary", true),
        dict.lookupOrDefault<label>("printLevel", 1),
        dict.lookupOrDefault<label>("logging", 1),
        dict.lookupOrDefault<label>("coarsenType", 8),
        dict.lookupOrDefault<label>("interpType", 6),
        dict.lookupOrDefault<label>("relaxType", 18),
        dict.lookupOrDefault<label>("relaxOrder", 0),
        dict.lookupOrDefault<label>("numSweeps", 1),
        dict.lookupOrDefault<label>("maxLevels", 25),
        dict.lookupOrDefault<label>("kDim", 30),
        dict.lookupOrDefault<scalar>("strongThreshold", 0.25)
    };
    if (c.matrixTolerance < 0 || c.printLevel < 0 || c.logging < 0
     || c.coarsenType < 0 || c.interpType < 0 || c.relaxType < 0
     || (c.relaxOrder != 0 && c.relaxOrder != 1)
     || c.numSweeps < 1 || c.maxLevels < 1 || c.kDim < 1
     || c.strongThreshold <= 0 || c.strongThreshold > 1)
    {
        FatalIOErrorInFunction(dict) << "Invalid hypre controls"
            << Foam::exit(Foam::FatalIOError);
    }
    return c;
}

void snapshotCoefficients
(
    HypreSolverContext& ctx,
    const Foam::lduMatrix& matrix,
    const Foam::FieldField<Foam::Field, scalar>& interfaceBouCoeffs,
    const Foam::lduInterfaceFieldPtrsList& interfaces
)
{
    ctx.coefficientSnapshot.clear();
    ctx.coefficientSnapshot.reserve
    (
        matrix.diag().size() + matrix.lower().size() + matrix.upper().size()
    );
    ctx.coefficientSnapshot.insert(ctx.coefficientSnapshot.end(),
        matrix.diag().begin(), matrix.diag().end());
    ctx.coefficientSnapshot.insert(ctx.coefficientSnapshot.end(),
        matrix.lower().begin(), matrix.lower().end());
    ctx.coefficientSnapshot.insert(ctx.coefficientSnapshot.end(),
        matrix.upper().begin(), matrix.upper().end());
    forAll(interfaces, patchi)
    {
        if (interfaces.set(patchi))
        {
            ctx.coefficientSnapshot.insert
            (
                ctx.coefficientSnapshot.end(),
                interfaceBouCoeffs[patchi].begin(),
                interfaceBouCoeffs[patchi].end()
            );
        }
    }
}

bool coefficientsEquivalent
(
    const HypreSolverContext& ctx,
    const Foam::lduMatrix& matrix,
    const Foam::FieldField<Foam::Field, scalar>& interfaceBouCoeffs,
    const Foam::lduInterfaceFieldPtrsList& interfaces,
    const Controls& controls,
    const std::uint64_t newHash
)
{
    if (ctx.coefficientHash == newHash) return true;
    if (!controls.verifyCoefficients) return false;
    if (controls.matrixTolerance <= 0) return false;
    std::size_t i = 0;
    const auto compare = [&](const auto& values)
    {
        for (const scalar value : values)
        {
            if (i >= ctx.coefficientSnapshot.size()) return false;
            const scalar old = ctx.coefficientSnapshot[i++];
            const scalar scale = std::max(scalar(1), std::max(Foam::mag(old), Foam::mag(value)));
            if (Foam::mag(value - old) > controls.matrixTolerance*scale) return false;
        }
        return true;
    };
    if (!(compare(matrix.diag()) && compare(matrix.lower())
       && compare(matrix.upper())))
    {
        return false;
    }
    forAll(interfaces, patchi)
    {
        if (interfaces.set(patchi) && !compare(interfaceBouCoeffs[patchi]))
        {
            return false;
        }
    }
    return i == ctx.coefficientSnapshot.size();
}

void buildSparsity
(
    HypreSolverContext& ctx,
    const Foam::lduMatrix& matrix,
    const Foam::lduInterfaceFieldPtrsList& interfaces
)
{
    const label nRows = matrix.diag().size();
    const label nFaces = matrix.upper().size();
    const Foam::labelUList& owner = matrix.lduAddr().lowerAddr();
    const Foam::labelUList& neighbour = matrix.lduAddr().upperAddr();
    long long localRows = nRows;
    long long rowStart = 0;
    long long globalRows = 0;
    MPI_Exscan
    (
        &localRows,
        &rowStart,
        1,
        MPI_LONG_LONG,
        MPI_SUM,
        MPI_COMM_WORLD
    );
    int mpiRank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    if (mpiRank == 0) rowStart = 0;
    MPI_Allreduce
    (
        &localRows,
        &globalRows,
        1,
        MPI_LONG_LONG,
        MPI_SUM,
        MPI_COMM_WORLD
    );
    ctx.rowStart = HYPRE_BigInt(rowStart);
    ctx.rowEnd = ctx.rowStart + nRows - 1;
    ctx.globalRows = HYPRE_BigInt(globalRows);

    Foam::labelList globalCellIds(nRows);
    for (label row = 0; row < nRows; ++row)
    {
        globalCellIds[row] = label(ctx.rowStart + row);
    }

    ctx.interfaceColumns.clear();
    ctx.interfaceColumns.resize(interfaces.size());
    const label transferRequestStart = Foam::UPstream::nRequests();
    forAll(interfaces, patchi)
    {
        if (!interfaces.set(patchi)) continue;
        const Foam::lduInterface& interface = interfaces[patchi].interface();
        if
        (
            !Foam::isA<Foam::processorLduInterface>(interface)
         && !Foam::isA<Foam::cyclicLduInterface>(interface)
        )
        {
            FatalErrorInFunction
                << "hyprePCGParallel supports processor and conformal cyclic "
                << "interfaces; patch " << patchi << " has type "
                << interface.type() << Foam::exit(FatalError);
        }
        interface.initInternalFieldTransfer
        (
            Foam::Pstream::commsTypes::nonBlocking,
            globalCellIds
        );
    }
    Foam::UPstream::waitRequests(transferRequestStart);
    forAll(interfaces, patchi)
    {
        if (!interfaces.set(patchi)) continue;
        const Foam::lduInterface& interface = interfaces[patchi].interface();
        Foam::tmp<Foam::labelField> transferred =
            interface.internalFieldTransfer
            (
                Foam::Pstream::commsTypes::nonBlocking,
                globalCellIds
            );
        const Foam::labelField& neighbourIds = transferred();
        auto& remote = ctx.interfaceColumns[patchi];
        remote.resize(neighbourIds.size());
        forAll(neighbourIds, facei)
        {
            remote[facei] = HYPRE_BigInt(neighbourIds[facei]);
        }
    }

    ctx.rowSizes.assign(nRows, 1);
    for (label face = 0; face < nFaces; ++face)
    {
        ++ctx.rowSizes[owner[face]];
        ++ctx.rowSizes[neighbour[face]];
    }
    forAll(interfaces, patchi)
    {
        if (!interfaces.set(patchi)) continue;
        const Foam::labelUList& faceCells =
            interfaces[patchi].interface().faceCells();
        forAll(faceCells, facei)
        {
            ++ctx.rowSizes[faceCells[facei]];
        }
    }
    std::vector<std::size_t> offsets(nRows + 1, 0);
    for (label row = 0; row < nRows; ++row)
        offsets[row + 1] = offsets[row] + ctx.rowSizes[row];
    ctx.rows.resize(nRows);
    ctx.indices.resize(nRows);
    ctx.columns.resize(offsets.back());
    ctx.values.resize(offsets.back());
    std::vector<std::size_t> next(offsets.begin(), offsets.end() - 1);
    for (label row = 0; row < nRows; ++row)
    {
        ctx.rows[row] = ctx.rowStart + row;
        ctx.indices[row] = ctx.rowStart + row;
        ctx.columns[next[row]++] = ctx.rowStart + row;
    }
    for (label face = 0; face < nFaces; ++face)
    {
        ctx.columns[next[owner[face]]++] = ctx.rowStart + neighbour[face];
        ctx.columns[next[neighbour[face]]++] = ctx.rowStart + owner[face];
    }
    forAll(interfaces, patchi)
    {
        if (!interfaces.set(patchi)) continue;
        const Foam::labelUList& faceCells =
            interfaces[patchi].interface().faceCells();
        const auto& remote = ctx.interfaceColumns[patchi];
        forAll(faceCells, facei)
        {
            ctx.columns[next[faceCells[facei]]++] = remote[facei];
        }
    }
    HYPRE_BigInt minColumn = ctx.globalRows;
    HYPRE_BigInt maxColumn = -1;
    for (const HYPRE_BigInt column : ctx.columns)
    {
        minColumn = std::min(minColumn, column);
        maxColumn = std::max(maxColumn, column);
        if (column < 0 || column >= ctx.globalRows)
        {
            FatalErrorInFunction
                << "invalid global column " << column
                << " outside [0," << ctx.globalRows - 1 << ']'
                << Foam::exit(FatalError);
        }
    }
    for (label row = 0; row < nRows; ++row)
    {
        if (next[row] != offsets[row + 1])
        {
            FatalErrorInFunction
                << "row fill mismatch at local row " << row
                << Foam::exit(FatalError);
        }
    }
    ctx.rhsValues.resize(nRows);
    ctx.solutionValues.resize(nRows);
}

void fillMatrixValues
(
    HypreSolverContext& ctx,
    const Foam::lduMatrix& matrix,
    const Foam::FieldField<Foam::Field, scalar>& interfaceBouCoeffs,
    const Foam::lduInterfaceFieldPtrsList& interfaces
)
{
    const label nRows = matrix.diag().size();
    const label nFaces = matrix.upper().size();
    const Foam::labelUList& owner = matrix.lduAddr().lowerAddr();
    const Foam::labelUList& neighbour = matrix.lduAddr().upperAddr();
    std::vector<std::size_t> offsets(nRows + 1, 0);
    for (label row = 0; row < nRows; ++row)
        offsets[row + 1] = offsets[row] + ctx.rowSizes[row];
    std::vector<std::size_t> next(offsets.begin(), offsets.end() - 1);
    std::vector<Foam::scalarField> neighbourCoeffs(interfaces.size());
    const label coefficientRequestStart = Foam::UPstream::nRequests();
    forAll(interfaces, patchi)
    {
        if (!interfaces.set(patchi)) continue;
        const Foam::lduInterface& interface=interfaces[patchi].interface();
        if(Foam::isA<Foam::processorLduInterface>(interface))
        {
            Foam::refCast<const Foam::processorLduInterface>(interface).send
            (
                Foam::Pstream::commsTypes::nonBlocking,
                interfaceBouCoeffs[patchi]
            );
        }
    }
    Foam::UPstream::waitRequests(coefficientRequestStart);
    scalar localMaxInterfaceAsymmetry = 0;
    forAll(interfaces, patchi)
    {
        if (!interfaces.set(patchi)) continue;
        const Foam::lduInterface& interface=interfaces[patchi].interface();
        neighbourCoeffs[patchi].resize(interfaceBouCoeffs[patchi].size());
        if(Foam::isA<Foam::processorLduInterface>(interface))
        {
            Foam::refCast<const Foam::processorLduInterface>(interface).receive
            (
                Foam::Pstream::commsTypes::nonBlocking,
                neighbourCoeffs[patchi]
            );
        }
        else
        {
            const label neighbourPatch=
                Foam::refCast<const Foam::cyclicLduInterface>(interface)
               .neighbPatchID();
            neighbourCoeffs[patchi]=interfaceBouCoeffs[neighbourPatch];
        }
        forAll(interfaceBouCoeffs[patchi], facei)
        {
            const scalar local = interfaceBouCoeffs[patchi][facei];
            const scalar remote = neighbourCoeffs[patchi][facei];
            localMaxInterfaceAsymmetry = std::max
            (
                localMaxInterfaceAsymmetry,
                Foam::mag(local - remote)
               /std::max(scalar(1), std::max(Foam::mag(local), Foam::mag(remote)))
            );
        }
    }
    scalar globalMaxInterfaceAsymmetry = 0;
    MPI_Allreduce
    (
        &localMaxInterfaceAsymmetry,
        &globalMaxInterfaceAsymmetry,
        1,
        MPI_DOUBLE,
        MPI_MAX,
        MPI_COMM_WORLD
    );
    if (matrix.symmetric() && globalMaxInterfaceAsymmetry > 1e-12)
    {
        FatalErrorInFunction
            << "distributed matrix is not symmetric across processor "
            << "interfaces; max relative mismatch="
            << globalMaxInterfaceAsymmetry << Foam::exit(FatalError);
    }
    for (label row = 0; row < nRows; ++row)
        ctx.values[next[row]++] =
            HYPRE_Complex(ctx.matrixScale*matrix.diag()[row]);
    for (label face = 0; face < nFaces; ++face)
    {
        ctx.values[next[owner[face]]++] =
            HYPRE_Complex(ctx.matrixScale*matrix.upper()[face]);
        ctx.values[next[neighbour[face]]++] =
            HYPRE_Complex(ctx.matrixScale*matrix.lower()[face]);
    }
    forAll(interfaces, patchi)
    {
        if (!interfaces.set(patchi)) continue;
        const Foam::labelUList& faceCells =
            interfaces[patchi].interface().faceCells();
        const Foam::scalarField& coeffs = interfaceBouCoeffs[patchi];
        if (coeffs.size() != faceCells.size())
        {
            FatalErrorInFunction
                << "processor interface coefficient/addressing size mismatch"
                << Foam::exit(FatalError);
        }
        forAll(faceCells, facei)
        {
            ctx.values[next[faceCells[facei]]++] =
                HYPRE_Complex(-ctx.matrixScale*coeffs[facei]);
        }
    }
}

void createContext
(
    HypreSolverContext& ctx,
    const Foam::lduMatrix& matrix,
    const Foam::lduInterfaceFieldPtrsList& interfaces,
    const Controls& c,
    const scalar tolerance,
    const label maxIter,
    const bool flexGmres
)
{
    const label nRows = matrix.diag().size();
    buildSparsity(ctx, matrix, interfaces);
    checkHypre(HYPRE_IJMatrixCreate
    (
        MPI_COMM_WORLD,
        ctx.rowStart,
        ctx.rowEnd,
        ctx.rowStart,
        ctx.rowEnd,
        &ctx.ijMatrix), "IJMatrixCreate");
    checkHypre(HYPRE_IJMatrixSetObjectType(ctx.ijMatrix, HYPRE_PARCSR),
        "IJMatrixSetObjectType");
    checkHypre(HYPRE_IJMatrixSetRowSizes(ctx.ijMatrix, ctx.rowSizes.data()),
        "IJMatrixSetRowSizes");

    checkHypre(HYPRE_IJVectorCreate
    (
        MPI_COMM_WORLD,
        ctx.rowStart,
        ctx.rowEnd,
        &ctx.ijRhs
    ),
        "IJVectorCreate(rhs)");
    checkHypre(HYPRE_IJVectorSetObjectType(ctx.ijRhs, HYPRE_PARCSR),
        "IJVectorSetObjectType(rhs)");
    checkHypre(HYPRE_IJVectorCreate
    (
        MPI_COMM_WORLD,
        ctx.rowStart,
        ctx.rowEnd,
        &ctx.ijSolution
    ),
        "IJVectorCreate(solution)");
    checkHypre(HYPRE_IJVectorSetObjectType(ctx.ijSolution, HYPRE_PARCSR),
        "IJVectorSetObjectType(solution)");

    checkHypre(HYPRE_BoomerAMGCreate(&ctx.amg), "BoomerAMGCreate");
    checkHypre(HYPRE_BoomerAMGSetCoarsenType(ctx.amg, c.coarsenType), "AMG coarsenType");
    checkHypre(HYPRE_BoomerAMGSetInterpType(ctx.amg, c.interpType), "AMG interpType");
    checkHypre(HYPRE_BoomerAMGSetRelaxType(ctx.amg, c.relaxType), "AMG relaxType");
    checkHypre(HYPRE_BoomerAMGSetRelaxOrder(ctx.amg, c.relaxOrder), "AMG relaxOrder");
    checkHypre(HYPRE_BoomerAMGSetNumSweeps(ctx.amg, c.numSweeps), "AMG numSweeps");
    checkHypre(HYPRE_BoomerAMGSetMaxLevels(ctx.amg, c.maxLevels), "AMG maxLevels");
    checkHypre(HYPRE_BoomerAMGSetStrongThreshold(ctx.amg, c.strongThreshold), "AMG threshold");
    checkHypre(HYPRE_BoomerAMGSetMaxIter(ctx.amg, 1), "AMG maxIter");
    checkHypre(HYPRE_BoomerAMGSetTol(ctx.amg, 0), "AMG tol");

    ctx.flexGmres = flexGmres;
    if (flexGmres)
    {
        checkHypre
        (
            HYPRE_ParCSRFlexGMRESCreate(MPI_COMM_WORLD, &ctx.pcg),
            "FlexGMRESCreate"
        );
        checkHypre(HYPRE_ParCSRFlexGMRESSetTol(ctx.pcg, tolerance),
            "FlexGMRES tol");
        checkHypre(HYPRE_ParCSRFlexGMRESSetAbsoluteTol(ctx.pcg, 0),
            "FlexGMRES abs tol");
        checkHypre(HYPRE_ParCSRFlexGMRESSetMaxIter(ctx.pcg, maxIter),
            "FlexGMRES maxIter");
        checkHypre(HYPRE_ParCSRFlexGMRESSetKDim(ctx.pcg, c.kDim),
            "FlexGMRES kDim");
        checkHypre(HYPRE_ParCSRFlexGMRESSetPrintLevel(ctx.pcg, c.printLevel),
            "FlexGMRES print");
        checkHypre(HYPRE_ParCSRFlexGMRESSetLogging(ctx.pcg, c.logging),
            "FlexGMRES logging");
        checkHypre
        (
            HYPRE_ParCSRFlexGMRESSetPrecond
            (
                ctx.pcg,
                HYPRE_BoomerAMGSolve,
                HYPRE_BoomerAMGSetup,
                ctx.amg
            ),
            "FlexGMRESSetPrecond"
        );
    }
    else
    {
        checkHypre(HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &ctx.pcg), "PCGCreate");
        checkHypre(HYPRE_ParCSRPCGSetTol(ctx.pcg, tolerance), "PCG tol");
        checkHypre(HYPRE_ParCSRPCGSetAbsoluteTol(ctx.pcg, 0), "PCG abs tol");
        checkHypre(HYPRE_ParCSRPCGSetMaxIter(ctx.pcg, maxIter), "PCG maxIter");
        checkHypre(HYPRE_ParCSRPCGSetPrintLevel(ctx.pcg, c.printLevel), "PCG print");
        checkHypre(HYPRE_ParCSRPCGSetLogging(ctx.pcg, c.logging), "PCG logging");
        checkHypre(HYPRE_ParCSRPCGSetPrecond(ctx.pcg, HYPRE_BoomerAMGSolve,
            HYPRE_BoomerAMGSetup, ctx.amg), "PCGSetPrecond");
    }
    ctx.nRows = nRows;
    ctx.nInternalFaces = matrix.upper().size();
    ctx.vectorsValid = true;
}

void assembleMatrix(HypreSolverContext& ctx)
{
    checkHypre(HYPRE_IJMatrixInitialize(ctx.ijMatrix), "IJMatrixInitialize");
    checkHypre(HYPRE_IJMatrixSetValues(ctx.ijMatrix, ctx.nRows,
        ctx.rowSizes.data(), ctx.rows.data(), ctx.columns.data(),
        ctx.values.data()), "IJMatrixSetValues");
    checkHypre(HYPRE_IJMatrixAssemble(ctx.ijMatrix), "IJMatrixAssemble");
    checkHypre(HYPRE_IJMatrixGetObject(ctx.ijMatrix,
        reinterpret_cast<void**>(&ctx.parMatrix)), "IJMatrixGetObject");
    ctx.matrixValid = true;
    ++ctx.nMatrixAssemblies;
}

void updateVector
(
    HYPRE_IJVector vector,
    HYPRE_ParVector& parVector,
    const label size,
    const std::vector<HYPRE_BigInt>& indices,
    const std::vector<HYPRE_Complex>& values,
    const char* name
)
{
    checkHypre(HYPRE_IJVectorInitialize(vector), name);
    checkHypre(HYPRE_IJVectorSetValues(vector, size, indices.data(),
        values.data()), name);
    checkHypre(HYPRE_IJVectorAssemble(vector), name);
    checkHypre(HYPRE_IJVectorGetObject(vector,
        reinterpret_cast<void**>(&parVector)), name);
}

std::string contextKey
(
    const Foam::lduMatrix& matrix,
    const word& field,
    const Controls& controls,
    const bool flexGmres
)
{
    if (!matrix.mesh().hasDb())
    {
        FatalErrorInFunction
            << "hyprePCGParallel cannot determine case/region identity for field "
            << field << Foam::exit(FatalError);
    }
    const Foam::objectRegistry& db = matrix.mesh().thisDb();
    const Foam::fileName casePath = db.time().globalPath();
    const word region(db.name());
    if (casePath.empty() || region.empty() || field.empty())
    {
        FatalErrorInFunction
            << "hyprePCGParallel incomplete context identity: case='" << casePath
            << "' region='" << region << "' field='" << field << '\''
            << Foam::exit(FatalError);
    }
    std::ostringstream os;
    os << casePath.c_str() << ':' << region.c_str() << ':' << field.c_str()
       << ":mesh@" << static_cast<const void*>(&matrix.mesh())
       << ":sym=" << matrix.symmetric() << ":comm=" << matrix.mesh().comm()
       << (flexGmres ? ":flexGMRES" : ":pcg")
       << ':' << modeName(controls.mode)
       << ':' << controls.coarsenType << ':' << controls.interpType
       << ':' << controls.relaxType << ':' << controls.maxLevels;
    return os.str();
}

}

Foam::hyprePCGParallel::hyprePCGParallel
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
:
    lduMatrix::solver(fieldName, matrix, interfaceBouCoeffs,
        interfaceIntCoeffs, interfaces, solverControls)
{
    const Controls controls = parseControls(controlDict_);
    if (controls.reportLifecycle)
    {
        Info<< "hyprePCGParallel lifecycle: constructed instance="
            << std::uintptr_t(reinterpret_cast<std::uintptr_t>(this))
            << " matrix=" << std::uintptr_t(reinterpret_cast<std::uintptr_t>(&matrix_))
            << " mesh=" << std::uintptr_t(reinterpret_cast<std::uintptr_t>(&matrix_.mesh()))
            << " field=" << fieldName_ << Foam::endl;
    }
}

Foam::hyprePCGParallel::~hyprePCGParallel()
{
    const Controls controls = parseControls(controlDict_);
    if (controls.reportLifecycle)
    {
        Info<< "hyprePCGParallel lifecycle: destroyed instance="
            << std::uintptr_t(reinterpret_cast<std::uintptr_t>(this))
            << " (persistent context retained)" << Foam::endl;
    }
}

Foam::hypreFlexGMRESParallel::hypreFlexGMRESParallel
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
:
    hyprePCGParallel
    (
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces,
        solverControls
    )
{}


Foam::solverPerformance Foam::hyprePCGParallel::solveKrylov
(
    scalarField& psi,
    const scalarField& source,
    const direction cmpt,
    const bool flexGmres
) const
{
    const char* const solverName =
        flexGmres ? "hypreFlexGMRESParallel" : "hyprePCGParallel";
    clockTime totalTimer;
    Timings timing;
    if (!flexGmres && !matrix_.symmetric())
        FatalErrorInFunction
            << "hyprePCGParallel supports symmetric matrices only; use "
            << "hypreFlexGMRESParallel for asymmetric scalar matrices"
            << exit(FatalError);
    const label nRows = matrix_.diag().size();
    if (nRows <= 0 || psi.size() != nRows || source.size() != nRows)
        FatalErrorInFunction << "hyprePCGParallel invalid matrix/vector sizes"
            << exit(FatalError);
    const Controls controls = parseControls(controlDict_);
    ContextRegistry& contexts = registry();
    ensureMeshLifecycle(matrix_.mesh());
    clockTime phase;
    const MatrixSignature sig =
        signature(matrix_, interfaceBouCoeffs_, interfaces_);
    timing.matrixInspection = phase.timeIncrement();
    const std::string key =
        contextKey(matrix_, fieldName_, controls, flexGmres);
    if (controls.reportLifecycle)
    {
        Info<< "hyprePCGParallel context key: " << key << Foam::endl;
    }
    HypreSolverContext& ctx = contexts.get(key);
    label meshTopologyEvent = 0;
    const Foam::objectRegistry& meshDb = matrix_.mesh().thisDb();
    if (meshDb.foundObject<Foam::IOdictionary>("hypreTopologyState"))
    {
        const Foam::IOdictionary& topologyState =
            meshDb.lookupObject<Foam::IOdictionary>("hypreTopologyState");
        meshTopologyEvent =
            topologyState.lookupOrDefault<label>("meshTopologyEvents", 0);
    }
    const scalar matrixScale = matrix_.diag()[0] < 0 ? -1 : 1;
    forAll(matrix_.diag(), rowi)
    {
        if (matrixScale*matrix_.diag()[rowi] <= 0)
        {
            FatalErrorInFunction
                << "hyprePCGParallel requires a consistently signed definite "
                << "matrix diagonal; row " << rowi << " has "
                << matrix_.diag()[rowi] << exit(FatalError);
        }
    }
    ctx.matrixScale = matrixScale;

    const bool first = !ctx.matrixValid;
    int localStructureChanged =
        (!first && ctx.structureHash != sig.structure) ? 1 : 0;
    int localCoefficientsChanged =
        (
            !first
         && !coefficientsEquivalent
            (
                ctx,
                matrix_,
                interfaceBouCoeffs_,
                interfaces_,
                controls,
                sig.coefficients
            )
        ) ? 1 : 0;
    int localTopologyEventAdvanced =
        (meshTopologyEvent > ctx.lastMeshTopologyEvent) ? 1 : 0;
    int structureChangedInt = 0;
    int coefficientsChangedInt = 0;
    int topologyEventAdvancedInt = 0;
    MPI_Allreduce
    (
        &localStructureChanged,
        &structureChangedInt,
        1,
        MPI_INT,
        MPI_MAX,
        MPI_COMM_WORLD
    );
    MPI_Allreduce
    (
        &localCoefficientsChanged,
        &coefficientsChangedInt,
        1,
        MPI_INT,
        MPI_MAX,
        MPI_COMM_WORLD
    );
    MPI_Allreduce
    (
        &localTopologyEventAdvanced,
        &topologyEventAdvancedInt,
        1,
        MPI_INT,
        MPI_MAX,
        MPI_COMM_WORLD
    );
    const bool structureChanged = structureChangedInt;
    const bool coefficientsChanged = coefficientsChangedInt;
    const bool topologyEventAdvanced = topologyEventAdvancedInt;
    if (controls.strictReuse && controls.mode != ReuseMode::none
     && (structureChanged || coefficientsChanged))
    {
        FatalErrorInFunction
            << "hyprePCGParallel strictReuse: matrix changed while reuse was requested"
            << exit(FatalError);
    }

    label callRebuilds = 0, callAssemblies = 0, callSetups = 0, callReuses = 0;
    label callTopologyRebuilds = 0, callStructuralRebuilds = 0;
    const label previousRows = ctx.nRows;
    if (first || structureChanged || controls.mode == ReuseMode::none)
    {
        if (!first) ctx.destroy();
        createContext
        (
            ctx,
            matrix_,
            interfaces_,
            controls,
            tolerance_,
            maxIter_,
            flexGmres
        );
        ++ctx.nContextRebuilds; ++callRebuilds;
        timing.contextCreation = phase.timeIncrement();
        fillMatrixValues(ctx, matrix_, interfaceBouCoeffs_, interfaces_);
        timing.matrixValueUpdate = phase.timeIncrement();
        assembleMatrix(ctx);
        ++callAssemblies;
        timing.matrixAssembly = phase.timeIncrement();
        ctx.amgSetupValid = false;
        if (structureChanged || topologyEventAdvanced)
        {
            ++ctx.nTopologyRebuilds;
            ++ctx.nStructuralRebuilds;
            ++callTopologyRebuilds;
            ++callStructuralRebuilds;
            Info<< solverName << ": matrix structure changed; rebuilding context"
                << " reason=topology oldRows=" << previousRows
                << " newRows=" << nRows
                << " meshTopologyEvent=" << meshTopologyEvent
                << Foam::endl;
        }
    }
    else if (controls.mode == ReuseMode::matrix || coefficientsChanged)
    {
        fillMatrixValues(ctx, matrix_, interfaceBouCoeffs_, interfaces_);
        timing.matrixValueUpdate = phase.timeIncrement();
        assembleMatrix(ctx);
        ++callAssemblies;
        timing.matrixAssembly = phase.timeIncrement();
        if (coefficientsChanged)
        {
            ++ctx.matrixRevision;
            if (controls.mode == ReuseMode::preconditioner && flexGmres)
            {
                Info<< solverName
                    << ": matrix changed; reusing BoomerAMG hierarchy"
                    << Foam::endl;
            }
            else
            {
                ctx.amgSetupValid = false;
                Info<< solverName
                    << ": matrix changed; rebuilding BoomerAMG hierarchy"
                    << Foam::endl;
            }
        }
    }

    ctx.structureHash = sig.structure;
    ctx.coefficientHash = sig.coefficients;
    ctx.lastMeshTopologyEvent = meshTopologyEvent;
    snapshotCoefficients
    (
        ctx,
        matrix_,
        interfaceBouCoeffs_,
        interfaces_
    );

    for (label i = 0; i < nRows; ++i)
        ctx.rhsValues[i] = HYPRE_Complex(ctx.matrixScale*source[i]);
    updateVector(ctx.ijRhs, ctx.parRhs, nRows, ctx.indices, ctx.rhsValues,
        "IJ rhs update");
    timing.rhsUpdate = phase.timeIncrement();
    for (label i = 0; i < nRows; ++i)
        ctx.solutionValues[i] = HYPRE_Complex(psi[i]);
    updateVector(ctx.ijSolution, ctx.parSolution, nRows, ctx.indices,
        ctx.solutionValues, "IJ solution update");
    timing.solutionUpdate = phase.timeIncrement();

    solveScalarField Apsi(nRows), psiSolve(psi), sourceSolve(source), work(nRows);
    matrix_.Amul(Apsi, psiSolve, interfaceBouCoeffs_, interfaces_, cmpt);
    const solveScalar normFactor =
        this->normFactor(psiSolve, sourceSolve, Apsi, work);
    const solveScalarField initialResidualField(sourceSolve - Apsi);
    const solveScalar initialResidual =
        gSumMag(initialResidualField, matrix_.mesh().comm())/normFactor;

    if (!ctx.amgSetupValid)
    {
        if (flexGmres)
        {
            checkHypre
            (
                HYPRE_ParCSRFlexGMRESSetup
                (
                    ctx.pcg, ctx.parMatrix, ctx.parRhs, ctx.parSolution
                ),
                "ParCSRFlexGMRESSetup"
            );
        }
        else
        {
            checkHypre
            (
                HYPRE_ParCSRPCGSetup
                (
                    ctx.pcg, ctx.parMatrix, ctx.parRhs, ctx.parSolution
                ),
                "ParCSRPCGSetup"
            );
        }
        timing.amgSetup = phase.timeIncrement();
        ctx.amgSetupValid = true;
        ++ctx.nAMGSetups; ++callSetups;
        ctx.cumulativeSetupTime += timing.amgSetup;
    }
    else
    {
        ++ctx.nAMGReuses; ++callReuses;
    }
    if (flexGmres)
    {
        checkHypre
        (
            HYPRE_ParCSRFlexGMRESSolve
            (
                ctx.pcg, ctx.parMatrix, ctx.parRhs, ctx.parSolution
            ),
            "ParCSRFlexGMRESSolve"
        );
    }
    else
    {
        checkHypre
        (
            HYPRE_ParCSRPCGSolve
            (
                ctx.pcg, ctx.parMatrix, ctx.parRhs, ctx.parSolution
            ),
            "ParCSRPCGSolve"
        );
    }
    timing.pcgSolve = phase.timeIncrement();
    ctx.cumulativeSolveTime += timing.pcgSolve;

    HYPRE_Int iterations = 0;
    HYPRE_Real finalHypreResidual = 0;
    if (flexGmres)
    {
        checkHypre
        (
            HYPRE_ParCSRFlexGMRESGetNumIterations(ctx.pcg, &iterations),
            "FlexGMRESGetNumIterations"
        );
        checkHypre
        (
            HYPRE_ParCSRFlexGMRESGetFinalRelativeResidualNorm
            (
                ctx.pcg, &finalHypreResidual
            ),
            "FlexGMRESGetFinalResidual"
        );
    }
    else
    {
        checkHypre
        (
            HYPRE_ParCSRPCGGetNumIterations(ctx.pcg, &iterations),
            "PCGGetNumIterations"
        );
        checkHypre
        (
            HYPRE_ParCSRPCGGetFinalRelativeResidualNorm
            (
                ctx.pcg, &finalHypreResidual
            ),
            "PCGGetFinalResidual"
        );
    }
    checkHypre(HYPRE_IJVectorGetValues(ctx.ijSolution, nRows,
        ctx.indices.data(), ctx.solutionValues.data()), "solution copy");
    for (label i = 0; i < nRows; ++i) psi[i] = scalar(ctx.solutionValues[i]);
    timing.solutionCopy = phase.timeIncrement();

    matrix_.Amul(Apsi, solveScalarField(psi), interfaceBouCoeffs_,
        interfaces_, cmpt);
    const solveScalarField finalResidualField(sourceSolve - Apsi);
    const solveScalar finalFoamResidual =
        gSumMag(finalResidualField, matrix_.mesh().comm())/normFactor;
    ++ctx.nSolveCalls;

    if (controls.mode == ReuseMode::none)
    {
        clockTime cleanup;
        ctx.destroy();
        timing.resourceCleanup = cleanup.timeIncrement();
    }
    timing.totalHypreSolve = totalTimer.elapsedTime();
    contexts.account
    (
        callRebuilds, callAssemblies, callSetups, callReuses,
        callTopologyRebuilds, callStructuralRebuilds, label(iterations),
        timing.amgSetup, timing.pcgSolve, timing, controls.reportSummary
    )
    ;

    solverPerformance perf(solverName, fieldName_);
    perf.initialResidual() = initialResidual;
    perf.finalResidual() = scalar(finalHypreResidual);
    perf.nIterations() = label(iterations);
    perf.checkConvergence(tolerance_, relTol_, log_);

    if (controls.reportTimings && controls.reportPerSolve)
    {
        Info<< solverName << " timings [s]: matrixInspection=" << timing.matrixInspection
            << " contextCreation=" << timing.contextCreation
            << " matrixValueUpdate=" << timing.matrixValueUpdate
            << " matrixAssembly=" << timing.matrixAssembly
            << " rhsUpdate=" << timing.rhsUpdate
            << " solutionUpdate=" << timing.solutionUpdate
            << " amgSetup=" << timing.amgSetup
            << " pcgSolve=" << timing.pcgSolve
            << " solutionCopy=" << timing.solutionCopy
            << " resourceCleanup=" << timing.resourceCleanup
            << " totalHypreSolve=" << timing.totalHypreSolve << nl;
    }
    if (controls.reportPerSolve)
    {
        Info<< solverName << " reuse: mode=" << modeName(controls.mode)
            << " solve=" << ctx.nSolveCalls
            << " matrixRevision=" << ctx.matrixRevision
            << " assemblies=" << ctx.nMatrixAssemblies
            << " setups=" << ctx.nAMGSetups << " reuses=" << ctx.nAMGReuses
            << " topologyRebuilds=" << ctx.nTopologyRebuilds
            << " structuralRebuilds=" << ctx.nStructuralRebuilds << nl
            << solverName << " residuals: OpenFOAM(initial)=" << initialResidual
            << " HYPRE(final relative L2)=" << finalHypreResidual
            << " OpenFOAM(final normalized L1)=" << finalFoamResidual
            << Foam::endl;
    }
    return perf;
}


Foam::solverPerformance Foam::hyprePCGParallel::solve
(
    scalarField& psi,
    const scalarField& source,
    const direction cmpt
) const
{
    return solveKrylov(psi, source, cmpt, false);
}


Foam::solverPerformance Foam::hypreFlexGMRESParallel::solve
(
    scalarField& psi,
    const scalarField& source,
    const direction cmpt
) const
{
    return solveKrylov(psi, source, cmpt, true);
}
