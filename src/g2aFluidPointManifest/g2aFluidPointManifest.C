#include "fvCFD.H"
#include "labelIOList.H"
#include "OSspecific.H"
#include <cstdint>
#include <map>
#include <iomanip>
#include <sstream>

using namespace Foam;

static std::uint64_t fnv1a(const std::string &value)
{
    std::uint64_t hash=1469598103934665603ULL;
    for (const unsigned char byte:value) { hash^=byte; hash*=1099511628211ULL; }
    return hash;
}

int main(int argc, char *argv[])
{
    argList::addOption("patch", "name", "interface patch name");
    argList::addOption("output", "path", "canonical TSV output (rank 0)");
    argList::addOption("geometry-contract-hash", "hash", "versioned interface geometry contract hash");
    argList::addOption("fluid-mesh-hash", "hash", "semantic fluid mesh hash");
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    const word patchName(args.getOrDefault<word>("patch", "interface"));
    const fileName output(args.getOrDefault<fileName>("output", "fluid-point-manifest.tsv"));
    if (!args.found("geometry-contract-hash") || !args.found("fluid-mesh-hash"))
        FatalErrorInFunction << "geometry-contract-hash and fluid-mesh-hash are required" << exit(FatalError);
    const word geometryContractHash(args.get<word>("geometry-contract-hash"));
    const word fluidMeshHash(args.get<word>("fluid-mesh-hash"));
    const label patchId=mesh.boundaryMesh().findPatchID(patchName);
    if (patchId<0) FatalErrorInFunction << "unknown patch " << patchName << exit(FatalError);
    const labelList &meshPoints=mesh.boundaryMesh()[patchId].meshPoints();
    labelList globalIds(mesh.nPoints());
    if (Pstream::parRun())
    {
        labelIOList addressing
        (
            IOobject("pointProcAddressing", mesh.facesInstance(), polyMesh::meshSubDir,
                     mesh, IOobject::MUST_READ, IOobject::NO_WRITE)
        );
        if (addressing.size()!=mesh.nPoints())
            FatalErrorInFunction << "pointProcAddressing size mismatch" << exit(FatalError);
        globalIds=addressing;
    }
    else forAll(globalIds, pointi) globalIds[pointi]=pointi;

    std::ostringstream local;
    local << std::setprecision(17);
    forAll(meshPoints, i)
    {
        const label pointi=meshPoints[i]; const point &p=mesh.points()[pointi];
        local << globalIds[pointi] << '\t' << p.x() << '\t' << p.y() << '\t' << p.z()
              << '\t' << Pstream::myProcNo() << '\n';
    }
    List<string> rankPayload(Pstream::nProcs()); rankPayload[Pstream::myProcNo()]=string(local.str());
    Pstream::gatherList(rankPayload);
    if (Pstream::master())
    {
        std::map<label, Tuple2<point,label>> records;
        forAll(rankPayload, ranki)
        {
            std::istringstream input(rankPayload[ranki]);
            while (true)
            {
                label id,owner; scalar x,y,z;
                if (!(input >> id >> x >> y >> z >> owner)) break;
                const point p(x,y,z); auto found=records.find(id);
                if (found==records.end()) records.emplace(id,Tuple2<point,label>(p,owner));
                else
                {
                    if (mag(found->second.first()-p)>1e-14)
                        FatalErrorInFunction << "duplicate global point " << id
                            << " has inconsistent coordinates" << exit(FatalError);
                    found->second.second()=min(found->second.second(),owner);
                }
            }
        }
        std::ostringstream body; body << std::setprecision(17);
        body << "globalFluidPointId\tx\ty\tz\towner\n";
        for (const auto &record:records)
            body << record.first << '\t' << record.second.first().x() << '\t'
                 << record.second.first().y() << '\t' << record.second.first().z()
                 << '\t' << record.second.second() << '\n';
        std::ostringstream checksum; checksum << std::hex << fnv1a(body.str());
        OFstream stream(output); stream.precision(17);
        stream << "# schemaVersion=1 geometryVersion=1 coordinateFrame=reference coordinateUnits=m"
               << " geometryContractHash=" << geometryContractHash
               << " fluidMeshHash=" << fluidMeshHash << " patch=" << patchName
               << " payloadChecksum=" << checksum.str() << nl << body.str();
        Info<< "G2A_QUERY_MANIFEST patch=" << patchName << " patchId=" << patchId
            << " points=" << records.size() << " output=" << output << nl;
    }
    return 0;
}
