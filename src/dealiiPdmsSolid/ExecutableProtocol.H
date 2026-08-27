#pragma once

#include "GeneralizedInterfaceLoad.H"

#include <arpa/inet.h>
#include <poll.h>
#include <sys/socket.h>
#include <sys/stat.h>
#include <sys/un.h>
#include <unistd.h>

#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

namespace executable_protocol
{
constexpr std::uint64_t max_payload=64ULL*1024ULL*1024ULL;
constexpr unsigned int schema_version=1;
const std::string magic="G2B1_EXEC";

enum class Type
{
  Hello,Capabilities,DofManifest,OperatorManifest,ForceMessage,TangentMessage,
  ActivateCorrectorState,ClearProvisionalState,Ready,Ack,Nack,AssemblyResult,
  RequestStructuralState,StructuralStateMessage,AcceptTimeStep,RejectTimeStep,
  RollbackToAcceptedState,Shutdown
};

inline std::string name(const Type type)
{
  static const char *names[]={"Hello","Capabilities","DofManifest","OperatorManifest",
    "ForceMessage","TangentMessage","ActivateCorrectorState","ClearProvisionalState",
    "Ready","Ack","Nack","AssemblyResult","RequestStructuralState",
    "StructuralStateMessage","AcceptTimeStep","RejectTimeStep",
    "RollbackToAcceptedState","Shutdown"};
  return names[static_cast<unsigned int>(type)];
}

inline Type type(const std::string &value)
{
  for (unsigned int i=0;i<=static_cast<unsigned int>(Type::Shutdown);++i)
    if (value==name(static_cast<Type>(i))) return static_cast<Type>(i);
  throw std::runtime_error("unknown messageType: "+value);
}

inline std::uint64_t fnv(const std::string &value)
{
  std::uint64_t hash=1469598103934665603ULL;
  generalized_interface::fnv_bytes(hash,value.data(),value.size());
  return hash;
}

struct Frame
{
  Type message_type=Type::Nack;
  std::uint64_t sequence=0,time_index=0,outer_corrector=0,operator_version=0;
  std::string producer="unknown",transfer_mode="dualConservative",payload;
  std::uint64_t payload_checksum=0,message_checksum=0;
};

inline std::string header_without_message_checksum(const Frame &f)
{
  std::ostringstream out;
  out << magic << ' ' << schema_version << ' ' << name(f.message_type) << ' '
      << f.sequence << ' ' << f.payload.size() << ' ' << f.payload_checksum << ' '
      << f.producer << ' ' << f.transfer_mode << ' ' << f.time_index << ' '
      << f.outer_corrector << ' ' << f.operator_version;
  return out.str();
}

inline std::string encode(Frame frame)
{
  if (frame.payload.size()>max_payload) throw std::runtime_error("payloadLength exceeds limit");
  frame.payload_checksum=fnv(frame.payload);
  const std::string prefix=header_without_message_checksum(frame);
  frame.message_checksum=fnv(prefix+"\n"+frame.payload);
  std::ostringstream out; out << prefix << ' ' << frame.message_checksum << '\n' << frame.payload;
  return out.str();
}

inline Frame parse_header(const std::string &line)
{
  std::istringstream in(line); Frame f; std::string frame_magic,type_name;
  unsigned int version=0; std::uint64_t length=0;
  in >> frame_magic >> version >> type_name >> f.sequence >> length
     >> f.payload_checksum >> f.producer >> f.transfer_mode >> f.time_index
     >> f.outer_corrector >> f.operator_version >> f.message_checksum;
  if (!in) throw std::runtime_error("malformed frame header");
  in >> std::ws;
  if (!in.eof()) throw std::runtime_error("trailing frame header fields");
  if (frame_magic!=magic) throw std::runtime_error("bad magic");
  if (version!=schema_version) throw std::runtime_error("schemaVersion mismatch");
  if (length>max_payload) throw std::runtime_error("payloadLength exceeds limit");
  f.message_type=type(type_name); f.payload.assign(static_cast<std::size_t>(length),'\0');
  return f;
}

inline void validate_complete(Frame &f)
{
  if (f.transfer_mode!="dualConservative") throw std::runtime_error("transferMode mismatch");
  if (fnv(f.payload)!=f.payload_checksum) throw std::runtime_error("payload checksum mismatch");
  if (fnv(header_without_message_checksum(f)+"\n"+f.payload)!=f.message_checksum)
    throw std::runtime_error("message checksum mismatch");
}

inline Frame decode(const std::string &bytes)
{
  const auto newline=bytes.find('\n');
  if (newline==std::string::npos) throw std::runtime_error("truncated frame header");
  Frame f=parse_header(bytes.substr(0,newline));
  const std::size_t length=f.payload.size();
  if (bytes.size()!=newline+1+length) throw std::runtime_error("truncated frame or trailing bytes");
  f.payload=bytes.substr(newline+1); validate_complete(f); return f;
}

class UnixChannel
{
public:
  UnixChannel(const std::string &path,const unsigned int timeout_seconds)
    : path(path),timeout_ms(1000*timeout_seconds)
  {
    bool bound=false;
    try
      {
        listener=::socket(AF_UNIX,SOCK_STREAM,0);
        if (listener<0) fail("socket");
        sockaddr_un address{}; address.sun_family=AF_UNIX;
        if (path.size()>=sizeof(address.sun_path)) throw std::runtime_error("socket path too long");
        std::strncpy(address.sun_path,path.c_str(),sizeof(address.sun_path)-1);
        struct stat existing{};
        if (::lstat(path.c_str(),&existing)==0) throw std::runtime_error("protocol endpoint already exists");
        if (errno!=ENOENT) fail("lstat protocol endpoint");
        if (::bind(listener,reinterpret_cast<sockaddr *>(&address),sizeof(address))<0) fail("bind");
        bound=true;
        if (::listen(listener,1)<0) fail("listen");
        pollfd p{listener,POLLIN,0}; const int ready=::poll(&p,1,timeout_ms);
        if (ready==0) throw std::runtime_error("protocol accept timeout");
        if (ready<0) fail("poll accept");
        connection=::accept(listener,nullptr,nullptr); if (connection<0) fail("accept");
      }
    catch (...)
      { if (connection>=0) ::close(connection); if (listener>=0) ::close(listener);
        if (bound) ::unlink(path.c_str());
        throw; }
  }
  ~UnixChannel()
  { if (connection>=0) ::close(connection); if (listener>=0) ::close(listener); ::unlink(path.c_str()); }
  Frame receive()
  {
    std::string line;
    for (;;)
      { const char c=read_byte(); if (c=='\n') break;
        if (line.size()>4096) throw std::runtime_error("frame header too long");
        line.push_back(c); }
    Frame f=parse_header(line);
    read_exact(f.payload.data(),f.payload.size());
    validate_complete(f); return f;
  }
  void send(const Frame &f)
  { const std::string bytes=encode(f); write_exact(bytes.data(),bytes.size()); }
private:
  [[noreturn]] void fail(const std::string &where)
  { throw std::runtime_error(where+": "+std::strerror(errno)); }
  char read_byte() { char c; read_exact(&c,1); return c; }
  void read_exact(char *data,std::size_t count)
  {
    while (count)
      { pollfd p{connection,POLLIN,0}; const int ready=::poll(&p,1,timeout_ms);
        if (ready==0) throw std::runtime_error("protocol read timeout");
        if (ready<0) fail("poll read");
        const ssize_t n=::read(connection,data,count);
        if (n==0) throw std::runtime_error("truncated protocol message");
        if (n<0) fail("read");
        data+=n; count-=static_cast<std::size_t>(n); }
  }
  void write_exact(const char *data,std::size_t count)
  { while (count) { const ssize_t n=::write(connection,data,count); if (n<0) fail("write");
      data+=n; count-=static_cast<std::size_t>(n); } }
  std::string path; int listener=-1,connection=-1,timeout_ms;
};

inline std::vector<std::string> lines(const std::string &payload)
{
  std::istringstream in(payload); std::vector<std::string> result; std::string line;
  while (std::getline(in,line)) result.push_back(line);
  return result;
}

inline generalized_interface::Stamp stamp(const Frame &f,std::istringstream &in)
{
  generalized_interface::Stamp s; s.time_index=f.time_index;
  s.outer_corrector=f.outer_corrector; s.operator_version=f.operator_version;
  s.transfer_mode=f.transfer_mode; s.valid=true;
  std::string key; in >> key >> s.z_version; if (key!="zVersion") throw std::runtime_error("missing zVersion");
  in >> key >> s.hash_graph; if (key!="hashGraph") throw std::runtime_error("missing hashGraph");
  in >> key >> s.hash_weights; if (key!="hashWeights") throw std::runtime_error("missing hashWeights");
  in >> key >> s.dof_manifest_hash; if (key!="dofManifestHash") throw std::runtime_error("missing dofManifestHash");
  return s;
}

inline generalized_interface::ForceMessage parse_force(const Frame &f)
{
  std::istringstream in(f.payload); generalized_interface::ForceMessage m; m.stamp=stamp(f,in);
  std::string key; std::size_t n=0; in >> key >> m.units; if (key!="units") throw std::runtime_error("missing force units");
  in >> key >> n; if (key!="entries" || n>100000000) throw std::runtime_error("invalid force entry count");
  generalized_interface::global_dof_index previous=0;
  for (std::size_t i=0;i<n;++i) { generalized_interface::DofValue v; in >> v.dof >> v.component >> v.value;
    if (!in || !std::isfinite(v.value)) throw std::runtime_error("invalid force scalar");
    if (i && previous>=v.dof) throw std::runtime_error("force entries not canonical/duplicate");
    previous=v.dof;
    m.values.push_back(v); }
  in >> key; if (key!="end") throw std::runtime_error("truncated force payload"); in >> std::ws;
  if (!in.eof()) throw std::runtime_error("trailing force bytes");
  m.checksum=generalized_interface::checksum(m); return m;
}

inline generalized_interface::TangentMessage parse_tangent(const Frame &f)
{
  std::istringstream in(f.payload); generalized_interface::TangentMessage m; m.stamp=stamp(f,in);
  std::string key; std::size_t n=0; in >> key >> m.units; if (key!="units") throw std::runtime_error("missing tangent units");
  in >> key >> n; if (key!="entries" || n>100000000) throw std::runtime_error("invalid tangent entry count");
  std::pair<generalized_interface::global_dof_index,generalized_interface::global_dof_index> previous{0,0};
  for (std::size_t i=0;i<n;++i) { generalized_interface::Triplet v; in >> v.row >> v.column >> v.value;
    if (!in || !std::isfinite(v.value)) throw std::runtime_error("invalid tangent scalar");
    const auto current=std::make_pair(v.row,v.column);
    if (i && !(previous<current)) throw std::runtime_error("tangent entries not canonical/duplicate");
    previous=current; m.values.push_back(v); }
  in >> key; if (key!="end") throw std::runtime_error("truncated tangent payload"); in >> std::ws;
  if (!in.eof()) throw std::runtime_error("trailing tangent bytes");
  m.checksum=generalized_interface::checksum(m); return m;
}
}
