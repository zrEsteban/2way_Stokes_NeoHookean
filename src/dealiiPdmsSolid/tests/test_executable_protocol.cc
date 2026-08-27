#include "ExecutableProtocol.H"

#include <iostream>
#include <stdexcept>

static void require(const bool condition)
{ if (!condition) throw std::runtime_error("executable protocol unit assertion failed"); }

int main()
{
  using namespace executable_protocol;
  Frame frame; frame.message_type=Type::ForceMessage; frame.sequence=17;
  frame.producer="unit"; frame.time_index=4; frame.outer_corrector=2;
  frame.operator_version=9; frame.payload="canonical\nbytes\n";
  const std::string wire=encode(frame); const Frame decoded=decode(wire);
  require(decoded.message_type==frame.message_type && decoded.sequence==17);
  require(decoded.payload==frame.payload && encode(decoded)==wire);
  for (const auto &mutation: {std::string("bad_magic"),std::string("checksum"),
                              std::string("truncated"),std::string("trailing")})
    {
      std::string damaged=wire;
      if (mutation=="bad_magic") damaged[0]='X';
      if (mutation=="checksum") damaged.back()^=1;
      if (mutation=="truncated") damaged.pop_back();
      if (mutation=="trailing") damaged.push_back('x');
      bool rejected=false; try { (void)decode(damaged); }
      catch (const std::exception &) { rejected=true; }
      require(rejected);
    }
  bool unknown=false;
  try { (void)type("Unknown"); } catch (const std::exception &) { unknown=true; }
  require(unknown);
  require(type("StructuralStateMessage")==Type::StructuralStateMessage);
  require(type("AcceptTimeStep")==Type::AcceptTimeStep);
  Frame too_large; too_large.payload.assign(max_payload+1,'x');
  bool limited=false; try { (void)encode(too_large); }
  catch (const std::exception &) { limited=true; }
  require(limited);
  std::cout << "PASS executable protocol canonical framing/checksums/limits\n";
}
