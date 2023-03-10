#include <fstream>

namespace ouin
{

inline void wtofile(const std::string& filename, const ::std::string& in)
{
    std::fstream outfile;
    outfile.open(filename, std::fstream::in | std::fstream::out | std::fstream::app);

    outfile << in << "\n";
    outfile.close();
}

} // namespace ouin