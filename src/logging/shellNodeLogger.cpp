#include "shellNodeLogger.h"

#include <utility>
#include "world.h"


shellNodeLogger::shellNodeLogger(std::string logfile_base, std::ofstream &df, int per) :
                             worldLogger("node", std::move(logfile_base), df, per)
{
}


shellNodeLogger::shellNodeLogger(std::string logfile_base, std::string logfile_suffix, std::ofstream &df, int per) :
        worldLogger("node", std::move(logfile_suffix), std::move(logfile_base), df, per)
{
}


shellNodeLogger::~shellNodeLogger() = default;

string shellNodeLogger::getLogHeader() {
    return "";
}

string shellNodeLogger::getLogData() {
    ostringstream log_data;
    log_data << world_ptr->getCurrentTime();
    for (const auto& shell_limb : world_ptr->soft_robots->shell_limbs) {
        for (int i = 0; i < shell_limb->nv; i++) {
            Vector3d v = shell_limb->getVertex(i);
            log_data << "," << v(0) << "," << v(1) << "," << v(2);
        }
    }
    return log_data.str();
}