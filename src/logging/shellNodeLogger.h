#ifndef SHELLNODE_LOGGER_H
#define SHELLNODE_LOGGER_H

#include "worldLogger.h"


class shellNodeLogger : public worldLogger {
public:
    shellNodeLogger(string logfile_base, ofstream& df, int per);
    shellNodeLogger(string logfile_base, string logfile_suffix, ofstream& df, int per);
    ~shellNodeLogger();
private:

    /*
     Populating the log header with hinge connectivity data. This makes it easy to visualize the log data using rendering
     tools outside of openGL. For example, running the simulation on a remote machine and copying over the log data
     that can then be rendered using the mesh connectivity information on a separate machine
    */
    string getLogHeader() override;
    string getLogData() override;
};


#endif
