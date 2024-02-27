#ifndef SHELLNODE_LOGGER_H
#define SHELLNODE_LOGGER_H

#include "worldLogger.h"


class shellNodeLogger : public worldLogger {
public:
    shellNodeLogger(string logfile_base, ofstream& df, int per);
    shellNodeLogger(string logfile_base, string logfile_suffix, ofstream& df, int per);
    ~shellNodeLogger();
private:
    string getLogHeader() override;
    string getLogData() override;
};


#endif
