#ifndef VFFILE_H
#define VFFILE_H

//Functions for files and strings of files.

#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>

class VFFile{
private:
    static int buf0;
    static double buf1, buf2;
    static char buf3;
    static char* buf30;
    static std::string buf31;
    static std::stringstream buf32;
    static FILE* buf40;
    static float buf01;
    static size_t pos1, pos2;

    VFFile();
public:
    static std::string loadFileToString(std::string filename);
    static double getParameterIni(std::string paramname, std::string inifile);
    static double stringToDouble(std::string str);
    static int tryFile(std::string name_of_file);
    static int tryReadFile(std::string name_of_file);
    static int getYNFromCin();
    static std::string convertDoubleToString(double number);
    static std::string getFilenameFromIni(std::string file_with_files, \
                                          std::string filename);
};

#endif // VFFILE_H
