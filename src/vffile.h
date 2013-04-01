#ifndef VFFILE_H
#define VFFILE_H

//Functions for files and strings of files.

#include <iostream>
#include <sstream>
#include <fstream>

#include "cad.h"

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

    VFFile();
    VFFile(VFFile&p){}
    VFFile& operator =(VFFile&){}
public:
    static std::string loadFileToString(std::string filename);
    static double getParameterIni(std::string paramname, std::string inifile);
    static double stringToDouble(std::string str);
    static int tryFile(std::string name_of_file);
    static int getYNFromCin();
    static std::string convertDoubleToString(double number);
};

#endif // VFFILE_H
