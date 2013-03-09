#ifndef VFRANDOM_H
#define VFRANDOM_H

/// Contains various helpful functions.

#include <iostream>
#include <sstream>
#include <cstdio>
#include "cad.h"

using namespace std;

class VFRandom{
private:
    static int buf0;
    static double buf1, buf2;
    static char buf30;
    static string buf31;
    static stringstream buf32;
    static FILE* buf40;

    VFRandom();
    VFRandom(VFRandom&p){}
    VFRandom& operator =(VFRandom&){}

public:
    static int getYNFromCin();
    static double checkNumberFromCin();
    static string convertDoubleToString(double number);
    static int tryFile(string name_of_file);
};

#endif // VFRANDOM_H
