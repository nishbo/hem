#include "vfrandom.h"

int VFRandom::buf0 = 0;
double VFRandom::buf1 = 0;
double VFRandom::buf2 = 0;
char VFRandom::buf30 = 'a';
string VFRandom::buf31 = "";
stringstream VFRandom::buf32;
FILE *VFRandom::buf40 = NULL;


int VFRandom::getYNFromCin(){
    // 1 for yes; 0 for no
    buf31 = "";
    buf30  = 'a';

    while(1){
        getline(cin, buf31);
        if (buf31.length() == 1){
            buf30 = buf31[0];
            if(buf30=='y' || buf30=='Y') return 1;
            if(buf30=='n' || buf30=='N') return 0;
        }
    }
}

double VFRandom::checkNumberFromCin(){
    //DOES NOT WORK. DO NOT USE
    //checks if valid number was entered. If true - returns it's value, else
    //returns -1
    buf1 = 0.0;
    buf31 = "";

    getline(cin, buf31);
    stringstream myStream(buf31);
    return buf1;
}

string VFRandom::convertDoubleToString(double number){
    buf32.str("");
    buf32 << number;    //add number to the stream
    return buf32.str(); //return a string with the contents of the stream
}

int VFRandom::tryFile(string name_of_file){
    // 1 if OK, 0 if not
    if(ENABLE_TEST)
        return 1;
    buf40 = fopen(name_of_file.c_str(), "r");
    if(buf40){
        cout<<"Overwrite "<<name_of_file<<"? (y/n)";
        buf0 = getYNFromCin();
        if(!buf0){
            fclose(buf40);
            return 0;
        }
    }
    fclose(buf40);
    return 1;
}
