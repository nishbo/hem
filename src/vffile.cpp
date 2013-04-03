#include "vffile.h"

using namespace std;

int VFFile::buf0 = 0;
float VFFile::buf01 = 0;
char VFFile::buf3 = 'a';
char* VFFile::buf30 = new char[300];
string VFFile::buf31 = "";
stringstream VFFile::buf32;
FILE *VFFile::buf40 = NULL;
size_t VFFile::pos1;
size_t VFFile::pos2;

int VFFile::getYNFromCin(){
    // 1 for yes; 0 for no
    buf31 = "";
    buf3  = 'a';

    while(1){
        getline(cin, buf31);
        if (buf31.length() == 1){
            buf3 = buf31[0];
            if(buf3=='y' || buf3=='Y') return 1;
            if(buf3=='n' || buf3=='N') return 0;
        }
    }
}

string VFFile::convertDoubleToString(double number){
    buf32.str("");
    buf32 << number;    //add number to the stream
    return buf32.str(); //return a string with the contents of the stream
}

string VFFile::loadFileToString(string filename){
    ifstream in(filename.c_str());
    ostringstream out;
    out << in.rdbuf();
    buf31 = out.str();
    in.close();
    return buf31;
}

double VFFile::getParameterIni(string paramname, string inifile){
    pos2 = inifile.find("=", inifile.find(paramname)) + 1; //Find start

    if(inifile.find(";", pos2) < inifile.find("\n", pos2))  //Find end
        pos1 = inifile.find(";", pos2) - 1;
    else
        pos1 = inifile.find("\n", pos2) - 1;

    while(!(isdigit(buf30[0] = inifile[pos2]) || buf30[0]=='.' || buf30[0]=='-') \
          && pos2<=pos1)
        pos2++;   //Find exact start - remove symbols

    inifile.copy(buf30, pos1 - pos2 + 1, pos2);

    buf0 = (int)(pos1 - pos2 + 1);
    buf30[buf0] = '\0';
    while(!isdigit(buf30[buf0-1])){ //Find exact end - remove symbols
        buf0--;
        buf30[buf0] = '\0';
    }

    buf01 = atof(buf30);

    return buf01;
}

int VFFile::tryReadFile(string name_of_file){
    //Return 1 if file exists
    buf40 = fopen(name_of_file.c_str(), "r");
    if(!buf40)
        return 0;
    fclose(buf40);
    return 1;
}

int VFFile::tryFile(string name_of_file){
    // 1 if OK, 0 if not

    if(tryReadFile(name_of_file)){
        cout<<"Overwrite "<<name_of_file<<"? (y/n)";
        buf0 = getYNFromCin();
        if(!buf0){
            fclose(buf40);
            return 0;
        }
    }

    return 1;
}

string VFFile::getFilenameFromIni(string file_with_files, string filename){
    buf31 = loadFileToString(file_with_files);

    pos2 = buf31.find("=", buf31.find(filename)) + 2;   //Find start
    if(buf31.find(";", pos2) < buf31.find("\n", pos2))  //Find end
        pos1 = buf31.find(";", pos2) - 1;
    else
        pos1 = buf31.find("\n", pos2) - 1;

    buf31.copy(buf30, pos1 - pos2 + 1, pos2);
    buf30[(int)(pos1 - pos2 + 1)] = '\0';

    return (string)buf30;
}

string VFFile::loadFileToString(string filename, string file_with_files){
    buf31 = getFilenameFromIni(file_with_files, filename);
    if(!(tryReadFile(buf31))){
        cout<<"Problems reading a file named "<<buf31<<", key "<<filename<<".";
        exit(74);
    }
    return loadFileToString(buf31);
}

string vf_file::loadFileToString(string filename){
    return VFFile::loadFileToString(filename);
}
double vf_file::getParameterIni(string paramname, string inifile){
    return VFFile::getParameterIni(paramname, inifile);
}
int vf_file::tryFile(string name_of_file){
    return VFFile::tryFile(name_of_file);
}
int vf_file::tryReadFile(string name_of_file){
    return VFFile::tryReadFile(name_of_file);
}
int vf_file::getYNFromCin(){
    return VFFile::getYNFromCin();
}
string vf_file::convertDoubleToString(double number){
    return VFFile::convertDoubleToString(number);
}
string vf_file::getFilenameFromIni(string file_with_files, string filename){
    return VFFile::getFilenameFromIni(file_with_files, filename);
}
string vf_file::loadFileToString(string filename, string file_with_files){
    return VFFile::loadFileToString(filename, file_with_files);
}
