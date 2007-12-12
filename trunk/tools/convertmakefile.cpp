#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;

int main(int argc, char* argv[]) {

    int iPos;
    ifstream inFile;
	ofstream outFile;
	string strLine;
	int iTasos;
    
    inFile.open("Makefile.win");

    if (!inFile) {
        cout << "Unable to open file";
        exit(1); // terminate with error
    }

	outFile.open("Makefile",ofstream::out);
    if (!outFile) {
        cout << "Unable to open file";
        exit(1); // terminate with error
    }
    
	while (!inFile.eof()) {
		getline(inFile,strLine,'\n');

		///////////////////////////////////////////////////////
		// delete
		//////////////////////////////////////////////////////

		//$(RES) -> N/A
		iPos=-1;
		iPos=strLine.find("$(RES)");
		if (iPos!=-1) {
			strLine.erase(iPos);
	    }

		//WINDRES -> N/A
		iPos=-1;
		iPos=strLine.find("WINDRES");
		if (iPos!=-1) {
			strLine.erase(iPos);
	    }

		//RES -> N/A
		iPos=-1;
		iPos=strLine.find("RES");
		if (iPos!=-1) {
			strLine.erase(iPos);
	    }

		//CXXINCS -> N/A
		iPos=-1;
		iPos=strLine.find("CXXINCS");
		if (iPos!=-1) {
			strLine.erase(iPos);
	    }

		//$(INCS) -> N/A
		iPos=-1;
		iPos=strLine.find("$(INCS)");
		if (iPos!=-1) {
			strLine.erase(iPos,7);
	    }

		//INCS -> N/A
		iPos=-1;
		iPos=strLine.find("INCS");
		if (iPos!=-1) {
			strLine.erase(iPos);
	    }
	    
		//$(CXXFLAGS) -> N/A
		iPos=-1;
		iPos=strLine.find("$(CXXFLAGS)");
		if (iPos!=-1) {
			strLine.erase(iPos);
	    }
	    
		//CXXFLAGS -> N/A
		iPos=-1;
		iPos=strLine.find("CXXFLAGS");
		if (iPos!=-1) {
			strLine.erase(iPos);
	    }

		//////////////////////////////////////////////////////
		// replace
		//////////////////////////////////////////////////////

		//# Makefile  -> convertmakefile comments
		iPos=-1;
		iPos=strLine.find("# Makefile");
		if (iPos!=-1) {
			strLine.erase(iPos);
			strLine.append("# Makefile converted by convertmakefile");
	    }

		//gcc.exe -> gcc
		iPos=-1;
		iPos=strLine.find("gcc.exe");
		if (iPos!=-1) {
			strLine.replace(iPos,7,"gcc");
		}

		//g++.exe -> gcc
		iPos=-1;
		iPos=strLine.find("g++.exe");
		if (iPos!=-1) {
			strLine.replace(iPos,7,"g++");
		}
		
		// .exe -> .x     this has to be done after gcc and g++ changes
		iPos=-1;
		iPos=strLine.find(".exe");
		if (iPos!=-1) {
			strLine.replace(iPos,4,".x");
		}
		
		// LIBS -> ...
		iPos=-1;
		iPos=strLine.find("LIBS =");
		if (iPos!=-1) {
			iTasos=-1;
            iTasos=strLine.find("libtasos.a");
			strLine.erase(iPos);
			strLine.append("LIBS = ");
			if (iTasos!=-1){
               strLine.append("mylibtasos/libtasos.a ");
            }
			strLine.append("-static -lgfortran -lgsl -lgslcblas -lm");
		}
				
		cout<<strLine<<endl;
		outFile<<strLine<<endl;
	}
    
    inFile.close();

	return 0;
}
