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
    
    inFile.open("Makefile.win");

    if (!inFile) {
        cout << "Unable to open file";
        exit(1); // terminate with error
    }

	outFile.open("Makefile2",ofstream::out);
    if (!outFile) {
        cout << "Unable to open file";
        exit(1); // terminate with error
    }
    
	while (!inFile.eof()) {
		getline(inFile,strLine,'\n');

		///////////////////////////////////////////////////////
		// delete
		//////////////////////////////////////////////////////

		//# Makefile  -> N/A
		iPos=-1;
		iPos=strLine.find("# Makefile");
		if (iPos!=-1) {
			strLine.erase(iPos);
	    }

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

		//CXXFLAGS -> N/A
		iPos=-1;
		iPos=strLine.find("CXXFLAGS");
		if (iPos!=-1) {
			strLine.erase(iPos);
	    }

		//////////////////////////////////////////////////////
		// replace
		//////////////////////////////////////////////////////

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

		//$(LINKOBJ) ... -> ...
		iPos=-1;
		iPos=strLine.find("$(LINKOBJ)");
		if (iPos!=-1) {
			strLine.erase(iPos);
			strLine.append("-static $(OBJ) -o mspms2.x $(LIBS) $(CFLAGS) -O2 -lgfortran -lgsl -lgslcblas -lm");
		}

		cout<<strLine<<endl;
		outFile<<strLine<<endl;
	}
    
    inFile.close();

	return 0;
}
