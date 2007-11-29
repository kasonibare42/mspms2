/**********************************************************
 * Filename:       inirw.h
 * Description:    Ini file reading and writing functions
 * Note:
 * Created by:     Yang Wang   10/28/04
 * Modified by:    
 * Copyright: 
 **********************************************************/

#ifndef _INIRW_H_
#define _INIRW_H_

// DEFINES ////////////////////////////////////////////////
// PROTOTYPES /////////////////////////////////////////////
// INI (input) file class
typedef class CINISETtag
{
    private:
	char FileName[80];	// file name
	int DataLen;		    	// file length
	char *Data;			// file content

	int IndexNum;			// number of index (number of [])
	int *IndexList;			// List of index positions

	int Point;			// current pointer
	int Line, Word;			// current line and column

    public:
	CINISETtag();
	CINISETtag(char *);		// Initial and open input file
	~CINISETtag();			// free memory
	char *GetData();		// return content of file
	int GetLines(int);		// return number of lines of the file

	bool Open(char *);		// open the file
	bool Save(char *filename=NULL);	// save the file

    private:
	void InitIndex();		// initial index
	int FindIndex(char *);		// find the positin of index
	int FindData(int, char *);	// return the position of data
	int GotoNextLine(int); 		// to the next line
	char *ReadDataName(int &);	// read the data name in certain position
	char *ReadText(int);		// read string in certain position

	bool AddIndex(char *);		// add a index
	bool AddData(int, char *, char *);	// add data in current position
	bool ModityData(int, char *, char *); // modify data in current position

    public:
	int ReadInt(char *, char *);	// read a integer
	double ReadDouble(char *, char *); // read a double number
	double ReadDoubleMany(char *, char *, int, double *); // read serveral double number
	char *ReadText(char *, char *);	// read a string
	int ReadInt(char *, int );	// read a integer in certain line
	double ReadDouble(char *, int );// read a double number at a certain line
	char *ReadText(char *, int);	// read a string in certain line
	char *ReadData(char *, int);	// read a data name in certain line

	bool SearchData(char *, char *); // search a data name

	bool WriteInt(char *, char *, int);	// write out a integer
	bool WriteDouble(char *, char *, double);  // write out a double number
	bool WriteText(char *, char *, char *);	// write out a string

	int GetContinueDataNum(char *);	// return number of continued lines

	int CheckDataName(char *, char *); // check if a name under certain index 
} INISET, *PINISET;


int fnGetFileLength(char *name);

#endif
