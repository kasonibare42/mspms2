/**********************************************************
 * Filename:       inirw.cpp
 * Description:    File reading functions
 * Note:
 * Created by:     Yang Wang   10/28/04
 * Modified by:    
 * Copyright: 
 **********************************************************/

// INCLUDES ///////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <iostream>
#include "inirw.h"

#define _DELETE(X)              if((X)!=NULL) {delete (X); (X)=NULL;}
#define _FREE(X)                if((X)!=NULL) { free(X); (X)=NULL; }
#define ERROR_DATA              -99999999
#define SUCCESS                 0
#define FAIL                    1

using namespace std;

// FUNCTIONS //////////////////////////////////////////////
// initialization
CINISETtag::CINISETtag()
{
    DataLen=0;
    Data=NULL;
    IndexNum=0;
    IndexList=NULL;
}

// initialization
CINISETtag::CINISETtag(char *filename)
{
    DataLen=0;
    Data=NULL;
    IndexNum=0;
    IndexList=NULL;
    Open(filename);
}

// destruction
CINISETtag::~CINISETtag()
{
    if( DataLen != 0 && Data != NULL )
    {
	_DELETE( Data );
    }

    if( IndexNum != 0 && IndexList != NULL )
    {
	_DELETE( IndexList );
    }
}

// read in file
bool CINISETtag::Open(char *filename)
{
    strcpy(FileName, filename);

    _DELETE( Data );

    DataLen = fnGetFileLength(filename);// get the length of the file

    if( DataLen > 0 )	// if file exists 
    {
	Data=new char[DataLen];

	FILE *fp;
	fp=fopen(filename, "rb");
	fread(Data, DataLen, 1, fp);	// read data
	fclose(fp);

	// initialize index
	InitIndex();
    }
    else
    {
	DataLen=1;
	Data=new char[DataLen];
	memset(Data, 0, 1);
	InitIndex();
    }

    return false;
}

// Save file
bool CINISETtag::Save(char *filename)
{
    if( filename==NULL )
	filename=FileName;

    FILE *fp;
    fp=fopen(filename, "wb");
    if( fp==NULL )
    {
	cout << "ERROR: Cannot save " << filename << endl;
	return false;
    }

    fwrite(Data, DataLen, 1, fp);
    fclose(fp);

    return true;
}

// calculate positions of all the index
void CINISETtag::InitIndex()
{
    int i;

    IndexNum=0;

    for(i=0; i<DataLen; i++)
    {
	// find the line beginning with [, and is seperated by return
	// from the previous line
	if( Data[i]=='[' && (i==0||(i>0&&Data[i-1]=='\n')  ) )
	{
	    IndexNum++;
	}
    }

    // allocate the memory
    _DELETE( IndexList );
    if( IndexNum>0 )
	IndexList=new int[IndexNum];

    int n=0;
    // initialize
    for(i=0; i<DataLen; i++)
    {
	if( Data[i]=='[' && ((i>0&&Data[i-1]=='\n' )|| i==0) )
	{
	    // position saved is the next char to the [
	    IndexList[n]=i+1;
	    n++;
	}
    }
}

// return the position of a certian index
int CINISETtag::FindIndex(char *string)
{
    for(int i=0; i<IndexNum; i++)
    {
	char *str=ReadText( IndexList[i] );	// go through the index list
	if( strcmp(string, str) == 0 )		// compare, if they are the same, return the position
	{
	    _FREE( str );
	    return IndexList[i];		// return the position, next char to the [
	}
	_FREE( str );
    }
    return -1;
}

// return the position of a certain data
int CINISETtag::FindData(int index, char *string)
{
    int p=index;	// pointer
    char *name;

    while(1)
    {					// skip the line of index
	p=GotoNextLine(p);		// the position of next line
	name=ReadDataName(p);		// read a data name, p is the pointer
	if( strcmp(string, name)==0 )	// if they are the same, return the position
	{
	    _FREE( name );
	    return p;
	}

	if (name[0] == '[' || strcmp(name, "") ==0)	// if ecounter next index or an empty line, that means
	{			// the name does not exist under the
	    _FREE( name );	// required index
	    return -1;
	}
	    
	_FREE( name );
	if( p>=DataLen ) return -1;	// no proper name found
    }
    return -1;
}

// to the next line
int CINISETtag::GotoNextLine(int p)
{
    int i;

    for(i=p; i<DataLen; i++)
    {
	if( Data[i]=='\n' )		// return to next line
	    return i+1;

    }
    return i;
}

// read a data name at a certain position
char *CINISETtag::ReadDataName(int &p)	// reference, p is the position of the data
{
    char chr;
    char *Ret;
    int m=0;
    int i;

    Ret=(char *)malloc(256);
    memset(Ret, 0, 256);

    for(i=p; i<DataLen; i++)
    {
	chr=Data[i];

	// end
	if( chr == '\r' || chr == '\n' )
	{
	    p=i+1;			
	    return Ret;
	}

	// end
	if( chr == '=' || chr == ';' )
	{
	    p=i+1;			// content of data, skip delimeter
	    return Ret;
	}

	Ret[m]=chr;
	m++;
    }
    return Ret;
}

// read string at a certain position
char *CINISETtag::ReadText(int p)
{
    char chr;
    char *Ret;
    int n=p, m=0;
    int i;

    int EndLine=GotoNextLine(p);	// next line
    Ret=(char *)malloc(EndLine-p+1);	
    memset(Ret, 0, EndLine-p+1);		

    for(i=0; i<DataLen-p; i++)
    {
	chr=Data[n];

	// end, ; return tab or ]
	if( chr == ';' || chr == '\r' || chr == '\t' || chr == ']' || chr == '\n' || chr == ',')
	{
	    return Ret;
	}

	Ret[m]=chr; // read in to buffer
	m++;
	n++;
    }
    return Ret;
}

// read a string via index and data name
char *CINISETtag::ReadText(char *index, char *name)
{
    int n=FindIndex(index);	// find index position
    if( n == -1 )
    {
	printf("ERROR: Cannot find [%s] in file '%s'!\n", index, FileName);		
	return "";
    }

    int m=FindData(n, name);	// find position of data content via index pos and data name
    if( m==-1 )
    {
	printf("ERROR: Cannot find [%s]-'%s' in file '%s'!\n", index, name, FileName);
	return "";
    }

    return ReadText(m);		// read the content
}

// read string at a certain line, does not care the data name
// useful for reading a bunch of data
char *CINISETtag::ReadText(char *index, int lines)
{
    int n=FindIndex(index);	// find index
    if( n == -1 )
    {
	printf("ERROR: Cannot find [%s] in file '%s'!\n", index, FileName);
	return "";
    }

    // go to certain line
    n=GotoNextLine(n);
    for(int i=0; i<lines; i++)
    {
	if( n<DataLen )
	    n=GotoNextLine(n);
    }

    // read data
    while( n<=DataLen )
    {
	if( Data[n] == '=' )	// find delimeter =
	{
	    n++;		// move to the next char to =
	    return ReadText(n);	// read content
	}
	if( Data[n] == '\r' || Data[n] == '\n' )
	{
	    return "";		// nothing found
	}
	n++;
    }

    return "";
}

// read integer via index and data name
int CINISETtag::ReadInt(char *index, char *name)
{
    int n=FindIndex(index);
    if( n == -1 )
    {
	printf("ERROR: Cannot find [%s] in file '%s'!\n", index, FileName);
	return ERROR_DATA;
    }

    int m=FindData(n, name);
    if( m==-1 )
    {
	printf("ERROR: Cannot find [%s]-'%s' in file '%s'!\n", index, name, FileName);
	return ERROR_DATA;
    }

    char *str=ReadText(m);
    int ret=atoi(str);		// convert to integer
    _FREE(str);
    return ret;
}

// read integer at a certain line
// useful read in bunch of data
int CINISETtag::ReadInt(char *index, int lines)
{
    int n=FindIndex(index);
    if( n == -1 )
    {
	printf("ERROR: Cannot find [%s] in file '%s'!\n", index, FileName);
	return ERROR_DATA;
    }

    // go to next line
    n=GotoNextLine(n);
    for(int i=0; i<lines; i++)
    {
	if( n<DataLen )
	    n=GotoNextLine(n);
    }

    // read data
    while( n<DataLen )
    {
	if( Data[n] == '=' )
	{
	    n++;
	    char *str=ReadText(n);
	    int ret=atoi(str);
	    _FREE(str);
	    return ret;
	}
	if( Data[n] == '\r' || Data[n] == '\n' )
	{
	    return ERROR_DATA;
	}
	n++;
    }

    return ERROR_DATA;
}

// read a double number via index and data name
double CINISETtag::ReadDouble(char *index, char *name)
{
    int n=FindIndex(index);
    if( n == -1 )
    {
	printf("ERROR: Cannot find [%s] in file '%s'!\n", index, FileName);
	return ERROR_DATA;
    }

    int m=FindData(n, name);
    if( m==-1 )
    {
	printf("ERROR: Cannot find [%s]-'%s' in file '%s'!\n", index, name, FileName);
	return ERROR_DATA;
    }

    char *str=ReadText(m);
    double ret=atof(str);	// convert to double
    _FREE(str);
    return ret;
}

// read serveral double number via index and data name
double CINISETtag::ReadDoubleMany(char *index, char *name, int num_data, double data_holder[])
{
    char *str;
    int n=FindIndex(index);
    if( n == -1 )
    {
	printf("ERROR: Cannot find [%s] in file '%s'!\n", index, FileName);
	return ERROR_DATA;
    }

    int m=FindData(n, name);
    if( m==-1 )
    {
	printf("ERROR: Cannot find [%s]-'%s' in file '%s'!\n", index, name, FileName);
	return ERROR_DATA;
    }

    for (int ii=0;ii<num_data;ii++)
    {
       	str=ReadText(m);
       	data_holder[ii]=atof(str);	// convert to double
	m += strlen(str)+1;
       	_FREE(str);
    }

    return SUCCESS;
}

// read a double number at a certain line
// useful read in bunch of data
double CINISETtag::ReadDouble(char *index, int lines)
{
    int n=FindIndex(index);
    if( n == -1 )
    {
	printf("ERROR: Cannot find [%s] in file '%s'!\n", index, FileName);
	return ERROR_DATA;
    }

    // go to next line
    n=GotoNextLine(n);
    for(int i=0; i<lines; i++)
    {
	if( n<DataLen )
	    n=GotoNextLine(n);
    }

    // read data
    while( n<DataLen )
    {
	if( Data[n] == '=' )
	{
	    n++;
	    char *str=ReadText(n);
	    double ret=atof(str);
	    _FREE(str);
	    return ret;
	}
	if( Data[n] == '\r' || Data[n] == '\n' )
	{
	    return ERROR_DATA;
	}
	n++;
    }

    return ERROR_DATA;
}

// read data name at a certain line
char *CINISETtag::ReadData(char *index, int lines)
{
    int n=FindIndex(index);
    if( n == -1 )
    {
	printf("ERROR: Cannot find [%s] in file '%s'!\n", index, FileName);
	return NULL;
    }

    // go to next line
    n=GotoNextLine(n);
    for(int i=0; i<lines; i++)
    {
	if( n<DataLen )
	    n=GotoNextLine(n);
    }

    return ReadDataName(n);	// return the name
}

// search a data name under a index
bool CINISETtag::SearchData(char *index, char *name)
{
    int n=FindIndex(index); // find index position
    if (n == -1)
    {
	printf("ERROR: Cannot find [%s] in file '%s'!\n", index, FileName);		
	return false;
    }

    int m=FindData(n, name);	// find position of data content via index pos and data name
    if (m == -1)
	return false;
    else
	return true;
}

// write a string via index, name, content
bool CINISETtag::WriteText(char *index, char *name, char *string)
{
    int n=FindIndex(index);
    if( n == -1 ) // new index
    {
	AddIndex(index);		// add a index
	n=FindIndex(index);		// find the position
	AddData(n, name, string);	// add a new data
	return true;
    }

    // exsit index
    int m=FindData(n, name);
    if( m==-1 )		// new data
    {
	AddData(n, name, string);	// add a new data
	return true;
    }

    // exist data
    ModityData(n, name, string);	// modify data

    return true;
}

// write a integer via index, name, number
bool CINISETtag::WriteInt(char *index, char *name, int num)
{
    char string[32];
    sprintf(string, "%d", num);

    int n=FindIndex(index);
    if( n == -1 )			// new index
    {
	AddIndex(index);
	n=FindIndex(index);
	AddData(n, name, string);	// add a new data
	return true;
    }

    // exist index
    int m=FindData(n, name);
    if( m==-1 )				// new data
    {
	AddData(n, name, string);	// add a new data
	return true;
    }

    // exist data
    ModityData(n, name, string);	// modify data

    return true;
}

// write a integer via index, name, number
bool CINISETtag::WriteDouble(char *index, char *name, double num)
{
    char string[32];
    sprintf(string, "%f", num);

    int n=FindIndex(index);
    if( n == -1 )			// new index
    {
	AddIndex(index);
	n=FindIndex(index);
	AddData(n, name, string);	// add a new data
	return true;
    }

    // exist index
    int m=FindData(n, name);
    if( m==-1 )				// new data
    {
	AddData(n, name, string);	// add a new data
	return true;
    }

    // exist data
    ModityData(n, name, string);	// modify data

    return true;
}

/////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////

// add a index
bool CINISETtag::AddIndex(char *index)
{
    char str[256];
    memset(str, 0, 256);
    int n=FindIndex(index);	// look up for same index

    if( n == -1 )	// new index
    {
	sprintf(str,"\r\n[%s]\r\n",index);
	Data=(char *)realloc(Data, DataLen+strlen(str));// reallocate memory
	sprintf(&Data[DataLen], "%s", str);
	DataLen+=strlen(str);				// renew length

	InitIndex();					// rebuild the index
	return true;
    }

    return false; // exist index
}

// add a data at current position
bool CINISETtag::AddData(int p, char *name, char *string)
{
    char *str;
    int len=strlen(string);
    str=new char[len+256];
    memset(str, 0, len+256);
    sprintf(str,"%s=%s\r\n",name,string);
    len=strlen(str);

    p=GotoNextLine(p);
    Data=(char *)realloc(Data, DataLen+len);

    char *temp=new char[DataLen-p];
    memcpy(temp, &Data[p], DataLen-p);
    memcpy(&Data[p+len], temp, DataLen-p);
    memcpy(&Data[p], str, len);
    DataLen+=len;

    _DELETE( temp );
    _DELETE( str );
    return true;
}

// modify data at current position
bool CINISETtag::ModityData(int p, char *name, char *string)
{
    int n=FindData(p, name);	// find data

    char *t=ReadText(n);	// read data
    p=n+strlen(t);		// move to next position to the end of the data
    _FREE(t);

    int newlen=strlen(string);	// new data length
    int oldlen=p-n;		// old data length

    Data=(char *)realloc(Data, DataLen+newlen-oldlen);	// reallocate memory

    char *temp=new char[DataLen-p];
    memcpy(temp, &Data[p], DataLen-p);
    memcpy(&Data[n+newlen], temp, DataLen-p);
    memcpy(&Data[n], string, newlen);
    DataLen+=newlen-oldlen;

    _DELETE( temp );
    return true;
}

// return file content
char *CINISETtag::GetData()
{
    return Data;
}

// return number of lines of file
int CINISETtag::GetLines(int cur)
{
    int n=1;
    for(int i=0; i<cur; i++)
    {
	if( Data[i]=='\n' )
	    n++;
    }
    return n;
}

// return the continued lines
int CINISETtag::GetContinueDataNum(char *index)
{
    int num=0;
    int n=FindIndex(index);
    n=GotoNextLine(n);
    while(1)
    {
	if( Data[n] == '\r' || Data[n] == -3 || Data[n] == EOF 
		|| Data[n] == ' ' || Data[n] == '/' || Data[n] == '\t' || Data[n] == '\n' )
	{
	    return num;
	}
	else
	{
	    num++;
	    n=GotoNextLine(n);
	    if( n >= DataLen )	
		return num;
	}
    }
}

// Check if a DataName exist under certain index
int CINISETtag::CheckDataName(char *cIdx, char *cName)
{
    int p = FindIndex(cIdx);
    if (p == -1)
    {
	printf("Cannot find [%s] in file '%s'!\n", cIdx, FileName);		
	return FAIL;
    }
    char *name;

    while(1)
    {					// skip the line of index
	p=GotoNextLine(p);		// the position of next line
	name=ReadDataName(p);		// read a data name, p is the pointer
	if( strcmp(cName, name)==0 )	// if they are the same, return the position
	{
	    _FREE( name );
	    return SUCCESS;
	}

	if (name[0] == '[' || strcmp(name, "") ==0)	// if ecounter next index or an empty line, that means
	{			// the name does not exist under the
	    _FREE( name );	// required index
	    return FAIL;
	}
	    
	_FREE( name );
	if( p>=DataLen ) return FAIL;	// no proper name found
    }
    return FAIL;
}

///////////////////////////////////////////////////////////
int fnGetFileLength(char *name)
{
    int nbytes = 0;
    FILE *fp;

    if ((fp = fopen( name, "rb" )) == NULL) // read only
	return -1;
    else
    {
	while(!feof(fp))
	{
	    getc(fp);
	    nbytes++;
	}
    }

    // close file
    fclose(fp);

    // return length of file	
    return nbytes;
}

