#ifndef MACROS_H_
#define MACROS_H_

#define _safealloc(pt,num,size)         pt=calloc(num,size); assert(pt!=NULL)
#define _safefree(pt)	if ((pt)!=NULL) {free(pt); (pt)=NULL;}

#endif /*MACROS_H_*/
