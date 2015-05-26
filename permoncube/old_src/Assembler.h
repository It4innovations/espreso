#ifndef ASSEMBLER_H_
#define ASSEMBLER_H_

#include "utility.h"
#include "KSparse.h"
#include "BSparse.h"
#include "RSparse.h"


class CAssembler {

public:
	static void stima_subdomain(double *f ,double *f_internal, 
                              double *u, double *du,
                              CSolid45EqNumbGlobal * stif_loc_numbering,
                              CFem * fem, CDomain *domainG,
                              CKSparse * Ksparse);
};

#endif /* ASSEMBLER_H_ */
