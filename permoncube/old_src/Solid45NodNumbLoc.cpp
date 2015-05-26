/*
 * Solid45NodNumbLoc.cpp
 *
 *  Created on: Oct 17, 2013
 *      Author: mila
 */

#include "Solid45NodNumbLoc.h"

CSolid45NodNumbLoc::CSolid45NodNumbLoc() {
	//int inod_loc[8];
	memset(inod_loc,0, 8 * sizeof(int));// ZeroMemm2(inod_loc, 8 * sizeof(int));

}

CSolid45NodNumbLoc::~CSolid45NodNumbLoc() {
	// TODO Auto-generated destructor stub
}

