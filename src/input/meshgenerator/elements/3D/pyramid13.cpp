
#include "pyramid13.h"

using namespace espreso::input;

size_t Pyramid13::subelements = 6;
size_t Pyramid13::subnodes[] = { 3, 3, 3 };

void Pyramid13::addElements(std::vector<Element*> &elements, const eslocal indices[], const eslocal params[])
{
	eslocal pyramid[13];
	pyramid[0] = indices[100];
	pyramid[1] = indices[104];
	pyramid[2] = indices[4];
	pyramid[3] = indices[0];
	pyramid[4] = indices[62];

	pyramid[5] = indices[102];
	pyramid[6] = indices[54];
	pyramid[7] = indices[2];
	pyramid[8] = indices[50];
	pyramid[9] = indices[81];
	pyramid[10] = indices[83];
	pyramid[11] = indices[33];
	pyramid[12] = indices[31];
	elements.push_back(new espreso::Pyramid13(pyramid, 13, params));

	pyramid[0] = indices[104];
	pyramid[1] = indices[124];
	pyramid[2] = indices[24];
	pyramid[3] = indices[4];

	pyramid[5] = indices[114];
	pyramid[6] = indices[74];
	pyramid[7] = indices[14];
	pyramid[8] = indices[54];
	pyramid[9] = indices[83];
	pyramid[10] = indices[93];
	pyramid[11] = indices[43];
	pyramid[12] = indices[33];
	elements.push_back(new espreso::Pyramid13(pyramid, 13, params));

	pyramid[0] = indices[124];
	pyramid[1] = indices[120];
	pyramid[2] = indices[20];
	pyramid[3] = indices[24];

	pyramid[5] = indices[122];
	pyramid[6] = indices[70];
	pyramid[7] = indices[22];
	pyramid[8] = indices[74];
	pyramid[9] = indices[93];
	pyramid[10] = indices[91];
	pyramid[11] = indices[41];
	pyramid[12] = indices[43];
	elements.push_back(new espreso::Pyramid13(pyramid, 13, params));

	pyramid[0] = indices[120];
	pyramid[1] = indices[100];
	pyramid[2] = indices[0];
	pyramid[3] = indices[20];

	pyramid[5] = indices[110];
	pyramid[6] = indices[50];
	pyramid[7] = indices[10];
	pyramid[8] = indices[70];
	pyramid[9] = indices[91];
	pyramid[10] = indices[81];
	pyramid[11] = indices[31];
	pyramid[12] = indices[41];
	elements.push_back(new espreso::Pyramid13(pyramid, 13, params));

	pyramid[0] = indices[100];
	pyramid[1] = indices[120];
	pyramid[2] = indices[124];
	pyramid[3] = indices[104];

	pyramid[5] = indices[110];
	pyramid[6] = indices[122];
	pyramid[7] = indices[114];
	pyramid[8] = indices[102];
	pyramid[9] = indices[81];
	pyramid[10] = indices[91];
	pyramid[11] = indices[93];
	pyramid[12] = indices[83];
	elements.push_back(new espreso::Pyramid13(pyramid, 13, params));

	pyramid[0] = indices[4];
	pyramid[1] = indices[24];
	pyramid[2] = indices[20];
	pyramid[3] = indices[0];

	pyramid[5] = indices[14];
	pyramid[6] = indices[22];
	pyramid[7] = indices[10];
	pyramid[8] = indices[2];
	pyramid[9] = indices[33];
	pyramid[10] = indices[43];
	pyramid[11] = indices[41];
	pyramid[12] = indices[31];
	elements.push_back(new espreso::Pyramid13(pyramid, 13, params));
}



