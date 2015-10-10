
#include "pyramid13.h"

using namespace esinput;

size_t Pyramid13::subelements = 6;
size_t Pyramid13::subnodes[] = { 3, 3, 3 };

void Pyramid13::addElements(std::vector<mesh::Element*> &elements, const eslocal indices[])
{
	eslocal pyramid[20];
	pyramid[0] = indices[100];
	pyramid[1] = indices[104];
	pyramid[2] = indices[4];
	pyramid[3] = indices[0];
	pyramid[4] = indices[62];
	pyramid[5] = indices[62];
	pyramid[6] = indices[62];
	pyramid[7] = indices[62];

	pyramid[8] = indices[102];
	pyramid[9] = indices[54];
	pyramid[10] = indices[2];
	pyramid[11] = indices[50];
	pyramid[12] = indices[62];
	pyramid[13] = indices[62];
	pyramid[14] = indices[62];
	pyramid[15] = indices[62];
	pyramid[16] = indices[81];
	pyramid[17] = indices[83];
	pyramid[18] = indices[33];
	pyramid[19] = indices[31];
	elements.push_back(new mesh::Pyramid13(pyramid));

	pyramid[0] = indices[104];
	pyramid[1] = indices[124];
	pyramid[2] = indices[24];
	pyramid[3] = indices[4];

	pyramid[8] = indices[114];
	pyramid[9] = indices[74];
	pyramid[10] = indices[14];
	pyramid[11] = indices[54];
	pyramid[16] = indices[83];
	pyramid[17] = indices[93];
	pyramid[18] = indices[43];
	pyramid[19] = indices[33];
	elements.push_back(new mesh::Pyramid13(pyramid));

	pyramid[0] = indices[124];
	pyramid[1] = indices[120];
	pyramid[2] = indices[20];
	pyramid[3] = indices[24];

	pyramid[8] = indices[122];
	pyramid[9] = indices[70];
	pyramid[10] = indices[22];
	pyramid[11] = indices[74];
	pyramid[16] = indices[93];
	pyramid[17] = indices[91];
	pyramid[18] = indices[41];
	pyramid[19] = indices[43];
	elements.push_back(new mesh::Pyramid13(pyramid));

	pyramid[0] = indices[120];
	pyramid[1] = indices[100];
	pyramid[2] = indices[0];
	pyramid[3] = indices[20];

	pyramid[8] = indices[110];
	pyramid[9] = indices[50];
	pyramid[10] = indices[10];
	pyramid[11] = indices[70];
	pyramid[16] = indices[91];
	pyramid[17] = indices[81];
	pyramid[18] = indices[31];
	pyramid[19] = indices[41];
	elements.push_back(new mesh::Pyramid13(pyramid));

	pyramid[0] = indices[100];
	pyramid[1] = indices[120];
	pyramid[2] = indices[124];
	pyramid[3] = indices[104];

	pyramid[8] = indices[110];
	pyramid[9] = indices[122];
	pyramid[10] = indices[114];
	pyramid[11] = indices[102];
	pyramid[16] = indices[81];
	pyramid[17] = indices[91];
	pyramid[18] = indices[93];
	pyramid[19] = indices[83];
	elements.push_back(new mesh::Pyramid13(pyramid));

	pyramid[0] = indices[4];
	pyramid[1] = indices[24];
	pyramid[2] = indices[20];
	pyramid[3] = indices[0];

	pyramid[8] = indices[14];
	pyramid[9] = indices[22];
	pyramid[10] = indices[10];
	pyramid[11] = indices[2];
	pyramid[16] = indices[33];
	pyramid[17] = indices[43];
	pyramid[18] = indices[41];
	pyramid[19] = indices[31];
	elements.push_back(new mesh::Pyramid13(pyramid));
}



