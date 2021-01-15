
#include "math/math.h"

using namespace espreso;

void MATH::vecSum(esint size, double *result, double alpha, double *a, double beta, double *b)
{
	for (esint i = 0; i < size; ++i) {
		result[i] = alpha * a[i] + beta * b[i];
	}
}

void MATH::vecAddToSparse(esint size, double *result, double alpha, esint *indices, double *other)
{
	for (esint i = 0; i < size; ++i) {
		result[i] += alpha * other[indices[i]];
	}
}

void MATH::CSRTranspose(esint rows, esint cols, esint *aRows, esint *aCols, double *aVals, esint *bRows, esint *bCols, double *bVals)
{
	esint indexing = aRows[0];

	for (esint c = 0; c <= cols; c++) {
		bRows[c] = 0;
	}
	for (esint i = 0; i < aRows[rows] - indexing; i++) {
		++bRows[aCols[i] - indexing + 1];
	}

	for (esint c = 1, sum = 0; c <= cols; c++) {
		esint num = bRows[c];
		bRows[c] = sum + indexing;
		sum += num;
	}
	bRows[0] = indexing;

	for (esint r = 0; r < rows; r++) {
		for (esint c = aRows[r]; c < aRows[r + 1]; c++) {
			bCols[bRows[aCols[c - indexing]] - indexing] = r + indexing;
			bVals[bRows[aCols[c - indexing]] - indexing] = aVals[c - indexing];
			++bRows[aCols[c - indexing]];
		}
	}
}

void MATH::DenseTranspose(esint rows, esint cols, double *vals)
{
	for (esint r = 0; r < rows; ++r) {
		for (esint c = 0; c < cols; ++c) {
			if (r < c) {
				double tmp = vals[r * rows + c];
				vals[r * rows + c] = vals[c * rows + r];
				vals[c * rows + r] = tmp;
			}
		}
	}
}

void MATH::CSRRemoveLower(esint rows, esint cols, esint *aRows, esint *aCols, double *aVals)
{
	esint indexing = aRows[0], total = 0;
	for (esint r = 0, i = 0; r < rows; ++r) {
		esint sum = 0;
		for (esint c = aRows[r]; c < aRows[r + 1]; ++c) {
			if (r <= aCols[c - indexing] - indexing) {
				aVals[i] = aVals[c - indexing];
				aCols[i] = aCols[c - indexing];
				++sum;
				++i;
			}
		}
		aRows[r] = total + indexing;
		total += sum;
	}
	aRows[rows] = total + indexing;
}


