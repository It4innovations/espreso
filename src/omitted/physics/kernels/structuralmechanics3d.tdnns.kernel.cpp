
#include "structuralmechanics3d.tdnns.kernel.h"

#include "esinfo/stepinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "physics/system/builder/builder.h"
#include "basis/containers/point.h"
#include "basis/evaluator/evaluator.h"
#include "config/ecf/physics/heattransfer.h"
#include "math/matrix.dense.h"
#include "math/vector.dense.h"

#include <cmath>

using namespace espreso;

StructuralMechanics3DTDNNSKernel::StructuralMechanics3DTDNNSKernel(StructuralMechanics3DTDNNSKernel *previous, PhysicsConfiguration &physics, StructuralMechanicsGlobalSettings &gsettings, StructuralMechanicsLoadStepConfiguration &configuration)
: StructuralMechanics3DBaseKernel(previous, physics, gsettings, configuration)
{
//    solutions.push_back(VectorDense(iterator.displacement.output.data->data.size(), iterator.displacement.output.data->data.data()));
}


StructuralMechanics3DTDNNSKernel::~StructuralMechanics3DTDNNSKernel()
{

}

void StructuralMechanics3DTDNNSKernel::assembleLinearElasticMaterialMatrix(esint node, double *coordinates, const MaterialBaseConfiguration *mat, double time, double temp, MatrixDense &K) const
{

}

void StructuralMechanics3DTDNNSKernel::assembleHyperElasticMaterialMatrix(const MaterialBaseConfiguration *mat, MatrixDense &F, MatrixDense &C, MatrixDense &S) const
{

}

void StructuralMechanics3DTDNNSKernel::processElement(const Builder &builder, const ElasticityElementIterator &iterator, InstanceFiller &filler) const
{
    int size = iterator.element->nodes;

    // base functions
    // const std::vector<MatrixDense> &N = *(iterator.element->N);

    printf("element-%d-nodes:\n", (int)iterator.element->nodes);
    for (int n = 0; n < size; n++) {
//        printf(" nIDs: %2d [%6.2f, %6.2f, %6.2f]\n", iterator.nIDs[n], iterator.coordinates[3 * n + 0], iterator.coordinates[3 * n + 1], iterator.coordinates[3 * n + 2]);
    }

    if (builder.matrices & (Builder::Request::K | Builder::Request::R)) {
        filler.Ke.resize(3 * size, 3 * size);
        filler.Ke.fill(0);
    } else {
        filler.Ke.resize(0, 0);
    }
    if (builder.matrices & Builder::Request::C) {
        filler.Ce.resize(3 * size, 3 * size);
        filler.Ce.fill(0);
    } else {
        filler.Ce.resize(0, 0);
    }
    if (builder.matrices & (Builder::Request::M | Builder::Request::R)) {
        filler.Me.resize(3 * size, 3 * size);
        filler.Me.fill(0);
    } else {
        filler.Me.resize(0, 0);
    }
    if (builder.matrices & Builder::Request::R) {
        filler.Re.resize(3 * size);
        filler.Re.fill(0);
    } else {
        filler.Re.resize(0);
    }
    if (builder.matrices & Builder::Request::f) {
        filler.Fe.resize(3 * size);
        filler.Fe.fill(0);
    } else {
        filler.Fe.resize(0);
    }
}

void StructuralMechanics3DTDNNSKernel::processFace(const Builder &builder, const ElasticityBoundaryIterator &iterator, InstanceFiller &filler) const
{

}

void StructuralMechanics3DTDNNSKernel::processEdge(const Builder &builder, const ElasticityBoundaryIterator &iterator, InstanceFiller &filler) const
{

}

void StructuralMechanics3DTDNNSKernel::processNode(const Builder &builder, const ElasticityBoundaryIterator &iterator, InstanceFiller &filler) const
{

}

void StructuralMechanics3DTDNNSKernel::elementSolution(ElasticityElementIterator &iterator)
{

}

void StructuralMechanics3DTDNNSKernel::nodeSolution(ElasticityElementIterator &iterator)
{

}

