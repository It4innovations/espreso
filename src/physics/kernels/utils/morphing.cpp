
#include "morphing.h"

#include "mesh/mesh.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/elementsregionstore.h"

#include "config/ecf/meshmorphing.h"
#include "config/reader/tokenizer.h"
#include "esinfo/eslog.hpp"
#include "esinfo/envinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"
#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"
#include "basis/evaluator/evaluator.h"
#include "wrappers/mpi/communication.h"
#include "basis/utilities/utils.h"

namespace espreso {
namespace morphing {


static void createTranslationMatrix2D(std::vector<double> &m, const double &x, const double &y)
{
	m.clear();
	m.resize(9);
	m[0 * 3 + 0] = m[1 * 3 + 1] = m[2 * 3 + 2] = 1;
	m[0 * 3 + 2] = x;
	m[1 * 3 + 2] = y;
}

static void createTranslationMatrix3D(std::vector<double> &m, const double &x, const double &y, const double &z)
{
	m.clear();
	m.resize(16);
	m[0 * 4 + 0] = m[1 * 4 + 1] = m[2 * 4 + 2] = m[3 * 4 + 3] = 1;
	m[0 * 4 + 3] = x;
	m[1 * 4 + 3] = y;
	m[2 * 4 + 3] = z;
}

static void createScalingMatrix2D(std::vector<double> &m, const double &x, const double &y)
{
	m.clear();
	m.resize(9);
	m[0 * 3 + 0] = x;
	m[1 * 3 + 1] = y;
	m[2 * 3 + 2] = 1;
}

static void createScalingMatrix3D(std::vector<double> &m, const double &x, const double &y, const double &z)
{
	m.clear();
	m.resize(16);
	m[0 * 4 + 0] = x;
	m[1 * 4 + 1] = y;
	m[2 * 4 + 2] = z;
	m[3 * 4 + 3] = 1;
}

static void createTranslationMatrixToCenter(const CoordinateSystemConfiguration &csystem, std::vector<double> &m)
{
	switch (*csystem.dimension) {
	case DIMENSION::D3: {
		double x,y,z;
		csystem.center.x.evaluator->evalVector(1, Evaluator::Params(), &x);
		csystem.center.y.evaluator->evalVector(1, Evaluator::Params(), &y);
		csystem.center.z.evaluator->evalVector(1, Evaluator::Params(), &z);
		createTranslationMatrix3D(m, x, y, z);
	} break;
	case DIMENSION::D2: {
		double x,y;
		csystem.center.x.evaluator->evalVector(1, Evaluator::Params(), &x);
		csystem.center.y.evaluator->evalVector(1, Evaluator::Params(), &y);
		createTranslationMatrix2D(m, x, y);
	} break;
	default:
		eslog::globalerror("ESPRESO internal error: unsupported operation.\n");
	}
}


static void createTranslationMatrixToZero(const CoordinateSystemConfiguration &csystem, std::vector<double> &m)
{
	switch (*csystem.dimension) {
	case DIMENSION::D3: {
		double x,y,z;
		csystem.center.x.evaluator->evalVector(1, Evaluator::Params(), &x);
		csystem.center.y.evaluator->evalVector(1, Evaluator::Params(), &y);
		csystem.center.z.evaluator->evalVector(1, Evaluator::Params(), &z);
		createTranslationMatrix3D(m, -x, -y, -z);
	} break;
	case DIMENSION::D2: {
		double x,y;
		csystem.center.x.evaluator->evalVector(1, Evaluator::Params(), &x);
		csystem.center.y.evaluator->evalVector(1, Evaluator::Params(), &y);
		createTranslationMatrix2D(m, -x, -y);
	} break;
	default:
		eslog::globalerror("ESPRESO internal error: unsupported operation.\n");
	}
}

static void multiplyTransformationMatrices(const CoordinateSystemConfiguration &csystem, std::vector<double> &left, std::vector<double> &result)
{
	std::vector<double> right;
	right = result;
	result.clear();
	switch (*csystem.dimension) {
	case DIMENSION::D3: {
		result.resize(16);
		MATH::DenseMatDenseMatRowMajorProduct(1, false, 4, 4, left.data(), false, 4, 4, right.data(), 0, result.data());
	} break;
	case DIMENSION::D2: {
		result.resize(9);
		MATH::DenseMatDenseMatRowMajorProduct(1, false, 3, 4, left.data(), false, 3, 3, right.data(), 0, result.data());
	} break;
	default:
		eslog::globalerror("ESPRESO internal error: unsupported operation.\n");
	}
}

static Point applyTransformation(const CoordinateSystemConfiguration &csystem, std::vector<double> &m, const Point &p)
{
	Point result;
	switch (*csystem.dimension) {
	case DIMENSION::D3: {
		result.x = m[0 * 4 + 0] * p.x + m[0 * 4 + 1] * p.y + m[0 * 4 + 2] * p.z + m[0 * 4 + 3];
		result.y = m[1 * 4 + 0] * p.x + m[1 * 4 + 1] * p.y + m[1 * 4 + 2] * p.z + m[1 * 4 + 3];
		result.z = m[2 * 4 + 0] * p.x + m[2 * 4 + 1] * p.y + m[2 * 4 + 2] * p.z + m[2 * 4 + 3];
		result /=  m[3 * 4 + 0] * p.x + m[3 * 4 + 1] * p.y + m[3 * 4 + 2] * p.z + m[3 * 4 + 3];
	} break;
	case DIMENSION::D2: {
		result.x = m[0 * 3 + 0] * p.x + m[0 * 3 + 1] * p.y + m[0 * 3 + 2];
		result.y = m[1 * 3 + 0] * p.x + m[1 * 3 + 1] * p.y + m[1 * 3 + 2];
		result /=  m[2 * 3 + 0] * p.x + m[2 * 3 + 1] * p.y + m[2 * 3 + 2];
	} break;
	default:
		eslog::globalerror("ESPRESO internal error: unsupported operation.\n");
	}
	return result;
}


static void createRotationMatrix(const CoordinateSystemConfiguration &csystem, std::vector<double> &m)
{
	auto d2r = [] (double degree) -> double { return M_PI * degree / 180; };

	Point rPoint;
	switch (*csystem.dimension) {
	case DIMENSION::D3: {
		csystem.rotation.x.evaluator->evalVector(1, Evaluator::Params(), &(rPoint.x));
		csystem.rotation.y.evaluator->evalVector(1, Evaluator::Params(), &(rPoint.y));
		csystem.rotation.z.evaluator->evalVector(1, Evaluator::Params(), &(rPoint.z));
	} break;
	case DIMENSION::Z:
	case DIMENSION::D2: {
		csystem.rotation.z.evaluator->evalVector(1, Evaluator::Params(), &(rPoint.z));
	} break;
	default:
		eslog::globalerror("ESPRESO internal error: unsupported operation.\n");
	}

	switch (*csystem.dimension) {
		case DIMENSION::D3: {
			m.clear();
			m.resize(16,0);
			Point sin, cos;
			switch (csystem.type) {
			case CoordinateSystemConfiguration::TYPE::CARTESIAN:{
				cos.x = std::cos(d2r(rPoint.x));
				cos.y = std::cos(d2r(rPoint.y));
				cos.z = std::cos(d2r(rPoint.z));

				sin.x = std::sin(d2r(rPoint.x));
				sin.y = std::sin(d2r(rPoint.y));
				sin.z = std::sin(d2r(rPoint.z));
			} break;
			default:
				eslog::globalerror("ESPRESO internal error: unsupported operation.\n");
			}
			m[0 * 4 + 0] = cos.y * cos.z;                         m[0 * 4 + 1] = cos.y * sin.z;                         m[0 * 4 + 2] = -sin.y;
			m[1 * 4 + 0] = cos.z * sin.x * sin.y - cos.x * sin.z; m[1 * 4 + 1] = cos.x * cos.z + sin.x * sin.y * sin.z; m[1 * 4 + 2] = cos.y * sin.x;
			m[2 * 4 + 0] = sin.x * sin.z + cos.x * cos.z * sin.y; m[2 * 4 + 1] = cos.x * sin.y * sin.z - cos.z * sin.x; m[2 * 4 + 2] = cos.x * cos.y;
			m[3 * 4 + 3] = 1;
		} break;

		case DIMENSION::Z:
		case DIMENSION::D2: {
			double cos = 0, sin = 0;
			m.clear();
			m.resize(9);

			switch (csystem.type) {
			case CoordinateSystemConfiguration::TYPE::CARTESIAN: {
				cos = std::cos(d2r(rPoint.z));
				sin = std::sin(d2r(rPoint.z));
			} break;
			default:
				eslog::globalerror("ESPRESO internal error: unsupported operation.\n");
			}
			m[0 * 3 + 0] =  cos; m[0 * 3 + 1] = sin;
			m[1 * 3 + 0] = -sin; m[1 * 3 + 1] = cos;
			m[2 * 3 + 2] = 1;
		} break;
		default:
			eslog::globalerror("ESPRESO internal error: unsupported operation.\n");
	}
}

void rbf()
{
	for (auto it = info::ecf->mesh_morphing.rbf.begin(); it != info::ecf->mesh_morphing.rbf.end(); ++it) {
		morphRBF(it->first, it->second, info::mesh->dimension);
	}
}

/*
Applies the defined tranformations in the configuration file and constructs the transformed set of points
*/
void processMorpher(const RBFTargetTransformationConfiguration &target, int dimension, std::vector<Point> &sPoints, size_t startPoint, std::vector<double> &sDisplacement)
{
	size_t pointsToProcess = sPoints.size() - startPoint;

	switch (target.transformation) {
		case MORPHING_TRANSFORMATION::FIXED: {
			sDisplacement.insert(sDisplacement.end(), dimension * sPoints.size(), 0);
		} break;

		case MORPHING_TRANSFORMATION::TRANSLATION: {
			Evaluator::Params params;
			params.coords(3, reinterpret_cast<double*>(sPoints.data() + startPoint));
			sDisplacement.resize(sDisplacement.size() + dimension * pointsToProcess);
			target.translation.x.evaluator->evalVector(pointsToProcess, dimension, params, sDisplacement.data() + dimension * startPoint + 0);
			target.translation.y.evaluator->evalVector(pointsToProcess, dimension, params, sDisplacement.data() + dimension * startPoint + 1);
			if (dimension == 3) {
				target.translation.z.evaluator->evalVector(pointsToProcess, dimension, params, sDisplacement.data() + dimension * startPoint + 2);
			}
		} break;

		case MORPHING_TRANSFORMATION::OFFSET: {
		} break;

		case MORPHING_TRANSFORMATION::SCALING: {

			std::vector<double> result;
			createTranslationMatrixToCenter(target.coordinate_system, result);
			std::vector<double> scaling;
			std::vector<double> sv(3);
			sv[0] = target.scaling.x.evaluator->eval({});
			sv[1] = target.scaling.y.evaluator->eval({});
			sv[2] = target.scaling.z.evaluator->eval({});
			switch (*target.coordinate_system.dimension) {
			case DIMENSION::D3: createScalingMatrix3D(scaling, sv[0], sv[1], sv[2]); break;
			case DIMENSION::D2: createScalingMatrix2D(scaling, sv[0], sv[1]); break;
			default: break;
			}

			std::vector<double> TtoZero;
			createTranslationMatrixToZero(target.coordinate_system, TtoZero);
			multiplyTransformationMatrices(target.coordinate_system, scaling, result);
			multiplyTransformationMatrices(target.coordinate_system, TtoZero, result);

			for(auto it = sPoints.begin() + startPoint; it != sPoints.end(); it++) {
				const Point& p1 = *it;
				Point p2= applyTransformation(target.coordinate_system, result, p1);
				sDisplacement.push_back(p2.x - p1.x);
				sDisplacement.push_back(p2.y - p1.y);
				if (dimension == 3) {
					sDisplacement.push_back(p2.z - p1.z);
				}
			}
		} break;

		case MORPHING_TRANSFORMATION::ROTATION: {
			std::vector<double> result;
			createTranslationMatrixToCenter(target.coordinate_system, result);
			std::vector<double> rotation;
			createRotationMatrix(target.coordinate_system, rotation);
			std::vector<double> TtoZero;
			createTranslationMatrixToZero(target.coordinate_system, TtoZero);
			multiplyTransformationMatrices(target.coordinate_system, rotation, result);
			multiplyTransformationMatrices(target.coordinate_system, TtoZero, result);

			for(auto it = sPoints.begin() + startPoint; it != sPoints.end(); it++) {
				const Point& p1 = *it;
				Point p2 = applyTransformation(target.coordinate_system, result, p1);
				sDisplacement.push_back(p2.x - p1.x);
				sDisplacement.push_back(p2.y - p1.y);
				if (dimension == 3) {
					sDisplacement.push_back(p2.z - p1.z);
				}
			}
		} break;

		default:
			eslog::internalFailure("implement mesh morphing tranformation.\n");
	}
}

esint prepareMatrixM(std::vector<Point> &rPoints,
		std::vector<double> &rDisplacement,
		int dimension, 
		const RBFTargetConfiguration &configuration,
		std::vector<double> &M_values,
		bool use_x, bool use_y, bool use_z
) {

	size_t rowsFromCoordinates = rPoints.size();
	size_t realsize = rowsFromCoordinates;

	M_values.clear();

	for (size_t r = 0; r < rowsFromCoordinates; r++) {
		for(size_t rr = 0; rr <= r; rr++) {
			M_values.push_back(configuration.function.evaluator->evaluate((rPoints[r] - rPoints[rr]).length()));
		}
	}
	if (use_x) {
		for (size_t r = 0; r < rowsFromCoordinates; r++) {
			M_values.push_back(rPoints[r].x);
		}
		M_values.push_back(0);
		realsize++;
	}

	if (use_y) {
		for (size_t r = 0; r < rowsFromCoordinates; r++) {
			M_values.push_back(rPoints[r].y);
		}
		M_values.push_back(0);
		realsize++;
		if (use_x)	M_values.push_back(0);
	}

	if (dimension == 3 && use_z) {
		for (size_t r = 0; r < rowsFromCoordinates; r++) {
			M_values.push_back(rPoints[r].z);
		}
		if (use_x) M_values.push_back(0);
		if (use_y) M_values.push_back(0);
		M_values.push_back(0);
		realsize++;
	}

	for (size_t r = 0; r < rowsFromCoordinates; r++) {
		M_values.push_back(1);
	}
	if (use_x) M_values.push_back(0);
	if (use_y) M_values.push_back(0);
	if (dimension == 3 && use_z) M_values.push_back(0);
	M_values.push_back(0);
	realsize++;
	
	return realsize;
}

void readExternalFile(
		const RBFTargetConfiguration &configuration, int dimension,
		std::map<std::string, std::vector<Point>> &external_data) {

	Tokenizer tokenizer(configuration.external_ffd.path);
	Tokenizer::Token token;

	auto myNextToken = [] (Tokenizer &tokenizer) -> Tokenizer::Token {
		Tokenizer::Token token;
		do {
			token = tokenizer.next();
		}while(token==Tokenizer::Token::LINE_END);
		return token;
	};
	auto myExpectToken =
			[] (Tokenizer &tokenizer, Tokenizer::Token real, Tokenizer::Token expected, const std::string &expectedToken) -> void {
				if (real!=expected) {
					eslog::internalFailure("mesh morphing.\n");
				}
			};
	auto myConvert =
			[] (Tokenizer &tokenizer) -> double {
				size_t size;
				double result = std::stod(tokenizer.value(), &size);
				if (size != tokenizer.value().size()) {
					eslog::internalFailure("mesh morphing.\n");
				}
				return result;
			};

	token = myNextToken(tokenizer);
	while (token != Tokenizer::Token::END) {
		myExpectToken(tokenizer, token, Tokenizer::Token::STRING,
				" a region name.");
		std::string name = tokenizer.value();
		if (external_data.find(name) != external_data.end()) {
			eslog::internalFailure("mesh morphing.\n");
		}
		token = myNextToken(tokenizer);
		myExpectToken(tokenizer, token, Tokenizer::Token::OBJECT_OPEN,
				" symbol \"{\".");
		token = myNextToken(tokenizer);
		while (token != Tokenizer::Token::OBJECT_CLOSE) {
			Point p;
			myExpectToken(tokenizer, token, Tokenizer::Token::STRING,
					" a number.");
			p.x = myConvert(tokenizer);
			token = myNextToken(tokenizer);
			myExpectToken(tokenizer, token, Tokenizer::Token::DELIMITER,
					" symbol \",\".");

			token = myNextToken(tokenizer);
			myExpectToken(tokenizer, token, Tokenizer::Token::STRING,
					" a number.");
			p.y = myConvert(tokenizer);

			token = myNextToken(tokenizer);

			if (dimension == 3) {
				myExpectToken(tokenizer, token, Tokenizer::Token::DELIMITER,
						" symbol \",\".");

				token = myNextToken(tokenizer);
				myExpectToken(tokenizer, token, Tokenizer::Token::STRING,
						" a number.");
				p.z = myConvert(tokenizer);

				token = myNextToken(tokenizer);
			}
			myExpectToken(tokenizer, token, Tokenizer::Token::EXPRESSION_END,
					" symbol \";\".");
			token = myNextToken(tokenizer);

			external_data[name].push_back(p);
		}
		token = myNextToken(tokenizer);
	}
}

void morphRBF(const std::string &name, const RBFTargetConfiguration &configuration, int dimension)
{
	size_t threads = info::env::OMP_NUM_THREADS;
	MATH::setNumberOfThreads(threads);

	if (info::mesh->nodes->originCoordinates == NULL) {
		info::mesh->nodes->originCoordinates = new serializededata<esint, Point>(*info::mesh->nodes->coordinates);
	}
	// else{
		// delete info::mesh->nodes->originCoordinates;
		// info::mesh->nodes->originCoordinates = new serializededata<esint, Point>(*info::mesh->nodes->coordinates);
	// }

	std::vector<Point> sPoints, rPoints;
	std::vector<double> sDisplacement, rDisplacement;

	size_t nSize = 0;
	// printf("dimension: %d\n", dimension);
	for(auto it = configuration.morphers.begin(); it != configuration.morphers.end(); ++it) {
		nSize += info::mesh->bregion(it->first)->nodeInfo.size;
		
		// printf("	[%2d] [%s] nSize: %10d\n", info::mpi::rank, it->first.c_str(), info::mesh->bregion(it->first)->nodeInfo.size);
	}

	sPoints.reserve(nSize);
	sDisplacement.reserve(dimension * nSize);

	for(auto it = configuration.morphers.begin(); it != configuration.morphers.end(); ++it) {
		const BoundaryRegionStore *region = info::mesh->bregion(it->first);
		const auto &nodes = region->nodes->datatarray();
		// const auto &coordinates = info::mesh->nodes->coordinates->datatarray();
		const auto &coordinates = info::mesh->nodes->originCoordinates->datatarray();

		size_t prevsize = sPoints.size();
		for (esint n = 0; n < region->nodeInfo.size; ++n) {
			sPoints.push_back(coordinates[nodes[n + region->nodeInfo.nhalo]]);
		}
		processMorpher(it->second, dimension, sPoints, prevsize, sDisplacement);
	}
	

	if (info::mpi::rank == 0) {
		if (configuration.external_ffd.path != "") {
			std::map<std::string, std::vector<Point>> external_data;

			readExternalFile(configuration, dimension, external_data);

			for(auto it = configuration.external_ffd.morphers.begin(); it != configuration.external_ffd.morphers.end(); ++it) {
				size_t prevsize = sPoints.size();
				if (external_data.find(it->first)==external_data.end()) {
					eslog::error("MORPHING error\n");
				}

				for (auto n = external_data[it->first].begin(); n != external_data[it->first].end(); ++n) {
					sPoints.push_back(*n);
				}
				processMorpher(it->second, dimension, sPoints, prevsize, sDisplacement);
			}
		}
	}

	// printf("	[%2d] nPoints: %10d, nDisplacements: %10d\n", info::mpi::rank, (int)sPoints.size(), (int)sDisplacement.size());
	
	// printf("	[%2d] Point[0]: [%10.8f, %10.8f, %10.8f]\n", info::mpi::rank, sPoints[0].x, sPoints[0].y, sPoints[0].z);
	

	if (!Communication::gatherUnknownSize(sPoints, rPoints)) {
		eslog::internalFailure("gather morphed points.\n");
	}
	// printf("	[%2d] nPoints: %10d, nDisplacements: %10d\n", info::mpi::rank, (int)rPoints.size(), (int)rDisplacement.size());

	if (!Communication::broadcastUnknownSize(rPoints)) {
		eslog::internalFailure("broadcast points.\n");
	}
	// printf("	[%2d] nPoints: %10d, nDisplacements: %10d\n", info::mpi::rank, (int)rPoints.size(), (int)rDisplacement.size());

	if (!Communication::gatherUnknownSize(sDisplacement, rDisplacement)) {
		eslog::internalFailure("gather morphed displacement.\n");
	}
	// printf("	[%2d] nPoints: %10d, nDisplacements: %10d\n", info::mpi::rank, (int)rPoints.size(), (int)rDisplacement.size());
	
	// if(info::mpi::rank == 0){
		// for(int i = 0; i < rPoints.size(); ++i){
			// printf("	    [%10.8f, %10.8f, %10.8f], displacement  [%10.8f, %10.8f, %10.8f]\n", rPoints[i].x, rPoints[i].y, rPoints[i].z, rDisplacement[dimension * i], rDisplacement[dimension * i+1], rDisplacement[dimension * i+2]);
		// }
	// }
	

	std::vector<double> rhs_values;
	std::vector<double> wq_values;
	
	std::vector<double> rhs_values_aca;
	std::vector<double> wq_values_aca;
	
	double aca_eta = configuration.aca_eta;
	double aca_eps = configuration.aca_precision;
	esint base_tree_size = configuration.aca_cluster_tree_leaf_size;
	
	
	bool use_transform_translate = configuration.use_transform_translate;
	bool use_transform_scale = configuration.use_transform_scale;
	bool use_transform_rotate = configuration.use_transform_rotate;
	int regularization_poly_degree = configuration.polynomial_regularization_degree;
	eslog::info("\n =================================== MORPHING SYSTEM ASSEMBLY ================================ \n");
	double aca_start = eslog::time();
	MorphingMatrix M_ACA(
		dimension,
		use_transform_translate,
		use_transform_scale,
		use_transform_rotate,
		regularization_poly_degree,
		sPoints,
		sDisplacement,
		configuration,
		aca_eps,
		aca_eta,
		base_tree_size,
		rhs_values_aca
	);
	eslog::info("   ACA MATRIX DIMENSION                                                         %8d \n", M_ACA.getNRows());
	eslog::info("   ASSEMBLY TIME                                                                %8.3f [s] \n", eslog::time() - aca_start);
	double aca_compr = M_ACA.getCompressionRatio();
	eslog::info("   COMPRESSION                                                                  %8.3f [%%] \n\n", aca_compr * 100.0f);
	// eslog::info("   # of rows           %10d\n", M_ACA.getNRows());
	// eslog::info("   # of columns        %10d\n", M_ACA.getNCols());
	wq_values_aca.resize(dimension * M_ACA.getNCols());	
	
	// M_ACA.print(rhs_values_aca);

	esint itercount_aca = 0;
	std::vector<esint> iters_dim(dimension);
	std::vector<const char*> dim_name = {"X", "Y", "Z", "T"};
	for(esint d = 0 ; d < dimension; d++) {
		esint mm__;
		double sol_time = eslog::time();
		iters_dim[d] += MATH::SOLVER::GMRESolverInternal_ACA(
			M_ACA,
			&rhs_values_aca[d * M_ACA.getNRows()], 
			&wq_values_aca[d * M_ACA.getNCols()],
			configuration.solver_precision, 
			configuration.solver_max_iter, 
			mm__
		) - 1;
		itercount_aca += iters_dim[d];
		eslog::info("   SOLUTION FOR %s                             %6d [iterations]               %8.3f [s] \n", dim_name[d], iters_dim[d], eslog::time() - sol_time);
	}
	
	// if (info::mpi::rank == 0){
		// printf("ACA solution found in %6d iterations (%4d for x, %4d for y, %4d for z)\n", itercount_aca,iters_dim[0],iters_dim[1],iters_dim[2]);
		// printf("sol=[\n");
		// for(esint r = 0; r < M_ACA.getNCols(); ++r){
			// for(esint c = 0; c < dimension; ++c){
				// printf("%10.6f ", wq_values_aca[M_ACA.getNCols() * c + r]);
			// }
			// printf(";\n");
		// }
		// printf("];\n");
	// }
	std::vector<double> error_morph(3);
	std::vector<double> orthogonality_error(3);
	M_ACA.calculateMorphingError(
		wq_values_aca, 
		rhs_values_aca,
		error_morph[0],
		error_morph[1],
		error_morph[2],
		orthogonality_error[0],
		orthogonality_error[1],
		orthogonality_error[2]
	);
	eslog::info("\n");
	for(esint d = 0; d < dimension; ++d){
		eslog::info("   MORPHING ERROR FOR %s                                                           %10.8f [%] \n", dim_name[d], error_morph[d]);
	}
	// eslog::info("\n");
	// for(esint d = 0; d < dimension; ++d){
	// 	eslog::info("   ORTHOGONALITY ERROR FOR %s                                                      %10.8f  \n", dim_name[d], orthogonality_error[d]);
	// }
	
	/*
	std::vector<double> M_values;
	
	if (info::mpi::rank == 0 ||
			(configuration.solver == MORPHING_RBF_SOLVER::ITERATIVE &&
			 info::mpi::rank < dimension &&
			 info::mpi::size >= dimension)) {

		size_t M_size = rPoints.size() + dimension + 1;

		size_t realSize = prepareMatrixM(rPoints, rDisplacement, dimension, configuration, M_values);
		if (realSize != M_size) {
			eslog::internalFailure("error while building matrix M.\n");
		}
		
		switch (configuration.solver) {
		case MORPHING_RBF_SOLVER::ITERATIVE: {


			if (info::mpi::rank == 0) {

				rhs_values.resize(M_size*dimension);
				wq_values.resize(M_size*dimension);

				for (int d = 0; d < dimension; d++) {
					size_t r;
					for (r = 0; r < rPoints.size(); r++) {
						rhs_values[M_size*d + r] = rDisplacement[r * dimension + d];
					}
					for ( ; r < M_size; r++) {
						rhs_values[M_size*d + r] = 0;
					}
				}
				
				// std::vector<double> diff(rhs_values_aca);
				// espreso::MATH::vecAdd(diff.size(), &diff[0], -1.0f, &rhs_values[0]);
				// double err = espreso::MATH::vecDot(diff.size(), &diff[0]);
				// double ref = espreso::MATH::vecDot(rhs_values.size(), &rhs_values[0]);
				// printf("rhs error: %f [%%]\n", 100.0f * err / ref);
			}

			std::vector<esint> iterations;

			if (info::mpi::rank == 0 && info::mpi::size < dimension) {
				for(esint d = 0 ; d < dimension; d++) {
					
					esint itercount;
					MATH::SOLVER::GMRESUpperSymetricColumnMajorMat(
							M_size, &M_values[0],
							&rhs_values[d * M_size], &wq_values[d * M_size],
							configuration.solver_precision, configuration.solver_max_iter, itercount);
					iterations.push_back(itercount);
					
				}
			}else {

				if (info::mpi::rank == 0) {

					for(esint d = 1 ; d < dimension; d++) {
						MPI_Send(&rhs_values[M_size*d], M_size, MPI_DOUBLE,
									d, 0, info::mpi::comm);
					}

				}else{

					rhs_values.resize(M_size);
					wq_values.resize(M_size);

					MPI_Recv(rhs_values.data(),M_size,
							MPI_DOUBLE, 0, 0, info::mpi::comm, MPI_STATUS_IGNORE);

				}

				esint itercount;
				MATH::SOLVER::GMRESUpperSymetricColumnMajorMat(
						M_size, &M_values[0],
						&rhs_values[0], &wq_values[0],
						configuration.solver_precision, configuration.solver_max_iter, itercount);

				iterations.push_back(itercount);
				
				if (info::mpi::rank == 0) {

					for (int d = 1; d < dimension; d++) {
						MPI_Recv(&wq_values[d*M_size],M_size, MPI_DOUBLE,
								d, 0, info::mpi::comm, MPI_STATUS_IGNORE);
						MPI_Recv(&itercount,1, MPI_INT,
										d, 0, info::mpi::comm, MPI_STATUS_IGNORE);
						iterations.push_back(itercount);
					}

				}else {
					MPI_Send(wq_values.data(), M_size, MPI_DOUBLE,
							0, 0, info::mpi::comm);
					MPI_Send(&itercount, 1, MPI_INT,
								0, 0, info::mpi::comm);
				}

			}

		} break;

		case MORPHING_RBF_SOLVER::DIRECT: {

			wq_values.resize(M_size*dimension);

			for (int d = 0; d < dimension; d++) {
				size_t r;
				for (r = 0; r < rPoints.size(); r++) {
					wq_values[M_size*d + r] = rDisplacement[r * dimension + d];
				}
				for ( ; r < M_size; r++) {
					wq_values[M_size*d + r] = 0;
					wq_values_new[M_size*d + r] = 0;
				}
			}

			int result= MATH::SOLVER::directUpperSymetricIndefiniteColumnMajor(
					M_size,  &M_values[0],
					dimension, &wq_values[0]);

					

					

			if (result!=0) {

				bool use_x = true,use_y = true,use_z = true;

				std::vector<double> xs(rPoints.size()), ys(rPoints.size()), zs(rPoints.size());
				for(size_t i=0; i<rPoints.size(); i++) {
					Point &p = rPoints[i];
					xs[i] = p.x;
					ys[i] = p.y;
					if (dimension == 3) {
						zs[i] = p.z;
					}
				}

				auto checkAllSame = [] (std::vector<double> &points) -> bool {
					if (points.size()==0) return false;
					for(auto it = points.begin()+1; it!=points.end();++it) {
						if (points[0]!=*it) {
							return false;
						}
					}
					return true;
				};

				if (checkAllSame(xs)) {
					use_x=false;
				}
				if (checkAllSame(ys)) {
					use_y=false;
				}
				if (dimension == 3) {
					if (checkAllSame(zs)) {
						use_z=false;
					}
				}

				esint realSize = prepareMatrixM(rPoints, rDisplacement, dimension, configuration, M_values, use_x, use_y, use_z);
				
				wq_values.clear();
				wq_values.resize(realSize*dimension);

				for (int d = 0; d < dimension; d++) {
					size_t r;
					for (r = 0; r < rPoints.size(); r++) {
						wq_values[realSize*d + r] = rDisplacement[r * dimension + d];
					}
					for ( ; r < M_size; r++) {
						wq_values[realSize*d + r] = 0;
						wq_values_new[realSize*d + r] = 0;
					}
				}

				int result= MATH::SOLVER::directUpperSymetricIndefiniteColumnMajor(
						realSize,  &M_values[0],
						dimension, &wq_values[0]);
						
				if (result!=0) {
					eslog::error("ESPRESO error: Dense solver is unable to solve the system. Try iterative solver.\n");
				}

				auto insertRowInWQ =
						[&wq_values,dimension] (esint row) -> void {
							esint size = wq_values.size()/dimension;
							wq_values.insert(wq_values.begin()+row,0);
							wq_values.insert(wq_values.begin()+ size +1 + row,0);
							if (dimension ==3) {
								wq_values.insert(wq_values.begin()+ (size +1)*2 + row,0);
							}
						};

				if (!use_x) insertRowInWQ(rPoints.size());
				if (!use_y) insertRowInWQ(rPoints.size()+1);
				if (dimension == 3 && !use_z) insertRowInWQ(rPoints.size()+2);

			}
		} break;

		default:
			eslog::internalFailure("implement mesh morphing solver.\n");
		}
	}

	if (!Communication::broadcastUnknownSize(wq_values)) {
		eslog::internalFailure("broadcast WQ.\n");
	}

	// error of the new approach
	if(info::mpi::rank == 0){
		double diff_ = 0.0f;
		double n2_ = 0.0f;
		
		double diff_max = 0.0f;
		double n_max = 1.0f;
		
		for(esint i_ = 0; i_ < wq_values_aca.size(); ++i_){
			double d_ = wq_values_aca[i_] - wq_values[i_];
			diff_ += d_ * d_;
			
			if(std::abs(wq_values[i_]) > 0.000000005){
				if( std::abs(d_ / wq_values[i_]) > diff_max){
					diff_max = std::abs(d_ / wq_values[i_]);
				}
			}
			
			n2_ += wq_values[i_] * wq_values[i_];
		}
		
		diff_ = std::sqrt(diff_ / n2_);
		
		printf("[ACA] relative error of the solution: %f [%%]\n", diff_ * 100.0f);
	}
	

	// error of the new approach
	if(info::mpi::rank == 0){
		double diff_ = 0.0f;
		double n2_ = 0.0f;
		
		double diff_max = 0.0f;
		double n_max = 1.0f;
		
		for(esint i_ = 0; i_ < wq_values_aca.size(); ++i_){
			double d_ = wq_values_aca[i_] - wq_values[i_];
			diff_ += d_ * d_;
			
			if(std::abs(wq_values[i_]) > 0.000000005){
				if( std::abs(d_ / wq_values[i_]) > diff_max){
					diff_max = std::abs(d_ / wq_values[i_]);
				}
			}
			
			n2_ += wq_values[i_] * wq_values[i_];
		}
		
		diff_ = std::sqrt(diff_ / n2_);
		
		printf("[ACA] relative error of the solution: %f [%%]\n", diff_ * 100.0f);
	}
	*/
	
	
	ElementsRegionStore *tregion = info::mesh->eregion(configuration.target);
	
	// std::vector<Point> coords_mem;
	// coords_mem.reserve(info::mesh->nodes->size);
	// for (esint n = 0; n < info::mesh->nodes->size; ++n) {
		// coords_mem.emplace_back(info::mesh->nodes->coordinates->datatarray()[n]);
	// }

	// {
		// for (esint n = 0; n < info::mesh->nodes->size; ++n) {
			// Point &morphed = info::mesh->nodes->coordinates->datatarray()[n];
			// Point origin = morphed;
			
			// M_ACA.transformPoint(morphed, origin, wq_values_aca, rPoints, configuration, wq_points, wq_points2);
		
		// }
	// }
	
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (auto n = tregion->nodes->begin(t)->begin(); n != tregion->nodes->end(t)->begin(); ++n) {
			Point &morphed = info::mesh->nodes->coordinates->datatarray()[*n];
			// Point morphed = info::mesh->nodes->originCoordinates->datatarray()[*n];

			Point origin = morphed;
			// Point origin = info::mesh->nodes->originCoordinates->datatarray()[*n];
			// info::mesh->nodes->originCoordinates->datatarray()[*n] = origin;
			
			M_ACA.transformPoint(morphed, origin, wq_values_aca, rPoints, configuration);
			
			info::mesh->nodes->coordinates->datatarray()[*n] = morphed;
			
		}
	}
	

	NodeData *morphing = info::mesh->nodes->appendData(3, NamedData::DataType::VECTOR, "RBF_MORPHING");
	double morphing_time = eslog::time();
	#pragma omp parallel
	{
		const auto &origin = info::mesh->nodes->originCoordinates->datatarray();
		const auto &morphed = info::mesh->nodes->coordinates->datatarray();
		
		#pragma omp for
		for (esint n = 0; n < info::mesh->nodes->size; ++n) {
			Point p = morphed[n] - origin[n];
			
			morphing->data[3 * n + 0] = p.x;
			morphing->data[3 * n + 1] = p.y;
			morphing->data[3 * n + 2] = p.z;
			
		}
	}
	eslog::info(" \n");
	eslog::info("   TIME TO MORPH ALL COORDINATES                                              %10.3f [s] \n", eslog::time() - morphing_time);
	eslog::info(" ============================================================================================= \n");

	eslog::checkpoint("MESH: COORDINATES MORPHED");
	eslog::param("morpher", name.c_str());
	eslog::ln();

	MATH::setNumberOfThreads(1);
}

}
}
