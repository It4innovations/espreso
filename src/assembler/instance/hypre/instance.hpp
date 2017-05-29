
#include "instance.h"

#include "../../../mesh/structures/mesh.h"
#include "../../../mesh/elements/element.h"
#include "../../../mesh/structures/coordinates.h"
#include "../../../mesh/settings/property.h"
#include "../../../basis/matrices/denseMatrix.h"

#define TEST 0

namespace espreso {

template <class TPhysics, class TConfiguration>
void HypreInstance<TPhysics, TConfiguration>::init()
{
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Status stat;
//------------------------------------------------------------------------------
	_physics.prepareMeshStructures();
	if (_output.settings || _output.results) {
		// _store.storeGeometry();
	}

	if (_output.settings) {
		_physics.saveMeshProperties(_store);
	}
//------------------------------------------------------------------------------
if(rank==TEST)	std::cout << "  > HYPRE INIT in "<<_physics.pointDOFs.size()<<"D\n";

	// Set the matrix storage type to HYPRE
	char **paramStrings = new char*[1];
	paramStrings[0] = new char[100];
	strcpy(paramStrings[0], "externalSolver HYPRE");
	feiPtr.parameters(1, paramStrings);
	delete paramStrings[0];
	delete [] paramStrings;

	const std::vector<Element*> &elements = _mesh.elements();
	const std::vector<Element*> &nodes = _mesh.nodes();	
	const std::vector<Property> DOFs = _physics.pointDOFs;

	// Set FEI fields = set degree of freedom
	const int nFields = 1;
	int *fieldSizes = new int[nFields];  fieldSizes[0] = DOFs.size();
	int *fieldIDs = new int[nFields];    fieldIDs[0] = 0;

	// Pass the field information to the FEI
	feiPtr.initFields(nFields, fieldSizes, fieldIDs);

	// Elements are grouped into blocks (in this case one block), and we
	// have to describe the number of elements in the block (nElems) as
	// well as the fields (unknowns) per element.
	const int elemBlkID = 0;
	const int nElems = _mesh.elements().size();
	const int elemNNodes = elements[0]->nodes();//TODO elements in one block must have same number of nodes

	const int elemNFields = 0; 						// number of (non-shared) fields per element
	int *elemFieldIDs = NULL; 					// element-fields IDs

	int *nodeNFields = new int[elemNNodes];		// fields per node
	int **nodeFieldIDs = new int*[elemNNodes];	// node-fields IDs
	for (size_t i = 0; i < elemNNodes; i++)
	{
		nodeNFields[i] = 1;
		nodeFieldIDs[i] = new int[nodeNFields[i]];
		nodeFieldIDs[i][0] = fieldIDs[0];		//reference to field that sets the degree of freedom
	}

	// Pass the block information to the FEI. The interleave parameter
	// controls how different fields are ordered in the element matrices.
	const int interleave = 0;
	feiPtr.initElemBlock(elemBlkID, nElems, elemNNodes, nodeNFields, nodeFieldIDs, elemNFields, elemFieldIDs, interleave);
//------------------------------------------------------------------------------
if(rank==TEST)	std::cout << "  > FILL MESS\n";

	// List the global indexes (IDs) of the nodes in each element
	int **elemConn = new int*[elements.size()];
	for (size_t i = 0; i < elements.size(); i++)
	{
		elemConn[i] = new int[elements[i]->nodes()];
		for (size_t j = 0; j < elements[i]->nodes(); j++)
		{
			elemConn[i][j] = _mesh.coordinates().globalIndex(elements[i]->node(j));
//if(rank==TEST)			std::cout << elements[i]->node(j) << " ("<<_mesh.coordinates().globalIndex(elements[i]->node(j)) << ") ";
		}
//if(rank==TEST)		std::cout << "\n";
	}

	// Pass the element topology information to the FEI
	for (size_t i = 0; i < elements.size(); i++)
		feiPtr.initElem(elemBlkID, i, elemConn[i]);

//------------------------------------------------------------------------------
if(rank==TEST)	std::cout << "  > SET SHARED\n";	//size_t coarseNodes()
	
	std::vector <int> SharedIDs;
	std::vector <int> SharedLengs;
	std::vector <int *> SharedProcs;

	for (size_t i=0; i<nodes.size(); i++)
	{
		std::vector<int> shared = nodes[i]->clusters();							//TODO domains/clusters
		if (shared.size() > 1) //more clusters = shared node
		{
			SharedIDs.push_back(_mesh.coordinates().globalIndex(i));
			SharedLengs.push_back(shared.size());
			int * procs = new int[shared.size()];
			for (int j=0; j<shared.size();j++)
				procs[j] = shared[j];
			SharedProcs.push_back(procs);
		}
	}
	// Pass the shared nodes information to the FEI
	if (SharedIDs.size() > 0)
		feiPtr.initSharedNodes(SharedIDs.size(), &SharedIDs[0], &SharedLengs[0], &SharedProcs[0]);

	feiPtr.initComplete();

/*
if(rank==TEST)
{
	std::cout << "SH ("<<rank<<"): "<<SharedIDs;
	for (int i=0; i<SharedIDs.size(); i++)
	{
		std::cout << SharedIDs[i]<< "] ";
		for (int j=0; j<SharedLengs[i]; j++)
			std::cout <<SharedProcs[i][j]<<" ";
		std::cout <<std::endl;
	}
}
// */
//------------------------------------------------------------------------------
	// Specify the boundary conditions
if(rank==TEST)	std::cout << "  > FILL BC\n";

	
	std::vector<int> BCEqn;
	std::vector<int> BCEid;
	for (size_t i = 0; i < nodes.size(); i++)
	{
		if (nodes[i]->hasProperty(DOFs[0], 0))
		{
			BCEqn.push_back(_mesh.coordinates().globalIndex(i));
			BCEid.push_back(i);
//if (rank==TEST) std::cout << _mesh.coordinates().globalIndex(i) << "\n";
		}
	}

	// The arrays alpha, beta and gamma specify the type of boundary
	// condition (essential, natural, mixed). The most general form
	// for Laplace problems is alpha U + beta dU/dn = gamma. In this
	// example we impose zero Dirichlet boundary conditions.
	double **alpha, **beta, **gamma;
	alpha = new double*[BCEqn.size()];
	beta  = new double*[BCEqn.size()];
	gamma = new double*[BCEqn.size()];
	for (size_t i = 0; i < BCEqn.size(); i++)
	{
		alpha[i] = new double[DOFs.size()];
		beta[i]  = new double[DOFs.size()];
		gamma[i] = new double[DOFs.size()];
		for (int j=0; j<DOFs.size(); j++)
		{
			alpha[i][j] = 1.0;
			beta[i][j]  = 0.0;
			gamma[i][j] = 0.0;
		}

		if (nodes[i]->hasProperty(Property::TEMPERATURE, 0))
		{
			const Point &p = _mesh.coordinates()[BCEid[i]];
			gamma[i][0] = nodes[BCEid[i]]->getProperty(Property::TEMPERATURE, BCEid[i], 0, 0, 0, 0);
		}
	}
	
	// Pass the boundary condition information to the FEI
	feiPtr.loadNodeBCs(BCEqn.size(), &BCEqn[0], fieldIDs[0], alpha, beta, gamma);

//------------------------------------------------------------------------------
if(rank==TEST)	std::cout << "  > FILL STIFFNESS\n";


	DenseMatrix Ke;
	std::vector<double> fe;
	double * rhs = new double[elements.size() * elemNNodes * DOFs.size()];
	double ***elemStiff = new double**[elements.size()];
	const std::vector<eslocal> &partition = _mesh.getPartition();
	size_t elemRHSn= 0;
	size_t elemCnt = 0;
	std::vector<eslocal> dofs;

	for (size_t subdomain = 0; subdomain < _mesh.parts(); subdomain++) 
	{
		for (eslocal e = partition[subdomain]; e < partition[subdomain + 1]; e++) //element matrix for each element
		{
			_physics.assembleStiffnessMatrix(elements[e], Ke, fe, dofs);
			
			elemStiff[elemCnt] = new double * [elements[e]->nodes() * DOFs.size()];

			for (size_t nx = 0; nx < elements[e]->nodes(); nx++) 
			{
				for (size_t dx = 0; dx < DOFs.size(); dx++) 
				{
					elemStiff[elemCnt][nx*DOFs.size()+dx] = new double [elements[e]->nodes() * DOFs.size()];
					size_t row = nodes[elements[e]->node(nx)]->DOFIndex(subdomain, dx);				//TODO row?
					for (size_t ny = 0; ny < elements[e]->nodes(); ny++) 
					{
						for (size_t dy = 0; dy < DOFs.size(); dy++)
						{
							elemStiff[elemCnt][nx*DOFs.size()+dx][ny*DOFs.size()+dy] = Ke(dx * elements[e]->nodes() + nx, dy * elements[e]->nodes() + ny);
							size_t column = nodes[elements[e]->node(ny)]->DOFIndex(subdomain, dy);	//TODO column?
//if (rank==TEST) std::cout << " " << std::setprecision(1) << Ke(dx * elements[e]->nodes() + nx, dy * elements[e]->nodes() + ny);
						}
					}
//					std::cout << " rhs: "<< fe[dx * elements[e]->nodes() + nx] << "\n";
					rhs[elemRHSn++] = fe[dx * elements[e]->nodes() + nx];
//if (rank==TEST) std::cout << fe[dx * elements[e]->nodes() + nx];
//if (rank==TEST) std::cout << std::endl;
				}
			}
			elemCnt++;
		}
	}


	// Assemble the matrix. The elemFormat parameter describes
	// the storage (symmetric/non-symmetric, row/column-wise)
	// of the element stiffness matrices.
	const int elemFormat = 0;
	for (size_t i = 0; i < _mesh.elements().size(); i++)
		feiPtr.sumInElem(elemBlkID, i, elemConn[i], elemStiff[i], &rhs[i*elements[i]->nodes()*DOFs.size()], elemFormat);

if(rank==TEST)	std::cout << "  > LOAD COMPLETE\n"; MPI_Barrier(MPI_COMM_WORLD);
	// Finish the FEI load phase
	feiPtr.loadComplete();

//------------------------------------------------------------------------------
if(rank==TEST)	std::cout << "  > CLEANING\n";

	for (size_t i = 0; i<elements.size(); i++)
	{
		for (size_t j=0; j<elemNNodes*DOFs.size(); j++)
			delete [] elemStiff[i][j];
		delete [] elemStiff[i];
	}
	delete [] elemStiff;
	delete [] rhs;
	
	for (size_t i = 0; i < SharedProcs.size(); i++)
		delete [] SharedProcs [i];

	for (size_t i = 0; i < BCEqn.size(); i++)
	{
		delete [] alpha[i];
		delete [] beta[i];
		delete [] gamma[i];
	}
	delete [] alpha;
	delete [] beta;
	delete [] gamma;
	
	for (size_t i = 0; i < elements.size(); i++)
		delete [] elemConn[i];
	delete [] elemConn;

	delete [] nodeNFields;
	for (size_t i = 0; i < elemNNodes; i++)
		delete [] nodeFieldIDs[i];
	delete []nodeFieldIDs;

	delete [] fieldSizes;
	delete [] fieldIDs;

}

template <class TPhysics, class TConfiguration>
void HypreInstance<TPhysics, TConfiguration>::solve(std::vector<std::vector<double> > &solution)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//------------------------------------------------------------------------------
// SET SOLVER PARAMS
	const int DOFsize = _physics.pointDOFs.size();
	const int nParams = 21;
	char **paramStrings = new char*[nParams];
	for (int i = 0; i < nParams; i++)
		paramStrings[i] = new char[100];

	strcpy(paramStrings[0],  "outputLevel 2");
	//	SOLVER
	switch(_configuration.solver)	//CG, GMRES, FGMRES, BICGS, BICGSTAB, TFQMR, SYMQMR, SUPERLU, SUPERLUX
	{
		case HYPRE_SOLVER::CG:
			strcpy(paramStrings[1], "solver cg");
			break;
		default:
		case HYPRE_SOLVER::GMRES:
			strcpy(paramStrings[1], "solver gmres");
			break;
		case HYPRE_SOLVER::FGMRES:
			strcpy(paramStrings[1], "solver fgmres");
			break;

		case HYPRE_SOLVER::BICGS:
			strcpy(paramStrings[1], "solver bicgs");
			break;
		case HYPRE_SOLVER::BICGSTAB:
			strcpy(paramStrings[1], "solver bicgstab");
			break;
		case HYPRE_SOLVER::TFQMR:
			strcpy(paramStrings[1], "solver tfqmr");
			break;
		case HYPRE_SOLVER::SYMQMR:
			strcpy(paramStrings[1], "solver symqmr");
			break;
		case HYPRE_SOLVER::SUPERLU:	//NOT PARALLEL
			strcpy(paramStrings[1], "solver superlu");
			break;
		case HYPRE_SOLVER::SUPERLUX: //NOT PARALLEL
			strcpy(paramStrings[1], "solver superlux");
			break;


	}
	//	PRECONDITIONER
	switch(_configuration.preconditioner) //DIAGONAL PILUT EUCLID PARASAILS BOOMERAMG POLY MLI
	{
		default:
		case HYPRE_PRECONDITIONER::DIAGONAL:
			strcpy(paramStrings[2], "preconditioner diagonal");
			break;
		case HYPRE_PRECONDITIONER::PILUT:
			strcpy(paramStrings[2], "preconditioner pilut");
			break;
		case HYPRE_PRECONDITIONER::EUCLID:
			strcpy(paramStrings[2], "preconditioner euclid");
			break;
		case HYPRE_PRECONDITIONER::PARASAILS:
			strcpy(paramStrings[2], "preconditioner parasails");
			break;
		case HYPRE_PRECONDITIONER::BOOMERAMG:
			strcpy(paramStrings[2], "preconditioner boomeramg");
			break;
		case HYPRE_PRECONDITIONER::POLY:
			strcpy(paramStrings[2], "preconditioner poly");
			break;
		case HYPRE_PRECONDITIONER::MLI:
			strcpy(paramStrings[2], "preconditioner mli");
			break;
	}
	strcpy(paramStrings[3],  ("maxIterations "+std::to_string(_configuration.iterations)).c_str());
	strcpy(paramStrings[4],  ("tolerance "+std::to_string(_configuration.epsilon)).c_str());
	
	strcpy(paramStrings[5],  "gmresDim 30");
	strcpy(paramStrings[6],  "amgNumSweeps 2");	// 1
	strcpy(paramStrings[7],  "amgCoarsenType falgout");
	strcpy(paramStrings[8],  "amgRelaxType hybridsym"); //hybridsym jacobi
	strcpy(paramStrings[9],  ("amgSystemSize "+std::to_string(DOFsize)).c_str());
	strcpy(paramStrings[20], "amgStrongThreshold 0.5"); //1D 0.25, 3D 0.5
	strcpy(paramStrings[10], "amgMaxLevels 20");
	strcpy(paramStrings[11], "amgMeasureType local"); // local / global
	strcpy(paramStrings[12], "MLI smoother HSGS");
	strcpy(paramStrings[13], "MLI numSweeps 1");
	strcpy(paramStrings[14], "MLI smootherWeight 1.0");
	strcpy(paramStrings[15], ("MLI nodeDOF "+std::to_string(DOFsize)).c_str());
	strcpy(paramStrings[16], "MLI nullSpaceDim 1");
	strcpy(paramStrings[17], "MLI minCoarseSize 50");
	strcpy(paramStrings[18], "MLI outputLevel 0");
	strcpy(paramStrings[19], "parasailsSymmetric outputLevel 0");
	

	feiPtr.parameters(nParams, paramStrings);

	for (int i = 0; i < nParams; i++)
		delete [] paramStrings[i];
	delete [] paramStrings;
//------------------------------------------------------------------------------
// SOLVE AND COPY RESULTS

if (rank==TEST)	std::cout << "  > SOLVER HYPRE\n";
	int status;
	feiPtr.solve(&status);
if (rank==TEST)	std::cout << "  > SOLVED - "<< status<<"\n";

	const std::vector<Element*> &nodes = _mesh.nodes();

	int * nodeIDList = new int[nodes.size()];
	int * solnOffsets = new int[nodes.size()];
	feiPtr.getBlockNodeIDList(0, nodes.size(), nodeIDList);

	//results are mixed - for visualisation we must put them to correct order
	size_t sortIndex = 1;
	while(sortIndex < nodes.size() && nodeIDList[sortIndex-1] < nodeIDList[sortIndex])
		sortIndex++;
	
	std::vector <int> sortIndexList;
	int baseIndex = 0;
	for (size_t i=0; i<nodes.size(); i++)
	{
		if (sortIndex == nodes.size() || nodeIDList[baseIndex] < nodeIDList[sortIndex])
			sortIndexList.push_back(baseIndex++);
		else
			sortIndexList.push_back(sortIndex++);
	}
	sortIndex = baseIndex+1;


	//get results	
	double * saver = new double [nodes.size()*DOFsize];
	feiPtr.getBlockNodeSolution(0, nodes.size(), nodeIDList, solnOffsets, saver);

	//solution resize
	solution.resize(_mesh.parts());
	for (size_t p = 0; p < _mesh.parts(); p++)
		solution[p].resize(_mesh.coordinates().localSize(p)*DOFsize);

	
	//solution sort
/*
	//predpoklad ze vsechny subdomeny maji stejny pocet elementu
	size_t partSize = _mesh.coordinates().localSize(0);
//XXX PRESUSPORADANI PRO LIBOVOLNE PROCESU a JEDNU SUBDOMENU
	for (size_t i=0; i<_mesh.coordinates().localSize(0); i++)
		for (int dof=0; dof<DOFsize; dof++)
		{
if (rank==TEST) std::cout << "["<<i*DOFsize+dof<< "] "<<sortIndexList[i]*DOFsize + dof<< " = " << saver[sortIndexList[i]*DOFsize + dof]<<"\n";
			solution[0][i*DOFsize+dof] = saver[sortIndexList[i]*DOFsize + dof];
		}

//XXX PREUSPORADANI JEN PRO JEDEN PROCES a JEDNU SUBDOMENU
	for (size_t p = 0; p < _mesh.parts(); p++)
	{
		for (size_t i=0; i<_mesh.coordinates().localSize(p)*DOFsize; i++)
			solution[p][i] = saver[p*partSize + i];
	}

//XXX PREUSPORADANI JEN PRO JEDEN PROCES A LIBOVOLNY POCET SUBDOMEN
	size_t idx;
	size_t transferIndex = 0;
	for (size_t i=0; i<nodes.size(); i++)
	{
		std::cout << i << "] \n";
		std::vector<int> shared = nodes[i]->domains();
		for (eslocal j=0; j<shared.size(); j++)
		{
			for (size_t k=0; k<DOFsize; k++)
			{
//				std::cout << partSize*DOFsize*shared[j] + nodes[i]->DOFIndex(j, k) << " ";
				idx = partSize*DOFsize*shared[j] + nodes[i]->DOFIndex(shared[j], k);
				solution[idx/(partSize*DOFsize)][idx%(partSize*DOFsize)] = saver[transferIndex];//saver[partSize*j + k];
				std::cout << "["<<idx/(partSize*DOFsize) << "] ["<<idx%(partSize*DOFsize)<<"] = "<< saver[transferIndex] << "\n";
				
				transferIndex++;
			}
			transferIndex -= DOFsize;
		}
		transferIndex += DOFsize;
	}
std::cout << " CLK " << transferIndex << std::endl;
*/

//XXX PRESUSPORADANI PRO LIBOVOLNE PROCESU a LIBOVOLNE SUBDOMEN
	for (size_t i=0; i<nodes.size(); i++)
	{
		std::vector<int> shared = nodes[i]->domains();
		for (eslocal j=0; j<shared.size(); j++)
			for (size_t dof=0; dof<DOFsize; dof++)
				solution[shared[j]][nodes[i]->DOFIndex(shared[j], dof)] = saver[sortIndexList[i]*DOFsize + dof];
	}

	delete [] saver;	
	delete [] nodeIDList;
	delete [] solnOffsets;

        if (_output.results) {
                _physics.saveMeshResults(_store, solution);
        }

	return;
}

}
