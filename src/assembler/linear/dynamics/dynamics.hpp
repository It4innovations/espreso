
#include "dynamics.h"

namespace espreso {

template <class TInput>
void TransientElasticity<TInput>::init()
{
	Linear<TInput>::init();

	TimeEvent timePrep("Prepare vectors");
	timePrep.startWithBarrier();

	 _u.resize(this->subdomains());
	 _v.resize(this->subdomains());
	 _w.resize(this->subdomains());

	 _un.resize(this->subdomains());
	 _vn.resize(this->subdomains());
	 _wn.resize(this->subdomains());

	 _b.resize(this->subdomains());
	 _tmp.resize(this->subdomains());

	cilk_for (size_t s = 0; s < this->subdomains(); s++) {
		_u[s].resize(this->_M[s].rows, 0);
		_v[s].resize(this->_M[s].rows, 0);
		_w[s].resize(this->_M[s].rows, 0);

		_un[s].resize(this->_M[s].rows, 0);
		_vn[s].resize(this->_M[s].rows, 0);
		_wn[s].resize(this->_M[s].rows, 0);

		_b[s].resize(this->_M[s].rows, 0);
		_tmp[s].resize(this->_M[s].rows, 0);
	}

	cilk_for (size_t s = 0; s < this->subdomains(); s++) {
		for (size_t i = 1; i < _w[s].size(); i += 3) {
			_w[s][i] = 1.0;
		}
	}

	_constantA = {
		  1.0 / (_beta * _deltaT * _deltaT),
		_gama / (_beta * _deltaT),
		  1.0 / (_beta * _deltaT),

		   1.0 / (_beta * 2) - 1.0,
		 _gama /  _beta      - 1.0,
		(_gama / (_beta * 2) - 1.0) * _deltaT,

		(1.0 - _gama) * _deltaT,
		       _gama  * _deltaT
	};

	timePrep.end();
	this->_timeStatistics.addEvent(timePrep);
}

template <class TInput>
void TransientElasticity<TInput>::pre_solve_update()
{
	ESINFO(PROGRESS1) << "Time: " << _time;

	cilk_for (size_t s = 0; s < this->subdomains(); s++) {
		for(size_t i = 0; i < _u[s].size(); i++) {
			_tmp[s][i] = _constantA[0] * _u[s][i] + _constantA[2] * _v[s][i] + _constantA[3] * _w[s][i];
		}
		this->_M[s].MatVec(_tmp[s], _b[s],'N');
	}
}

template <class TInput>
void TransientElasticity<TInput>::post_solve_update()
{
	cilk_for (size_t s = 0; s < this->subdomains(); s++) {
		for(size_t i = 0; i < _u[s].size(); i++) {
			_wn[s][i] = (_constantA[0] * (_un[s][i] - _u[s][i])) - (_constantA[2] * _v[s][i]) - (_constantA[3] * _w[s][i]);
			_vn[s][i] = _v[s][i] + (_constantA[6] * _w[s][i]) + (_constantA[7] * _wn[s][i]);

			_u[s][i] = _un[s][i];
			_v[s][i] = _vn[s][i];
			_w[s][i] = _wn[s][i];
		}
	}

//#ifdef CATALYST
//	unsigned int timeStep = tt;
//	double time = timeStep * dynamic_timestep;
//	Adaptor::CoProcess(input.mesh,l2g_vec, vec_u_n,  time, timeStep, timeStep == numberOfTimeSteps - 1);
//#endif

	std::stringstream ss;
	ss << "results_" << config::env::MPIrank << "_" << _time;

	output::VTK_Full vtk(this->_input.mesh, ss.str());
	vtk.store(_u, this->DOFs(), config::output::SUBDOMAINS_SHRINK_RATIO, config::output::CLUSTERS_SHRINK_RATIO);

//
//	_instance.mesh().store(mesh::VTK_FULL, ss.str(), vec_u_n, 0.95, 0.9);
//	saveVTK(ss.str().c_str(), vec_u_n, l2g_vec, _instance.localBoundaries(), _instance.globalBoundaries(), 0.95, 0.9);
//
//


	_time++;
}

template <class TInput>
void TransientElasticity<TInput>::solve(std::vector<std::vector<double> > &solution)
{
	TimeEvent timeLSrun("Linear Solver - runtime");
	timeLSrun.start();


	this->_lin_solver.Solve(_b, _un);

	timeLSrun.endWithBarrier();
	this->_timeStatistics.addEvent(timeLSrun);
}


}
