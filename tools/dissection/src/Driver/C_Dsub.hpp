/*! \file   C_Dsub.hpp
    \brief  routines for substiution of off-diagonal matrix with strips
    \author Atsushi Suzuki, Laboratoire Jacques-Louis Lions
    \date   Jun. 20th 2014
    \date   Jul. 12th 2015
    \date   Feb. 29th 2016
*/

// This file is part of Dissection
// 
// Dissection is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Dissection is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Dissection.  If not, see <http://www.gnu.org/licenses/>.

#include "Compiler/blas.hpp"
#include "Compiler/OptionLibrary.h"
#include "Driver/C_threads_tasks.hpp"

template<typename T>
void C_Dsub_task_exec(void *arg_); 

template<typename T>
void dsub_sym2sym_diag(C_Dsub_task<T> *arg);

template<typename T>
void dsub_unsym2unsym_diag(C_Dsub_task<T> *arg);

template<typename T>
void dsub_sym2sym(C_Dsub_task<T> *arg);

template<typename T>
void dsub_unsym2unsym(C_Dsub_task<T> *arg);

template<typename T>
void dsub_unsym2diag(C_Dsub_task<T> *arg);

template<typename T>
void dsub_sym2rct(C_Dsub_task<T> *arg);

template<typename T>
void dsub_unsym2rct(C_Dsub_task<T> *arg);

template<typename T>
void dsub_sym2sym_diag_two(C_Dsub_task<T> *arg);

template<typename T>
void dsub_unsym2unsym_diag_two(C_Dsub_task<T> *arg);

template<typename T>
void dsub_sym2sym_two(C_Dsub_task<T> *arg);

template<typename T>
void dsub_unsym2unsym_two(C_Dsub_task<T> *arg);

template<typename T>
void dsub_unsym2diag_two(C_Dsub_task<T> *arg);

template<typename T>
void dsub_sym2rct_two(C_Dsub_task<T> *arg);

template<typename T>
void dsub_unsym2rct_two(C_Dsub_task<T> *arg);

template<typename T>
void C_Dsub_queue(bool isSym, 
		  int father_id,
		  bool skip_flag,
		  vector<C_task *>& queue,
		  list <child_contribution<T> > &child_contrib,
		  vector<C_task *>* tasks_p, // _tasks_DSymmGEMM
		  vector<int>* tasks_p_indcol,
		  vector<C_task *>* tasks_q, // _tasks_DfillSymm
		  vector<C_task *>* tasks_r, // _tasks_SparseLocalSchur
		  vector<C_task *>* tasks_s, // _tasks_DSub[level + 1][(*it)]
		  vector<C_task *>* tasks_d, // _tasks_deallocateLocalSchur
		  vector<int>* tasks_d_indcol,
		  int level,
		  const bool verbose,
		  FILE *fp);

template<typename T>
void update_parents_list(list <int>& parents,
			 const int begin, const int end, 
			 SquareBlockMatrix<T>* mtrx);
