#ifndef FORTRAN_H_
#define FORTRAN_H_

//
// flag_fortran


/* FOTRAN lib.
 * 0  - - - library - off - - -
 * 1  MPI_Fint                            main
 * 2  init_general_                       main
 * 3  fe2feti_init_symbolic_              fe_symbolic_form_matrix
 * 4  fe2feti_symbolic_map_element_       fe_symbolic_form_matrix (in loop)
 * 5  fe2feti_symbolic_map_global_bc_     fe_symbolic_form_matrix
 * 6  fe2feti_symbolic_finalize_          fe_symbolic_form_matrix
 * 7  fe2feti_symbolic_factorize_         fe_symbolic_form_matrix
 * 8  fe2feti_init_numeric_               fe_symbolic_form_matrix ???
 * 9  fe2feti_numeric_map_element_           stima_subdomain         (in loop)
 * 10 fe2feti_numeric_factorize_          stima_subdomain
 * 11 fe2feti_numeric_map_global_bc_      definitionAndSolving
 * 12 fe2feti_map_global_rhs_             definitionAndSolving
 * 13 fe2feti_solve_                      definitionAndSolving
 * 14 fe2feti_finalize_                   main
 * 15 fe2feti_free_data_                  main
 */


extern "C" {
  void init_general_(const int &c_com);
  void fe2feti_init_symbolic_(const int &gneq, const int &n_elmt, const int &connectivity, const int &real_or_complex);
  void fe2feti_symbolic_map_element_(const int &neq, int *ieqs, const int &symmetric);
  void fe2feti_numeric_map_element_(const int &neq, int *ieqs ,double *element_matrix);
  void fin_general_();
  void fe2feti_symbolic_finalize_();
  void fe2feti_symbolic_factorize_();
  void fe2feti_numeric_factorize_();
  void fe2feti_symbolic_map_global_bc_(const int &neq,int *ig);
  void fe2feti_numeric_map_global_bc_(const int &neq,int *ig,double *vg);
  void fe2feti_map_global_rhs_(const int &neq,double *vg);
  void fe2feti_solve_();
  void fe2feti_numeric_factorize_();
  void fe2feti_init_numeric_(const int &real_or_complex);
  void fe2feti_gather_solution_(double *u, double *v);
  void fe2feti_finalize_();
  void fe2feti_free_data_();
}

#endif /* FORTRAN_H_ */
