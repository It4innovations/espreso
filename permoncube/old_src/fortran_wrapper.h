#ifndef FORTRAN_WRAPPER_H_
#define FORTRAN_WRAPPER_H_


extern "C" {
  void userpl_ (
      const int &elem,
      const int &intpt,
      const int &mat,
      const int &ncomp,
      const int &kfirst,
      const int &kfsteq,
      const double &e,
      const double &nu,
      const double &dens,
      double *prop,
      const double *d,
      const int &ktform,
      const double &timval,
      const double &timinc,
      const double &tem,
      const double &dtem,
      const double &toffst,
      const double &flu,
      const double &dflu,
      double *epel,
      double *eppl,
      double *statev,
      double *usvr,
      double &epeq,
      double &plwork,
      double &sigepl,
      double &sigrat,
      double &depeq,
      double *dt);
}

#endif /* FORTRAN_WRAPPER_H_ */
