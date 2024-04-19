// Dalibor Lukas, May 2023

using namespace std;

//
//#include "mex.h"
//#include "BEM3dLaplace.h"

#include "bem.quadrature.h"

//#define  M_PI  3.1415926535897932

namespace espreso {

static double *qPoints2, *qWeights2;
static double *qPoints3, *qWeights3;
static double *qPoints4, *qWeights4;
static int quadn;
static double *qPoints, *qWeights;
static int quadTrin;
static Point2D *qTriPoints;
static double *qTriWeights;


static double SingleLayerLaplace3d_idPanels (double *us, double *vs)
{
  int i;
  double u1, u2, u3, v1, v2, v3, J1, J2, J3, Jsq;
  double xy1, xy2, xy3, I, Ii, eta3;

  u1 = us[0]; u2 = us[1]; u3 = us[2];
  v1 = vs[0]; v2 = vs[1]; v3 = vs[2];
  J1 = u2*v3-u3*v2; J2 = u3*v1-u1*v3; J3 = u1*v2-u2*v1;
  Jsq = J1*J1+J2*J2+J3*J3;

  for (i=0, I=0.0; i<quadn; i++)
  {
    eta3 = qPoints[i];
    xy1 = eta3*u1+v1;
    xy2 = eta3*u2+v2;
    xy3 = eta3*u3+v3;
    Ii = 1.0/sqrt(xy1*xy1+xy2*xy2+xy3*xy3);
    xy1 = eta3*v1+u1;
    xy2 = eta3*v2+u2;
    xy3 = eta3*v3+u3;
    Ii += 1.0/sqrt(xy1*xy1+xy2*xy2+xy3*xy3);
    xy1 = eta3*(u1+v1)-v1;
    xy2 = eta3*(u2+v2)-v2;
    xy3 = eta3*(u3+v3)-v3;
    Ii += 1.0/sqrt(xy1*xy1+xy2*xy2+xy3*xy3);
    I += qWeights[i]*Ii;
  }

  return Jsq * I / M_PI / 12.0;
}


static double SingleLayerLaplace3d_commonEdge (double *us, double *vs, double *ws)
{
  int i, j;
  double u1, u2, u3, v1, v2, v3, w1, w2, w3, J1, J2, J3, J;
  double xy1, xy2, xy3, I, Iij, eta2, eta3, weta2, weta3;

  u1 = us[0]; u2 = us[1]; u3 = us[2];
  v1 = vs[0]; v2 = vs[1]; v3 = vs[2];
  w1 = ws[0]; w2 = ws[1]; w3 = ws[2];
  J1 = u2*v3-u3*v2; J2 = u3*v1-u1*v3; J3 = u1*v2-u2*v1;
  J = sqrt(J1*J1+J2*J2+J3*J3);
  J1 = u2*w3-u3*w2; J2 = u3*w1-u1*w3; J3 = u1*w2-u2*w1;
  J *= sqrt(J1*J1+J2*J2+J3*J3);
  for (i=0, I=0.0; i<quadn; i++)
  {
    eta2 = qPoints[i];
    weta2 = qWeights[i];
    for (j=0; j<quadn; j++)
    {
      eta3 = qPoints[j];
      weta3 = qWeights[j];
      xy1 = eta2*(eta3*(u1+w1)-w1)+v1;
      xy2 = eta2*(eta3*(u2+w2)-w2)+v2;
      xy3 = eta2*(eta3*(u3+w3)-w3)+v3;
      Iij = 1.0/sqrt(xy1*xy1+xy2*xy2+xy3*xy3);
      xy1 = -eta2*(u1+v1+eta3*w1)+v1;
      xy2 = -eta2*(u2+v2+eta3*w2)+v2;
      xy3 = -eta2*(u3+v3+eta3*w3)+v3;
      Iij += 1.0/sqrt(xy1*xy1+xy2*xy2+xy3*xy3);
      xy1 = eta2*(v1-eta3*(u1+v1))-w1;
      xy2 = eta2*(v2-eta3*(u2+v2))-w2;
      xy3 = eta2*(v3-eta3*(u3+v3))-w3;
      Iij += 1.0/sqrt(xy1*xy1+xy2*xy2+xy3*xy3);
      xy1 = -eta2*(w1+eta3*(u1+v1))+v1;
      xy2 = -eta2*(w2+eta3*(u2+v2))+v2;
      xy3 = -eta2*(w3+eta3*(u3+v3))+v3;
      Iij += 1.0/sqrt(xy1*xy1+xy2*xy2+xy3*xy3);
      Iij *= eta2;
      xy1 = eta2*(u1+w1)+eta3*v1-w1;
      xy2 = eta2*(u2+w2)+eta3*v2-w2;
      xy3 = eta2*(u3+w3)+eta3*v3-w3;
      Iij += 1.0/sqrt(xy1*xy1+xy2*xy2+xy3*xy3);
      I += weta2*weta3*Iij;
    }
  }

  return J * I / M_PI / 24.0;
}


static double SingleLayerLaplace3d_commonVertex (double *us, double *vs, double *ws, double *zs)
{
  int i, j, k;
  double u1, u2, u3, v1, v2, v3, w1, w2, w3, z1, z2, z3, J1, J2, J3, J;
  double xy1, xy2, xy3, I, Iijk, eta1, eta2, eta3, weta1, weta2, weta3;

  u1 = us[0]; u2 = us[1]; u3 = us[2];
  v1 = vs[0]; v2 = vs[1]; v3 = vs[2];
  w1 = ws[0]; w2 = ws[1]; w3 = ws[2];
  z1 = zs[0]; z2 = zs[1]; z3 = zs[2];
  J1 = u2*v3-u3*v2; J2 = u3*v1-u1*v3; J3 = u1*v2-u2*v1;
  J = sqrt(J1*J1+J2*J2+J3*J3);
  J1 = w2*z3-w3*z2; J2 = w3*z1-w1*z3; J3 = w1*z2-w2*z1;
  J *= sqrt(J1*J1+J2*J2+J3*J3);

  for (i=0, I=0.0; i<quadn; i++)
  {
    eta1 = qPoints[i];
    weta1 = qWeights[i];
    for (j=0; j<quadn; j++)
    {
      eta2 = qPoints[j];
      weta2 = qWeights[j];
      for (k=0; k<quadn; k++)
      {
        eta3 = qPoints[k];
    weta3 = qWeights[k];
    xy1 = u1+eta1*v1-eta2*(w1+eta3*z1);
    xy2 = u2+eta1*v2-eta2*(w2+eta3*z2);
    xy3 = u3+eta1*v3-eta2*(w3+eta3*z3);
        Iijk = 1.0/sqrt(xy1*xy1+xy2*xy2+xy3*xy3);
    xy1 = w1+eta1*z1-eta2*(u1+eta3*v1);
    xy2 = w2+eta1*z2-eta2*(u2+eta3*v2);
    xy3 = w3+eta1*z3-eta2*(u3+eta3*v3);
        Iijk += 1.0/sqrt(xy1*xy1+xy2*xy2+xy3*xy3);
    Iijk *= eta2;
    I += weta1*weta2*weta3*Iijk;
      }
    }
  }

  return J * I / M_PI / 12.0;
}


static double SingleLayerLaplace3d_disjointPanels
(double *As, double *us, double *vs, double *Bs, double *ws, double *zs)
{
  int i, j;
  double A1, A2, A3, u1, u2, u3, v1, v2, v3, B1, B2, B3, w1, w2, w3, z1, z2, z3;
  double J1, J2, J3, J;
  double x1, x2, x3, xy1, xy2, xy3, I, ksi1, ksi2, eta1, eta2, wksi, weta;

  A1 = As[0]; A2 = As[1]; A3 = As[2];
  u1 = us[0]; u2 = us[1]; u3 = us[2];
  v1 = vs[0]; v2 = vs[1]; v3 = vs[2];
  B1 = Bs[0]; B2 = Bs[1]; B3 = Bs[2];
  w1 = ws[0]; w2 = ws[1]; w3 = ws[2];
  z1 = zs[0]; z2 = zs[1]; z3 = zs[2];
  J1 = u2*v3-u3*v2; J2 = u3*v1-u1*v3; J3 = u1*v2-u2*v1;
  J = sqrt(J1*J1+J2*J2+J3*J3);
  J1 = w2*z3-w3*z2; J2 = w3*z1-w1*z3; J3 = w1*z2-w2*z1;
  J *= sqrt(J1*J1+J2*J2+J3*J3);

  for (i=0, I=0.0; i<quadTrin; i++)
  {
    ksi1 = qTriPoints[i].x;
    ksi2 = qTriPoints[i].y;
    wksi = qTriWeights[i];
    x1 = A1+ksi1*u1+ksi2*v1;
    x2 = A2+ksi1*u2+ksi2*v2;
    x3 = A3+ksi1*u3+ksi2*v3;
    for (j=0; j<quadTrin; j++)
    {
      eta1 = qTriPoints[j].x;
      eta2 = qTriPoints[j].y;
      weta = qTriWeights[j];
      xy1 = x1-B1-eta1*w1-eta2*z1;
      xy2 = x2-B2-eta1*w2-eta2*z2;
      xy3 = x3-B3-eta1*w3-eta2*z3;

      I += wksi*weta/sqrt(xy1*xy1+xy2*xy2+xy3*xy3);
    }
  }

  return J * I / M_PI / 4.0;
}


static void DoubleLayerLaplace3d_commonEdge(double *us, double *vs, double *ws,
                     double *Kij)
{
  int i, j;
  double u1, u2, u3, v1, v2, v3, w1, w2, w3, J1, J2, J3, J, JT, n1, n2, n3;
  double xy1, xy2, xy3, normxy, xyn, H, I1, I2, I3, I1ij, I2ij, I3ij;
  double eta1, eta2, eta3, weta2, weta3, eta12, eta123;

  u1 = us[0]; u2 = us[1]; u3 = us[2];
  v1 = vs[0]; v2 = vs[1]; v3 = vs[2];
  w1 = ws[0]; w2 = ws[1]; w3 = ws[2];
  J1 = u2*v3-u3*v2; J2 = u3*v1-u1*v3; J3 = u1*v2-u2*v1;
  J = sqrt(J1*J1+J2*J2+J3*J3);
  J1 = u2*w3-u3*w2; J2 = u3*w1-u1*w3; J3 = u1*w2-u2*w1;
  JT = sqrt(J1*J1+J2*J2+J3*J3);
  n1 = J1/JT; n2 = J2/JT; n3 = J3/JT;
  J *= JT;

  eta1 = 0.5;
  for (i=0, I1=0.0, I2=0.0, I3=0.0; i<quadn; i++)
  {
    eta2 = qPoints[i];
    weta2 = qWeights[i];
    eta12 = eta1*eta2;
    for (j=0; j<quadn; j++)
    {
      eta3 = qPoints[j];
      weta3 = qWeights[j];
      eta123 = eta12*eta3;
      xy1 = eta2*(eta3*(u1+w1)-w1)+v1;
      xy2 = eta2*(eta3*(u2+w2)-w2)+v2;
      xy3 = eta2*(eta3*(u3+w3)-w3)+v3;
      normxy = sqrt(xy1*xy1+xy2*xy2+xy3*xy3);
      xyn = xy1*n1+xy2*n2+xy3*n3;
      H = xyn/normxy/normxy/normxy;
      I1ij = qWeights2[0]*qPoints2[0]*(1.0-qPoints2[0]*(1.0-eta123));
      I1ij+= qWeights2[1]*qPoints2[1]*(1.0-qPoints2[1]*(1.0-eta123));
      I1ij *= H;
      I2ij = qWeights2[0]*qPoints2[0] * qPoints2[0]*(1.0-eta12);
      I2ij += qWeights2[1]*qPoints2[1] * qPoints2[1]*(1.0-eta12);
      I2ij *= H;
      I3ij = qWeights2[0]*qPoints2[0] * qPoints2[0]*eta12*(1.0-eta3);
      I3ij += qWeights2[1]*qPoints2[1] * qPoints2[1]*eta12*(1.0-eta3);
      I3ij *= H;

      xy1 = -eta2*(u1+v1+eta3*w1)+v1;
      xy2 = -eta2*(u2+v2+eta3*w2)+v2;
      xy3 = -eta2*(u3+v3+eta3*w3)+v3;
      normxy = sqrt(xy1*xy1+xy2*xy2+xy3*xy3);
      xyn = xy1*n1+xy2*n2+xy3*n3;
      H = xyn/normxy/normxy/normxy;
      I1ij += qWeights2[0]*qPoints2[0] * (1.0-qPoints2[0]) * H;
      I1ij += qWeights2[1]*qPoints2[1] * (1.0-qPoints2[1]) * H;
      I2ij += qWeights2[0]*qPoints2[0] * qPoints2[0]*(1.0-eta123) * H;
      I2ij += qWeights2[1]*qPoints2[1] * qPoints2[1]*(1.0-eta123) * H;
      I3ij += qWeights2[0]*qPoints2[0] * qPoints2[0]*eta123 * H;
      I3ij += qWeights2[1]*qPoints2[1] * qPoints2[1]*eta123 * H;

      xy1 = eta2*(v1-eta3*(u1+v1))-w1;
      xy2 = eta2*(v2-eta3*(u2+v2))-w2;
      xy3 = eta2*(v3-eta3*(u3+v3))-w3;
      normxy = sqrt(xy1*xy1+xy2*xy2+xy3*xy3);
      xyn = xy1*n1+xy2*n2+xy3*n3;
      H = xyn/normxy/normxy/normxy;
      I1ij += qWeights2[0]*qPoints2[0] * (1.0-qPoints2[0]) * H;
      I1ij += qWeights2[1]*qPoints2[1] * (1.0-qPoints2[1]) * H;
      I2ij += qWeights2[0]*qPoints2[0] * qPoints2[0]*(1.0-eta1) * H;
      I2ij += qWeights2[1]*qPoints2[1] * qPoints2[1]*(1.0-eta1) * H;
      I3ij += qWeights2[0]*qPoints2[0] * qPoints2[0]*eta1 * H;
      I3ij += qWeights2[1]*qPoints2[1] * qPoints2[1]*eta1 * H;

      xy1 = -eta2*(w1+eta3*(u1+v1))+v1;
      xy2 = -eta2*(w2+eta3*(u2+v2))+v2;
      xy3 = -eta2*(w3+eta3*(u3+v3))+v3;
      normxy = sqrt(xy1*xy1+xy2*xy2+xy3*xy3);
      xyn = xy1*n1+xy2*n2+xy3*n3;
      H = xyn/normxy/normxy/normxy;
      I1ij += qWeights2[0]*qPoints2[0] * (1.0-qPoints2[0]) * H;
      I1ij += qWeights2[1]*qPoints2[1] * (1.0-qPoints2[1]) * H;
      I2ij += qWeights2[0]*qPoints2[0] * qPoints2[0]*(1.0-eta12) * H;
      I2ij += qWeights2[1]*qPoints2[1] * qPoints2[1]*(1.0-eta12) * H;
      I3ij += qWeights2[0]*qPoints2[0] * qPoints2[0]*eta12 * H;
      I3ij += qWeights2[1]*qPoints2[1] * qPoints2[1]*eta12 * H;

      I1ij *= eta2;
      I2ij *= eta2;
      I3ij *= eta2;

      xy1 = eta2*(u1+w1)+eta3*v1-w1;
      xy2 = eta2*(u2+w2)+eta3*v2-w2;
      xy3 = eta2*(u3+w3)+eta3*v3-w3;
      normxy = sqrt(xy1*xy1+xy2*xy2+xy3*xy3);
      xyn = xy1*n1+xy2*n2+xy3*n3;
      H = xyn/normxy/normxy/normxy;
      I1ij += qWeights2[0]*qPoints2[0]*(1.0-qPoints2[0]*(1.0-eta12))
    * H;
      I1ij += qWeights2[1]*qPoints2[1]*(1.0-qPoints2[1]*(1.0-eta12))
    * H;
      I2ij += qWeights2[0]*qPoints2[0] * qPoints2[0]*(1.0-eta1) * H;
      I2ij += qWeights2[1]*qPoints2[1] * qPoints2[1]*(1.0-eta1) * H;
      I3ij += qWeights2[0]*qPoints2[0] * qPoints2[0]*eta1*(1.0-eta2)
    * H;
      I3ij += qWeights2[1]*qPoints2[1] * qPoints2[1]*eta1*(1.0-eta2)
    * H;

      I1 += weta2*weta3*I1ij;
      I2 += weta2*weta3*I2ij;
      I3 += weta2*weta3*I3ij;
    }
  }

  Kij[0] = J * I1 / 4.0 / M_PI;
  Kij[1] = J * I2 / 4.0 / M_PI;
  Kij[2] = J * I3 / 4.0 / M_PI;
}


static void DoubleLayerLaplace3d_commonVertex(double *us, double *vs, double *ws,
                       double *zs, double *Kij)
{
  int i, j, k;
  double u1, u2, u3, v1, v2, v3, w1, w2, w3, z1, z2, z3;
  double J1, J2, J3, J, JT, n1, n2, n3;
  double xy1, xy2, xy3, normxy, xyn, H, I1, I2, I3, I1ijk, I2ijk, I3ijk;
  double eta1, eta2, eta3, weta1, weta2, weta3;

  u1 = us[0]; u2 = us[1]; u3 = us[2];
  v1 = vs[0]; v2 = vs[1]; v3 = vs[2];
  w1 = ws[0]; w2 = ws[1]; w3 = ws[2];
  z1 = zs[0]; z2 = zs[1]; z3 = zs[2];
  J1 = u2*v3-u3*v2; J2 = u3*v1-u1*v3; J3 = u1*v2-u2*v1;
  J = sqrt(J1*J1+J2*J2+J3*J3);
  J1 = w2*z3-w3*z2; J2 = w3*z1-w1*z3; J3 = w1*z2-w2*z1;
  JT = sqrt(J1*J1+J2*J2+J3*J3);
  n1 = J1/JT; n2 = J2/JT; n3 = J3/JT;
  J *= JT;

  for (i=0, I1=0.0, I2=0.0, I3=0.0; i<quadn; i++)
  {
    eta1 = qPoints[i];
    weta1 = qWeights[i];
    for (j=0; j<quadn; j++)
    {
      eta2 = qPoints[j];
      weta2 = qWeights[j];
      for (k=0; k<quadn; k++)
      {
        eta3 = qPoints[k];
    weta3 = qWeights[k];
    xy1 = u1+eta1*v1-eta2*(w1+eta3*z1);
    xy2 = u2+eta1*v2-eta2*(w2+eta3*z2);
    xy3 = u3+eta1*v3-eta2*(w3+eta3*z3);
    normxy = sqrt(xy1*xy1+xy2*xy2+xy3*xy3);
    xyn = xy1*n1+xy2*n2+xy3*n3;
    H = xyn/normxy/normxy/normxy;
        I1ijk = qWeights2[0]*qPoints2[0]*(1.0-qPoints2[0]*eta2);
        I1ijk += qWeights2[1]*qPoints2[1]*(1.0-qPoints2[1]*eta2);
    I1ijk *= H;
        I2ijk = qWeights2[0]*qPoints2[0]*qPoints2[0]*eta2*(1.0-eta3);
        I2ijk += qWeights2[1]*qPoints2[1]*qPoints2[1]*eta2*(1.0-eta3);
    I2ijk *= H;
        I3ijk = qWeights2[0]*qPoints2[0]*qPoints2[0]*eta2*eta3;
        I3ijk += qWeights2[1]*qPoints2[1]*qPoints2[1]*eta2*eta3;
    I3ijk *= H;

    xy1 = eta2*(u1+eta3*v1)-(w1+eta1*z1);
    xy2 = eta2*(u2+eta3*v2)-(w2+eta1*z2);
    xy3 = eta2*(u3+eta3*v3)-(w3+eta1*z3);
    normxy = sqrt(xy1*xy1+xy2*xy2+xy3*xy3);
    xyn = xy1*n1+xy2*n2+xy3*n3;
    H = xyn/normxy/normxy/normxy;
        I1ijk += qWeights2[0]*qPoints2[0]*(1.0-qPoints2[0]) * H;
        I1ijk += qWeights2[1]*qPoints2[1]*(1.0-qPoints2[1]) * H;
        I2ijk += qWeights2[0]*qPoints2[0]*qPoints2[0]*(1.0-eta1) * H;
        I2ijk += qWeights2[1]*qPoints2[1]*qPoints2[1]*(1.0-eta1) * H;
        I3ijk += qWeights2[0]*qPoints2[0]*qPoints2[0]*eta1 * H;
        I3ijk += qWeights2[1]*qPoints2[1]*qPoints2[1]*eta1 * H;

    I1 += weta1*weta2*weta3 * eta2 * I1ijk;
    I2 += weta1*weta2*weta3 * eta2 * I2ijk;
    I3 += weta1*weta2*weta3 * eta2 * I3ijk;
      }
    }
  }

  Kij[0] = J * I1 / 4.0 / M_PI;
  Kij[1] = J * I2 / 4.0 / M_PI;
  Kij[2] = J * I3 / 4.0 / M_PI;
}


static void DoubleLayerLaplace3d_disjointPanels(double *As, double *us, double *vs,
                     double *Bs, double *ws, double *zs,
                     double *Kij)
{
  int i, j;
  double A1, A2, A3, u1, u2, u3, v1, v2, v3, B1, B2, B3, w1, w2, w3, z1, z2, z3;
  double J1, J2, J3, J, JT, n1, n2, n3;
  double x1, x2, x3, xy1, xy2, xy3, normxy, xyn, H, I1, I2, I3;
  double ksi1, ksi2, eta1, eta2, wksi, weta;

  A1 = As[0]; A2 = As[1]; A3 = As[2];
  u1 = us[0]; u2 = us[1]; u3 = us[2];
  v1 = vs[0]; v2 = vs[1]; v3 = vs[2];
  B1 = Bs[0]; B2 = Bs[1]; B3 = Bs[2];
  w1 = ws[0]; w2 = ws[1]; w3 = ws[2];
  z1 = zs[0]; z2 = zs[1]; z3 = zs[2];
  J1 = u2*v3-u3*v2; J2 = u3*v1-u1*v3; J3 = u1*v2-u2*v1;
  J = sqrt(J1*J1+J2*J2+J3*J3);
  J1 = w2*z3-w3*z2; J2 = w3*z1-w1*z3; J3 = w1*z2-w2*z1;
  JT = sqrt(J1*J1+J2*J2+J3*J3);
  n1 = J1/JT; n2 = J2/JT; n3 = J3/JT;
  J *= JT;

  for (i=0, I1=0.0, I2=0.0, I3=0.0; i<quadTrin; i++)
  {
    ksi1 = qTriPoints[i].x;
    ksi2 = qTriPoints[i].y;
    wksi = qTriWeights[i];
    x1 = A1+ksi1*u1+ksi2*v1;
    x2 = A2+ksi1*u2+ksi2*v2;
    x3 = A3+ksi1*u3+ksi2*v3;
    for (j=0; j<quadTrin; j++)
    {
      eta1 = qTriPoints[j].x;
      eta2 = qTriPoints[j].y;
      weta = qTriWeights[j];
      xy1 = x1-B1-eta1*w1-eta2*z1;
      xy2 = x2-B2-eta1*w2-eta2*z2;
      xy3 = x3-B3-eta1*w3-eta2*z3;
      normxy = sqrt(xy1*xy1+xy2*xy2+xy3*xy3);
      xyn = xy1*n1+xy2*n2+xy3*n3;
      H = xyn/normxy/normxy/normxy;
      I1 += wksi*weta * (1.0-eta1-eta2) * H;
      I2 += wksi*weta * eta1 * H;
      I3 += wksi*weta * eta2 * H;
    }
  }

  Kij[0] = J * I1 / M_PI / 4.0;
  Kij[1] = J * I2 / M_PI / 4.0;
  Kij[2] = J * I3 / M_PI / 4.0;
}


static void BEM3dLaplace (esint np, double *points, esint ne, esint *elemNodes,
           int order, double *V, double *K, double *D, double *M)
{
  esint i, j, k, l, cmnIdxi[2],cmnIdxj[2],cmnIdxSize, remIdxi, remIdxj, idx,tmp;

  /*
  if (V==NULL)
    V = new double[ne*ne];
  */

  memset(V,0,ne*ne*sizeof(double));
  memset(K,0,ne*np*sizeof(double));
  memset(D,0,np*np*sizeof(double));
  memset(M,0,ne*np*sizeof(double));

  quadn = LineQuadratureSizes[order];
  qPoints = LineQuadraturePoints[order];
  qWeights = LineQuadratureWeights[order];
  quadTrin = TriangleQuadratureSizes[order];
  qTriPoints = TriangleQuadraturePoints[order];
  qTriWeights = TriangleQuadratureWeights[order];
  qPoints2 = LineQuadraturePoints[1];
  qWeights2 = LineQuadratureWeights[1];
  qPoints3 = LineQuadraturePoints[2];
  qWeights3 = LineQuadratureWeights[2];
  qPoints4 = LineQuadraturePoints[3];
  qWeights4 = LineQuadratureWeights[3];

  esint faceNodesi[3], faceNodesj[3];
  double Pi[3][3], Pj[3][3], ni[3], Ji, nj[3], Jj;
  double u[3], v[3], w[3], z[3];
  double Kij[3];
  double a11, a12, a22, detA, b11, b21, b12, b22, b13, b23;
  double y11, y21, y12, y22, y13, y23;
  double curli1[3], curli2[3], curli3[3], curlj1[3], curlj2[3], curlj3[3];

  for (i=0, idx=0; i<ne; i++)
    {
      for (k=0; k<3; k++)
    faceNodesi[k] = elemNodes[3*i+k] + 1;
      for (k=0; k<3; k++)
    memcpy(Pi[k],points+3*(faceNodesi[k]-1),3*sizeof(double));
      ni[0] = (Pi[1][1]-Pi[0][1])*(Pi[2][2]-Pi[0][2]) -
    (Pi[1][2]-Pi[0][2])*(Pi[2][1]-Pi[0][1]);
      ni[1] = (Pi[1][2]-Pi[0][2])*(Pi[2][0]-Pi[0][0]) -
    (Pi[1][0]-Pi[0][0])*(Pi[2][2]-Pi[0][2]);
      ni[2] = (Pi[1][0]-Pi[0][0])*(Pi[2][1]-Pi[0][1]) -
    (Pi[1][1]-Pi[0][1])*(Pi[2][0]-Pi[0][0]);
      Ji = sqrt(ni[0]*ni[0]+ni[1]*ni[1]+ni[2]*ni[2]);
      for (k=0; k<3; k++)
    M[i+(faceNodesi[k]-1)*ne] += Ji/6.0;
      a11 = a12 = a22 = 0;
      for (k=0; k<3; k++)
    {
      ni[k] /= Ji;
      a11 += (Pi[1][k]-Pi[0][k]) * (Pi[1][k]-Pi[0][k]);
      a22 += (Pi[2][k]-Pi[0][k]) * (Pi[2][k]-Pi[0][k]);
      a12 += (Pi[1][k]-Pi[0][k]) * (Pi[2][k]-Pi[0][k]);
    }
      detA = a11*a22-a12*a12;
      b11 = -(Pi[1][1]-Pi[0][1])*ni[2] + (Pi[1][2]-Pi[0][2])*ni[1];
      b21 = -(Pi[2][1]-Pi[0][1])*ni[2] + (Pi[2][2]-Pi[0][2])*ni[1];
      b12 = (Pi[1][0]-Pi[0][0])*ni[2] - (Pi[1][2]-Pi[0][2])*ni[0];
      b22 = (Pi[2][0]-Pi[0][0])*ni[2] - (Pi[2][2]-Pi[0][2])*ni[0];
      b13 = -(Pi[1][0]-Pi[0][0])*ni[1] + (Pi[1][1]-Pi[0][1])*ni[0];
      b23 = -(Pi[2][0]-Pi[0][0])*ni[1] + (Pi[2][1]-Pi[0][1])*ni[0];
      y11 = (b11*a22-b21*a12)/detA; y21 = (a11*b21-a12*b11)/detA;
      y12 = (b12*a22-b22*a12)/detA; y22 = (a11*b22-a12*b12)/detA;
      y13 = (b13*a22-b23*a12)/detA; y23 = (a11*b23-a12*b13)/detA;
      curli1[0] = -y11-y21; curli1[1] = y11; curli1[2] = y21;
      curli2[0] = -y12-y22; curli2[1] = y12; curli2[2] = y22;
      curli3[0] = -y13-y23; curli3[1] = y13; curli3[2] = y23;

    for (j=0; j<ne; j++, idx++) {
      if (j==i) {
        for (k=0; k<3; k++)
        {
          u[k] = Pi[1][k]-Pi[0][k];
          v[k] = Pi[2][k]-Pi[1][k];
          K[i+(faceNodesi[k]-1)*ne] += 0.0;
        }
          V[idx] = SingleLayerLaplace3d_idPanels(u,v);
          for (k=0; k<3; k++)
        for (l=0; l<3; l++)
          D[faceNodesi[k]-1+(faceNodesi[l]-1)*np] +=
            curli1[k]*V[idx]*curli1[l] + curli2[k]*V[idx]*curli2[l] +
            curli3[k]*V[idx]*curli3[l];
          continue;
        }

      for (k=0; k<3; k++)
        faceNodesj[k] = elemNodes[3*j+k] + 1;
      for (k=0; k<3; k++)
        memcpy(Pj[k],points+3*(faceNodesj[k]-1),3*sizeof(double));

      nj[0] = (Pj[1][1]-Pj[0][1])*(Pj[2][2]-Pj[0][2]) -
        (Pj[1][2]-Pj[0][2])*(Pj[2][1]-Pj[0][1]);
      nj[1] = (Pj[1][2]-Pj[0][2])*(Pj[2][0]-Pj[0][0]) -
        (Pj[1][0]-Pj[0][0])*(Pj[2][2]-Pj[0][2]);
      nj[2] = (Pj[1][0]-Pj[0][0])*(Pj[2][1]-Pj[0][1]) -
        (Pj[1][1]-Pj[0][1])*(Pj[2][0]-Pj[0][0]);
      Jj = sqrt(nj[0]*nj[0]+nj[1]*nj[1]+nj[2]*nj[2]);
      a11 = a12 = a22 = 0;
      for (k=0; k<3; k++)
        {
          nj[k] /= Jj;
          a11 += (Pj[1][k]-Pj[0][k]) * (Pj[1][k]-Pj[0][k]);
          a22 += (Pj[2][k]-Pj[0][k]) * (Pj[2][k]-Pj[0][k]);
          a12 += (Pj[1][k]-Pj[0][k]) * (Pj[2][k]-Pj[0][k]);
        }
      detA = a11*a22-a12*a12;
      b11 = -(Pj[1][1]-Pj[0][1])*nj[2] + (Pj[1][2]-Pj[0][2])*nj[1];
      b21 = -(Pj[2][1]-Pj[0][1])*nj[2] + (Pj[2][2]-Pj[0][2])*nj[1];
      b12 = (Pj[1][0]-Pj[0][0])*nj[2] - (Pj[1][2]-Pj[0][2])*nj[0];
      b22 = (Pj[2][0]-Pj[0][0])*nj[2] - (Pj[2][2]-Pj[0][2])*nj[0];
      b13 = -(Pj[1][0]-Pj[0][0])*nj[1] + (Pj[1][1]-Pj[0][1])*nj[0];
      b23 = -(Pj[2][0]-Pj[0][0])*nj[1] + (Pj[2][1]-Pj[0][1])*nj[0];
      y11 = (b11*a22-b21*a12)/detA; y21 = (a11*b21-a12*b11)/detA;
      y12 = (b12*a22-b22*a12)/detA; y22 = (a11*b22-a12*b12)/detA;
      y13 = (b13*a22-b23*a12)/detA; y23 = (a11*b23-a12*b13)/detA;
      curlj1[0] = -y11-y21; curlj1[1] = y11; curlj1[2] = y21;
      curlj2[0] = -y12-y22; curlj2[1] = y12; curlj2[2] = y22;
      curlj3[0] = -y13-y23; curlj3[1] = y13; curlj3[2] = y23;

      for (k=0, cmnIdxSize=0; k<3; k++)
        {
          l = Find<esint>(faceNodesj[k],faceNodesi,3);
          if (l>0)
        {
          cmnIdxi[cmnIdxSize] = l-1;
          cmnIdxj[cmnIdxSize] = k;
          cmnIdxSize++;
        }
        }
      if (cmnIdxSize==2)
        if ((cmnIdxj[0]==0 && cmnIdxj[1]==2) ||
        (cmnIdxj[0]==1 && cmnIdxj[1]==0) ||
        (cmnIdxj[0]==2 && cmnIdxj[1]==1))
          {
        tmp = cmnIdxj[0]; cmnIdxj[0] = cmnIdxj[1]; cmnIdxj[1] = tmp;
        tmp = cmnIdxi[0]; cmnIdxi[0] = cmnIdxi[1]; cmnIdxi[1] = tmp;
          }

      switch (cmnIdxSize)
        {
        case 0: // disjoint panels
          for (k=0; k<3; k++)
        {
          u[k] = Pi[1][k]-Pi[0][k];
          v[k] = Pi[2][k]-Pi[0][k];
          w[k] = Pj[1][k]-Pj[0][k];
          z[k] = Pj[2][k]-Pj[0][k];
        }
          V[idx] = SingleLayerLaplace3d_disjointPanels(Pi[0],u,v,
                               Pj[0],w,z);
          DoubleLayerLaplace3d_disjointPanels(Pi[0],u,v,Pj[0],w,z,Kij);
          for (k=0; k<3; k++)
        K[i+(faceNodesj[k]-1)*ne] += Kij[k];
          break;
        case 1: // common vertex
          for (k=0; k<3; k++)
        {
          u[k] = Pi[(cmnIdxi[0]+1)%3][k]-Pi[cmnIdxi[0]][k];
          v[k] = Pi[(cmnIdxi[0]+2)%3][k]-Pi[(cmnIdxi[0]+1)%3][k];
          w[k] = Pj[(cmnIdxj[0]+1)%3][k]-Pj[cmnIdxj[0]][k];
          z[k] = Pj[(cmnIdxj[0]+2)%3][k]-Pj[(cmnIdxj[0]+1)%3][k];
        }
          V[idx] = SingleLayerLaplace3d_commonVertex(u,v,w,z);
          DoubleLayerLaplace3d_commonVertex(u,v,w,z,Kij);
          for (k=0; k<3; k++)
        K[i+(faceNodesj[(cmnIdxj[0]+k)%3]-1)*ne] += Kij[k];
          break;
        case 2: // common edge
          for (k=0; k<3; k++)
        {
          if (Find<esint>(k,cmnIdxi,2)==0)
            remIdxi = k;
          if (Find<esint>(k,cmnIdxj,2)==0)
            remIdxj = k;
        }
          for (k=0; k<3; k++)
        {
          u[k] = Pj[cmnIdxj[1]][k]-Pj[cmnIdxj[0]][k];
          v[k] = Pi[remIdxi][k]-Pj[cmnIdxj[1]][k];
          w[k] = Pj[remIdxj][k]-Pj[cmnIdxj[1]][k];
        }
          V[idx] = SingleLayerLaplace3d_commonEdge(u,v,w);
          DoubleLayerLaplace3d_commonEdge(u,v,w,Kij);
          K[i+(faceNodesj[cmnIdxj[0]]-1)*ne] += Kij[0];
          K[i+(faceNodesj[cmnIdxj[1]]-1)*ne] += Kij[1];
          K[i+(faceNodesj[remIdxj]-1)*ne] += Kij[2];
          break;
        }

      for (k=0; k<3; k++)
        for (l=0; l<3; l++)
          D[faceNodesi[k]-1+(faceNodesj[l]-1)*np] +=
        curli1[k]*V[idx]*curlj1[l] + curli2[k]*V[idx]*curlj2[l] +
        curli3[k]*V[idx]*curlj3[l];
    }
    }
}


static void BEM3dLaplaceVolume (double *nodes, esint ne, esint *elemNodes,
             esint np, double *points, int order,
             double *t, double *u, double *Vt_Wu)
{
  esint i, j, k, faceNodesj[3];
  double ksi1, ksi2, wksi, A1, A2, A3, u1, u2, u3, v1, v2, v3, n1, n2, n3, J, I;
  double x1, x2, x3, xy1, xy2, xy3, tj, uj1, uj2, uj3, normxy;

  quadTrin = TriangleQuadratureSizes[order];
  qTriPoints = TriangleQuadraturePoints[order];
  qTriWeights = TriangleQuadratureWeights[order];

  memset(Vt_Wu,0,np*sizeof(double));
  for (j=0; j<ne; j++)
    {
      tj = t[j];
      for (k=0; k<3; k++)
    faceNodesj[k] = elemNodes[3*j+k] + 1;
      uj1 = u[faceNodesj[0]-1];
      uj2 = u[faceNodesj[1]-1];
      uj3 = u[faceNodesj[2]-1];
      A1 = nodes[3*(faceNodesj[0]-1)];
      A2 = nodes[3*(faceNodesj[0]-1)+1];
      A3 = nodes[3*(faceNodesj[0]-1)+2];
      u1 = nodes[3*(faceNodesj[1]-1)]-A1;
      u2 = nodes[3*(faceNodesj[1]-1)+1]-A2;
      u3 = nodes[3*(faceNodesj[1]-1)+2]-A3;
      v1 = nodes[3*(faceNodesj[2]-1)]-A1;
      v2 = nodes[3*(faceNodesj[2]-1)+1]-A2;
      v3 = nodes[3*(faceNodesj[2]-1)+2]-A3;
      n1 = u2*v3-u3*v2; n2 = u3*v1-u1*v3; n3 = u1*v2-u2*v1;
      J = sqrt(n1*n1+n2*n2+n3*n3);
      n1 /= J; n2 /= J; n3 /= J;
      for (i=0; i<np; i++)
    {
      x1 = points[3*i];
      x2 = points[3*i+1];
      x3 = points[3*i+2];
      for (k=0, I=0.0; k<quadTrin; k++)
        {
          ksi1 = qTriPoints[k].x;
          ksi2 = qTriPoints[k].y;
          wksi = qTriWeights[k];
          xy1 = x1 - A1 - u1*ksi1 - v1*ksi2;
          xy2 = x2 - A2 - u2*ksi1 - v2*ksi2;
          xy3 = x3 - A3 - u3*ksi1 - v3*ksi2;
          normxy = sqrt(xy1*xy1+xy2*xy2+xy3*xy3);
          I += wksi * tj / normxy;
          I -= wksi * (uj1+(uj2-uj1)*ksi1+(uj3-uj1)*ksi2) *
        (xy1*n1+xy2*n2+xy3*n3) / normxy / normxy / normxy;
        }
      Vt_Wu[i] += J * I / 4.0 / M_PI;
    }
    }
 }


//static void invertMatrix (esint n, double *M, double *invM)
//{
//  esint i, j, k;
//  double pivot, factor;
//
//  memset(invM,0,n*n*sizeof(double));
//  for (i=0; i<n; i++)
//    invM[i*n+i] = 1.0;
//
//  // forward elimination
//  for (i=0; i<n; i++)
//    {
//      pivot = M[i*n+i];
//      for (j=i+1; j<n; j++)
//    {
//      factor = M[i*n+j] / pivot;
//          M[i*n+j] = 0.0;
//      for (k=i+1; k<n; k++)
//        M[k*n+j] -= factor * M[k*n+i];
//      for (k=0; k<n; k++)
//        invM[k*n+j] -= factor * invM[k*n+i];
//    }
//    }
//
//  // backward substitution
//  for (i=n-1; i>=0; i--)
//    {
//      pivot = M[i*n+i];
//      M[i*n+i] = 1.0;
//      for (k=0; k<n; k++)
//    invM[k*n+i] /= pivot;
//      for (j=i-1; j>=0; j--)
//    {
//      factor = M[i*n+j];
//      M[i*n+j] = 0.0;
//      for (k=0; k<n; k++)
//        invM[k*n+j] -= factor * invM[k*n+i];
//    }
//    }
//}


static void cholesky (int n, double *V)
{
  double x, r;
  esint i, j, k;

  /* Loop over columns */
  for (j=0; j<n; j++)
    {
      /* i = j */
      x = V[j*n+j]; /* V_jj */
      for (k=0; k<j; k++)
    x -= V[j*n+k] * V[j*n+k]; /* L_jk L_jk */
      x = sqrt(x);
      V[j*n+j] = x; /* L_jj */
      r = 1.0 / x;

      /* i != j */
      for (i=j+1; i<n; i++)
    {
      x = V[i*n+j]; /* V_ij */
      for (k=0; k<j; k++)
        x -= V[i*n+k] * V[j*n+k]; /* L_ik L_ij */
      V[i*n+j] = x * r; /* L_ij = x / L_jj */
    }
    }
}


static void choleskySolve (int n, double *L, double *rhs, double *u, int m=1)
{
  esint i, j, k;
  double *b;
  b = new double[n*m];
  memcpy(b,rhs,n*m*sizeof(double));
  for (i=0; i<n; i++)
    {
      for (k=0; k<m; k++)
    u[i+k*n] = b[i+k*n] / L[i*n+i];
      for (j=i+1; j<n; j++)
    for (k=0; k<m; k++)
      b[j+k*n] -= u[i+k*n]*L[j*n+i];
    }
  for (i=n-1; i>=0; i--)
    {
      for (k=0; k<m; k++)
    u[i+k*n] = u[i+k*n] / L[i*n+i];
      for (j=0; j<i; j++)
    for (k=0; k<m; k++)
      u[j+k*n] -= u[i+k*n]*L[i*n+j];
    }
  delete [] b;
}


//static void assembleSteklovPoincare (esint nElements, esint nNodes, double *invV,
//                  double *KK, double *D, double *K)
//{
//  esint i, j, k, l, idx;
//  memcpy(K,D,nNodes*nNodes*sizeof(double));
//  for (j=0, idx=0; j<nNodes; j++) // cols of K
//    for (i=0; i<nNodes; i++, idx++) // rows of K
//      for (k=0; k<nElements; k++) // rows of invV
//    for (l=0; l<nElements; l++) // cols of invV
//      K[idx] += KK[i*nElements+k] * invV[l*nElements+k] * KK[j*nElements+l];
//}


//void mexPrintfMatrix (int n, int m, double *mat, char *name=0)
//{
//  esint i, j;
//  if (name)
//    mexPrintf("%s = [",name);
//  else
//    mexPrintf("Mat = [");
//  for (i=0; i<n; i++)
//    {
//      for (j=0; j<m; j++)
//    mexPrintf("%1.5g ",mat[j*n+i]);
//      mexPrintf("; ");
//    }
//  mexPrintf("];\n");
//}

static void assembleSteklovPoincareCholesky (esint nElements, esint nNodes, double *V,
                      double *KK, double *D, double *K)
{
  esint i, j, k, idx;
  double *tmp;
  tmp = new double[nElements*nNodes];
  choleskySolve(nElements,V,KK,tmp,nNodes); // tmp := V \ KK
  // K := D + KK' * tmp
  for (j=0, idx=0; j<nNodes; j++) // cols of K
    for (i=0; i<nNodes; i++, idx++) // rows of K
      {
    K[idx] = D[idx];
    for (k=0; k<nElements; k++)
      K[idx] += KK[i*nElements+k] * tmp[j*nElements+k];
      }
  delete [] tmp;
}


//static void deleteMyBEMData (MyBEMData* &bem)
//{
//  bem->nNodes = 0;
//  if (bem->S)
//    {
//      delete [] bem->S; bem->S = 0;
//    }
//}


void BEM3DLaplace (double *K, esint np, double *points, esint ne, esint *elements, double conductivity)
{
  esint i, j, idx;
  double *V, *KK, *D, *M;
  /*
  deleteMyBEMData(bem);
  bem->nNodes = nNodes;
  */
  V = new double[ne*ne];
  KK = new double[ne*np];
  D = new double[np*np];
  M = new double[ne*np];
  BEM3dLaplace(np,points,ne,elements,7,V,KK,D,M);
  /*
  mexPrintfMatrix(nElements,nElements,V,"V");
  mexPrintfMatrix(nElements,nNodes,KK,"K");
  mexPrintfMatrix(nElements,nNodes,M,"M");
  mexPrintfMatrix(nNodes,nNodes,D,"D");
  */
  //invV = new double[nElements*nElements];
  //invertMatrix(nElements,V,invV);
  for (i=0, idx=0; i<ne; i++)
    for (j=0; j<np; j++, idx++)
      KK[idx] += 0.5*M[idx];
  //K = new double[nNodes*nNodes];
  //assembleSteklovPoincare(nElements,nNodes,invV,KK,D,K);
  ///*
  cholesky(ne,V);
  assembleSteklovPoincareCholesky(ne,np,V,KK,D,K);
  //*/
  for (i=0, idx=0; i<np; i++)
    for (j=0; j<np; j++, idx++)
      K[idx] *= conductivity;
  //bem->S = K;
  delete [] V;
  //delete [] invV;
  delete [] KK;
  delete [] D;
  delete [] M;
}


void BEM3DLaplaceEval(double *results, esint np, double *points, esint ne, esint *elements, esint ni, double *inner, double conductivity, double *dirichlet)
{
  esint i, j, idx;
  double *rhs, *neumann;
  double *V, *KK, *D, *M;
  /*
  deleteMyBEMData(bem);
  bem->nNodes = nNodes;
  */
  V = new double[ne*ne];
  KK = new double[ne*np];
  D = new double[np*np];
  M = new double[ne*np];
  BEM3dLaplace(np,points,ne,elements,7,V,KK,D,M);
  for (i=0, idx=0; i<ne; i++)
    for (j=0; j<np; j++, idx++)
      KK[idx] += 0.5*M[idx];

  rhs = new double[ne];
  memset(rhs,0,ne*sizeof(double));
  for (j=0, idx=0; j<np; j++)
    for (i=0; i<ne; i++, idx++)
      rhs[i] += KK[idx]*dirichlet[j];
  neumann = new double[ne];
  cholesky(ne,V);
  choleskySolve(ne,V,rhs,neumann);
  for (i=0; i<ne; i++)
    neumann[i] *= conductivity;
  BEM3dLaplaceVolume(points,ne,elements,ni,inner,7,neumann,dirichlet,results);
  delete [] neumann;
  delete [] rhs;
  delete [] V;
  delete [] KK;
  delete [] D;
  delete [] M;
}


//void mexFunction(int nlhs, mxArray *plhs[],
//                 int nrhs, const mxArray *prhs[] )
//{
//  int n, m, nP;
//  double *elemNodes;
//  double *nodes;
//  double *points;
//  double *dirichlet;
//  double *results;
//  int order;
//  double *V, *K, *D, *M;
//
//  n = mxGetN(prhs[0]);
//  nodes = mxGetPr(prhs[0]);
//
//  m = mxGetN(prhs[1]);
//  elemNodes = mxGetPr(prhs[1]);
//
//  /*
//  order = mxGetScalar(prhs[2]);
//
//  plhs[0] = mxCreateDoubleMatrix(m,m,mxREAL);
//  V = mxGetPr(plhs[0]);
//  plhs[1] = mxCreateDoubleMatrix(m,n,mxREAL);
//  K = mxGetPr(plhs[1]);
//  plhs[2] = mxCreateDoubleMatrix(n,n,mxREAL);
//  D = mxGetPr(plhs[2]);
//  plhs[3] = mxCreateDoubleMatrix(m,n,mxREAL);
//  M = mxGetPr(plhs[3]);
//  BEM3dLaplace(n,nodes,m,elemNodes,order,V,K,D,M);
//  */
//
//  /*
//  plhs[0] = mxCreateDoubleMatrix(n,n,mxREAL);
//  K = mxGetPr(plhs[0]);
//  getLaplace ( K, n, nodes, m, elemNodes, 1.0);
//  */
//
//  ///*
//  nP = mxGetN(prhs[2]);
//  points = mxGetPr(prhs[2]);
//
//  dirichlet = mxGetPr(prhs[3]);
//
//  plhs[0] = mxCreateDoubleMatrix(nP,1,mxREAL);
//  results = mxGetPr(plhs[0]);
//  evaluateLaplace( results, n, nodes, m, elemNodes, nP, points, 1.0, dirichlet);
//  //*/
//
//  return;
//}

}
