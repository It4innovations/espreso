
#include "w.bem.h"

// Dalibor Lukas, May 2023

using namespace std;

#include "bem.quadrature.h"

namespace espreso {

//#define  M_PI  3.1415926535897932

static double *qPoints2, *qWeights2;
static double *qPoints3, *qWeights3;
static double *qPoints4, *qWeights4;
static int quadn;
static double *qPoints, *qWeights;
static int quadTrin;
static Point2D *qTriPoints;
static double *qTriWeights;

static void SingleLayerElasticity3d_idPanels (double *us, double *vs, double &I,
                       double &I11, double &I22, double &I33,
                       double &I12, double &I13, double &I23)
{
  int i;
  double u1, u2, u3, v1, v2, v3, J1, J2, J3, Jsq;
  double xy1, xy2, xy3, normxy, normxy3, eta3;
  double Ii, I11i, I22i, I33i, I12i, I13i, I23i;

  u1 = us[0]; u2 = us[1]; u3 = us[2];
  v1 = vs[0]; v2 = vs[1]; v3 = vs[2];
  J1 = u2*v3-u3*v2; J2 = u3*v1-u1*v3; J3 = u1*v2-u2*v1;
  Jsq = J1*J1+J2*J2+J3*J3;

  I = I11 = I22 = I33 = I12 = I13 = I23 = 0.0;
  for (i=0; i<quadn; i++)
  {
    eta3 = qPoints[i];
    xy1 = eta3*u1+v1;
    xy2 = eta3*u2+v2;
    xy3 = eta3*u3+v3;
    normxy = sqrt(xy1*xy1+xy2*xy2+xy3*xy3);
    normxy3 = normxy*normxy*normxy;
    Ii = 1.0/normxy;
    I11i = xy1*xy1/normxy3;
    I22i = xy2*xy2/normxy3;
    I33i = xy3*xy3/normxy3;
    I12i = xy1*xy2/normxy3;
    I13i = xy1*xy3/normxy3;
    I23i = xy2*xy3/normxy3;
    xy1 = eta3*v1+u1;
    xy2 = eta3*v2+u2;
    xy3 = eta3*v3+u3;
    normxy = sqrt(xy1*xy1+xy2*xy2+xy3*xy3);
    normxy3 = normxy*normxy*normxy;
    Ii += 1.0/normxy;
    I11i += xy1*xy1/normxy3;
    I22i += xy2*xy2/normxy3;
    I33i += xy3*xy3/normxy3;
    I12i += xy1*xy2/normxy3;
    I13i += xy1*xy3/normxy3;
    I23i += xy2*xy3/normxy3;
    xy1 = eta3*(u1+v1)-v1;
    xy2 = eta3*(u2+v2)-v2;
    xy3 = eta3*(u3+v3)-v3;
    normxy = sqrt(xy1*xy1+xy2*xy2+xy3*xy3);
    normxy3 = normxy*normxy*normxy;
    Ii += 1.0/normxy;
    I11i += xy1*xy1/normxy3;
    I22i += xy2*xy2/normxy3;
    I33i += xy3*xy3/normxy3;
    I12i += xy1*xy2/normxy3;
    I13i += xy1*xy3/normxy3;
    I23i += xy2*xy3/normxy3;
    I += qWeights[i]*Ii;
    I11 += qWeights[i]*I11i;
    I22 += qWeights[i]*I22i;
    I33 += qWeights[i]*I33i;
    I12 += qWeights[i]*I12i;
    I13 += qWeights[i]*I13i;
    I23 += qWeights[i]*I23i;
  }

  I *= Jsq / M_PI / 12.0;
  I11 *= Jsq / M_PI / 12.0;
  I22 *= Jsq / M_PI / 12.0;
  I33 *= Jsq / M_PI / 12.0;
  I12 *= Jsq / M_PI / 12.0;
  I13 *= Jsq / M_PI / 12.0;
  I23 *= Jsq / M_PI / 12.0;
}


static void SingleLayerElasticity3d_commonEdge (double *us, double *vs, double *ws,
                       double &I, double &I11, double &I22,
                       double &I33, double &I12,
                       double &I13, double &I23)
{
  int i, j;
  double u1, u2, u3, v1, v2, v3, w1, w2, w3, J1, J2, J3, J;
  double xy1, xy2, xy3, normxy, normxy3, eta2, eta3, weta2, weta3;
  double Iij, I11ij, I22ij, I33ij, I12ij, I13ij, I23ij;

  u1 = us[0]; u2 = us[1]; u3 = us[2];
  v1 = vs[0]; v2 = vs[1]; v3 = vs[2];
  w1 = ws[0]; w2 = ws[1]; w3 = ws[2];
  J1 = u2*v3-u3*v2; J2 = u3*v1-u1*v3; J3 = u1*v2-u2*v1;
  J = sqrt(J1*J1+J2*J2+J3*J3);
  J1 = u2*w3-u3*w2; J2 = u3*w1-u1*w3; J3 = u1*w2-u2*w1;
  J *= sqrt(J1*J1+J2*J2+J3*J3);

  I = I11 = I22 = I33 = I12 = I13 = I23 = 0.0;
  for (i=0; i<quadn; i++)
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
      normxy = sqrt(xy1*xy1+xy2*xy2+xy3*xy3);
      normxy3 = normxy*normxy*normxy;
      Iij = 1.0/normxy;
      I11ij = xy1*xy1/normxy3;
      I22ij = xy2*xy2/normxy3;
      I33ij = xy3*xy3/normxy3;
      I12ij = xy1*xy2/normxy3;
      I13ij = xy1*xy3/normxy3;
      I23ij = xy2*xy3/normxy3;
      xy1 = -eta2*(u1+v1+eta3*w1)+v1;
      xy2 = -eta2*(u2+v2+eta3*w2)+v2;
      xy3 = -eta2*(u3+v3+eta3*w3)+v3;
      normxy = sqrt(xy1*xy1+xy2*xy2+xy3*xy3);
      normxy3 = normxy*normxy*normxy;
      Iij += 1.0/normxy;
      I11ij += xy1*xy1/normxy3;
      I22ij += xy2*xy2/normxy3;
      I33ij += xy3*xy3/normxy3;
      I12ij += xy1*xy2/normxy3;
      I13ij += xy1*xy3/normxy3;
      I23ij += xy2*xy3/normxy3;
      xy1 = eta2*(v1-eta3*(u1+v1))-w1;
      xy2 = eta2*(v2-eta3*(u2+v2))-w2;
      xy3 = eta2*(v3-eta3*(u3+v3))-w3;
      normxy = sqrt(xy1*xy1+xy2*xy2+xy3*xy3);
      normxy3 = normxy*normxy*normxy;
      Iij += 1.0/normxy;
      I11ij += xy1*xy1/normxy3;
      I22ij += xy2*xy2/normxy3;
      I33ij += xy3*xy3/normxy3;
      I12ij += xy1*xy2/normxy3;
      I13ij += xy1*xy3/normxy3;
      I23ij += xy2*xy3/normxy3;
      xy1 = -eta2*(w1+eta3*(u1+v1))+v1;
      xy2 = -eta2*(w2+eta3*(u2+v2))+v2;
      xy3 = -eta2*(w3+eta3*(u3+v3))+v3;
      normxy = sqrt(xy1*xy1+xy2*xy2+xy3*xy3);
      normxy3 = normxy*normxy*normxy;
      Iij += 1.0/normxy; Iij *= eta2;
      I11ij += xy1*xy1/normxy3; I11ij *= eta2;
      I22ij += xy2*xy2/normxy3; I22ij *= eta2;
      I33ij += xy3*xy3/normxy3; I33ij *= eta2;
      I12ij += xy1*xy2/normxy3; I12ij *= eta2;
      I13ij += xy1*xy3/normxy3; I13ij *= eta2;
      I23ij += xy2*xy3/normxy3; I23ij *= eta2;
      xy1 = eta2*(u1+w1)+eta3*v1-w1;
      xy2 = eta2*(u2+w2)+eta3*v2-w2;
      xy3 = eta2*(u3+w3)+eta3*v3-w3;
      normxy = sqrt(xy1*xy1+xy2*xy2+xy3*xy3);
      normxy3 = normxy*normxy*normxy;
      Iij += 1.0/normxy;
      I11ij += xy1*xy1/normxy3;
      I22ij += xy2*xy2/normxy3;
      I33ij += xy3*xy3/normxy3;
      I12ij += xy1*xy2/normxy3;
      I13ij += xy1*xy3/normxy3;
      I23ij += xy2*xy3/normxy3;
      I += weta2*weta3*Iij;
      I11 += weta2*weta3*I11ij;
      I22 += weta2*weta3*I22ij;
      I33 += weta2*weta3*I33ij;
      I12 += weta2*weta3*I12ij;
      I13 += weta2*weta3*I13ij;
      I23 += weta2*weta3*I23ij;
    }
  }

  I *= J / M_PI / 24.0;
  I11 *= J / M_PI / 24.0;
  I22 *= J / M_PI / 24.0;
  I33 *= J / M_PI / 24.0;
  I12 *= J / M_PI / 24.0;
  I13 *= J / M_PI / 24.0;
  I23 *= J / M_PI / 24.0;
}


static void SingleLayerElasticity3d_commonVertex (double *us, double *vs,
                         double *ws, double *zs, double &I,
                         double &I11, double &I22,
                         double &I33, double &I12,
                         double &I13, double &I23)

{
  int i, j, k;
  double u1, u2, u3, v1, v2, v3, w1, w2, w3, z1, z2, z3, J1, J2, J3, J;
  double xy1, xy2, xy3, normxy, normxy3, eta1, eta2, eta3, weta1, weta2, weta3;
  double Iijk, I11ijk, I22ijk, I33ijk, I12ijk, I13ijk, I23ijk;

  u1 = us[0]; u2 = us[1]; u3 = us[2];
  v1 = vs[0]; v2 = vs[1]; v3 = vs[2];
  w1 = ws[0]; w2 = ws[1]; w3 = ws[2];
  z1 = zs[0]; z2 = zs[1]; z3 = zs[2];
  J1 = u2*v3-u3*v2; J2 = u3*v1-u1*v3; J3 = u1*v2-u2*v1;
  J = sqrt(J1*J1+J2*J2+J3*J3);
  J1 = w2*z3-w3*z2; J2 = w3*z1-w1*z3; J3 = w1*z2-w2*z1;
  J *= sqrt(J1*J1+J2*J2+J3*J3);

  I = I11 = I22 = I33 = I12 = I13 = I23 = 0.0;
  for (i=0; i<quadn; i++)
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
    normxy3 = normxy*normxy*normxy;
        Iijk = 1.0/normxy;
    I11ijk = xy1*xy1/normxy3;
    I22ijk = xy2*xy2/normxy3;
    I33ijk = xy3*xy3/normxy3;
    I12ijk = xy1*xy2/normxy3;
    I13ijk = xy1*xy3/normxy3;
    I23ijk = xy2*xy3/normxy3;
    xy1 = w1+eta1*z1-eta2*(u1+eta3*v1);
    xy2 = w2+eta1*z2-eta2*(u2+eta3*v2);
    xy3 = w3+eta1*z3-eta2*(u3+eta3*v3);
    normxy = sqrt(xy1*xy1+xy2*xy2+xy3*xy3);
    normxy3 = normxy*normxy*normxy;
    Iijk += 1.0/normxy; Iijk *= eta2;
    I11ijk += xy1*xy1/normxy3; I11ijk *= eta2;
    I22ijk += xy2*xy2/normxy3; I22ijk *= eta2;
    I33ijk += xy3*xy3/normxy3; I33ijk *= eta2;
    I12ijk += xy1*xy2/normxy3; I12ijk *= eta2;
    I13ijk += xy1*xy3/normxy3; I13ijk *= eta2;
    I23ijk += xy2*xy3/normxy3; I23ijk *= eta2;
    I += weta1*weta2*weta3*Iijk;
    I11 += weta1*weta2*weta3*I11ijk;
    I22 += weta1*weta2*weta3*I22ijk;
    I33 += weta1*weta2*weta3*I33ijk;
    I12 += weta1*weta2*weta3*I12ijk;
    I13 += weta1*weta2*weta3*I13ijk;
    I23 += weta1*weta2*weta3*I23ijk;
      }
    }
  }

  I *= J / M_PI / 12.0;
  I11 *= J / M_PI / 12.0;
  I22 *= J / M_PI / 12.0;
  I33 *= J / M_PI / 12.0;
  I12 *= J / M_PI / 12.0;
  I13 *= J / M_PI / 12.0;
  I23 *= J / M_PI / 12.0;
}


static void SingleLayerElasticity3d_disjointPanels
(double *As, double *us, double *vs, double *Bs, double *ws, double *zs,
 double &I, double &I11, double &I22, double &I33, double &I12, double &I13,
 double &I23)
{
  int i, j;
  double A1, A2, A3, u1, u2, u3, v1, v2, v3, B1, B2, B3, w1, w2, w3, z1, z2, z3;
  double J1, J2, J3, J;
  double x1, x2, x3, xy1, xy2, xy3, normxy, normxy3;
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
  J *= sqrt(J1*J1+J2*J2+J3*J3);

  I = I11 = I22 = I33 = I12 = I13 = I23 = 0.0;
  for (i=0; i<quadTrin; i++)
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
      normxy3 = normxy*normxy*normxy;

      I += wksi*weta/normxy;
      I11 += wksi*weta*xy1*xy1/normxy3;
      I22 += wksi*weta*xy2*xy2/normxy3;
      I33 += wksi*weta*xy3*xy3/normxy3;
      I12 += wksi*weta*xy1*xy2/normxy3;
      I13 += wksi*weta*xy1*xy3/normxy3;
      I23 += wksi*weta*xy2*xy3/normxy3;
    }
  }

  I *= J / M_PI / 4.0;
  I11 *= J / M_PI / 4.0;
  I22 *= J / M_PI / 4.0;
  I33 *= J / M_PI / 4.0;
  I12 *= J / M_PI / 4.0;
  I13 *= J / M_PI / 4.0;
  I23 *= J / M_PI / 4.0;
}


static void DoubleLayerLaplace3d_commonEdge(double *us, double *vs, double *ws, double *Kij)
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


static void DoubleLayerLaplace3d_commonVertex(double *us, double *vs, double *ws, double *zs, double *Kij)
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


static void BEM3DElasticity (esint np, double *points, esint ne, esint *elemNodes, int order, double *V, double *K, double *D, double *M, double nu)
{
  esint i, j, k, l, cmnIdxi[2],cmnIdxj[2],cmnIdxSize, remIdxi, remIdxj, tmp;
  esint ne2 = 2*ne, ne3 = 3*ne, np2 = 2*np, np3 = 3*np;
  double factor = (1+nu)/2.0/ /* E/ */ (1-nu), fact = 3-4*nu;
  double factorK = 1.0/2.0/(1-nu), factK = 1-2*nu;
  double V11K, V12K, V13K, V22K, V23K, V33K;
  double mu = 1.0/2.0/(1+nu), mu2 = mu*mu;
  double V11, V12, V13, V22, V23, V33;

  memset(V,0,ne3*ne3*sizeof(double));
  memset(K,0,ne3*np3*sizeof(double));
  memset(D,0,np3*np3*sizeof(double));
  memset(M,0,ne3*np3*sizeof(double));

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
  double Vij, V11ij, V22ij, V33ij, V12ij, V13ij, V23ij, Kij[3];
  double a11, a12, a22, detA, b11, b21, b12, b22, b13, b23;
  double y11, y21, y12, y22, y13, y23;
  double curli1[3], curli2[3], curli3[3], curlj1[3], curlj2[3], curlj3[3];

  for (i=0; i<ne; i++)
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
    {
      M[i+(faceNodesi[k]-1)*ne3] += Ji/6.0;
      M[ne+i+(np+faceNodesi[k]-1)*ne3] += Ji/6.0;
      M[ne2+i+(np2+faceNodesi[k]-1)*ne3] += Ji/6.0;
    }
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

      for (j=0; j<ne; j++)
    {
      for (k=0; k<3; k++)
        faceNodesj[k] = elemNodes[3*j+k] + 1;
      if (j==i) // identical panels
        {
          memcpy(curlj1,curli1,3*sizeof(double));
          memcpy(curlj2,curli2,3*sizeof(double));
          memcpy(curlj3,curli3,3*sizeof(double));
          for (k=0; k<3; k++)
        {
          u[k] = Pi[1][k]-Pi[0][k];
          v[k] = Pi[2][k]-Pi[1][k];
        }
          SingleLayerElasticity3d_idPanels(u,v,Vij,V11ij,V22ij,V33ij,
                           V12ij,V13ij,V23ij);
          Kij[0] = Kij[1] = Kij[2] = 0;
        }

      else
        {
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
          SingleLayerElasticity3d_disjointPanels(Pi[0],u,v,Pj[0],w,z,
                             Vij,V11ij,V22ij,V33ij,
                             V12ij,V13ij,V23ij);
          DoubleLayerLaplace3d_disjointPanels(Pi[0],u,v,Pj[0],w,z,Kij);
          for (k=0; k<3; k++)
            {
              K[i+(faceNodesj[k]-1)*ne3] += Kij[k];
              K[ne+i+(np+faceNodesj[k]-1)*ne3] += Kij[k];
              K[ne2+i+(np2+faceNodesj[k]-1)*ne3] += Kij[k];
            }
          break;
        case 1: // common vertex
          for (k=0; k<3; k++)
            {
              u[k] = Pi[(cmnIdxi[0]+1)%3][k]-Pi[cmnIdxi[0]][k];
              v[k] = Pi[(cmnIdxi[0]+2)%3][k]-Pi[(cmnIdxi[0]+1)%3][k];
              w[k] = Pj[(cmnIdxj[0]+1)%3][k]-Pj[cmnIdxj[0]][k];
              z[k] = Pj[(cmnIdxj[0]+2)%3][k]-Pj[(cmnIdxj[0]+1)%3][k];
            }
          SingleLayerElasticity3d_commonVertex(u,v,w,z,Vij,V11ij,V22ij,
                               V33ij,V12ij,V13ij,V23ij);
          DoubleLayerLaplace3d_commonVertex(u,v,w,z,Kij);
          for (k=0; k<3; k++)
            {
              K[i+(faceNodesj[(cmnIdxj[0]+k)%3]-1)*ne3] += Kij[k];
              K[ne+i+(np+faceNodesj[(cmnIdxj[0]+k)%3]-1)*ne3] += Kij[k];
              K[ne2+i+(np2+faceNodesj[(cmnIdxj[0]+k)%3]-1)*ne3] +=
            Kij[k];
            }
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
          SingleLayerElasticity3d_commonEdge(u,v,w,Vij,V11ij,V22ij,
                             V33ij,V12ij,V13ij,V23ij);
          DoubleLayerLaplace3d_commonEdge(u,v,w,Kij);
          K[i+(faceNodesj[cmnIdxj[0]]-1)*ne3] += Kij[0];
          K[ne+i+(np+faceNodesj[cmnIdxj[0]]-1)*ne3] += Kij[0];
          K[ne2+i+(np2+faceNodesj[cmnIdxj[0]]-1)*ne3] += Kij[0];
          K[i+(faceNodesj[cmnIdxj[1]]-1)*ne3] += Kij[1];
          K[ne+i+(np+faceNodesj[cmnIdxj[1]]-1)*ne3] += Kij[1];
          K[ne2+i+(np2+faceNodesj[cmnIdxj[1]]-1)*ne3] += Kij[1];
          K[i+(faceNodesj[remIdxj]-1)*ne3] += Kij[2];
          K[ne+i+(np+faceNodesj[remIdxj]-1)*ne3] += Kij[2];
          K[ne2+i+(np2+faceNodesj[remIdxj]-1)*ne3] += Kij[2];
          break;
        }
        }

      // ***
      //  V
      // ***
      V[i+j*ne3] = V11 = factor * (fact*Vij + V11ij);
      V[i+(ne+j)*ne3] = V12 = factor * V12ij;
      V[i+(ne2+j)*ne3] = V13 = factor * V13ij;
      V[ne+i+j*ne3] = factor * V12ij;
      V[ne+i+(ne+j)*ne3] = V22 = factor * (fact*Vij + V22ij);
      V[ne+i+(ne2+j)*ne3] = V23 = factor * V23ij;
      V[ne2+i+j*ne3] = factor * V13ij;
      V[ne2+i+(ne+j)*ne3] = factor * V23ij;
      V[ne2+i+(ne2+j)*ne3] = V33 = factor * (fact*Vij + V33ij);

      // ***
      //  K
      // ***
      V11K = factorK * (factK*Vij + V11ij);
      V12K = factorK * V12ij;
      V13K = factorK * V13ij;
      V22K = factorK * (factK*Vij + V22ij);
      V23K = factorK * V23ij;
      V33K = factorK * (factK*Vij + V33ij);
      for (k=0; k<3; k++)
        { // K11, K12, K13, K21, K22, K23, K31, K32, K33 ... change signs!
          K[i+(faceNodesj[k]-1)*ne3] +=
        V12K*curlj3[k]-V13K*curlj2[k];
          K[i+(np+faceNodesj[k]-1)*ne3] +=
        -V11K*curlj3[k]+V13K*curlj1[k];
          K[i+(np2+faceNodesj[k]-1)*ne3] +=
        V11K*curlj2[k]-V12K*curlj1[k];
          K[ne+i+(faceNodesj[k]-1)*ne3] +=
        V22K*curlj3[k]-V23K*curlj2[k];
          K[ne+i+(np+faceNodesj[k]-1)*ne3] +=
        -V12K*curlj3[k]+V23K*curlj1[k];
          K[ne+i+(np2+faceNodesj[k]-1)*ne3] +=
        V12K*curlj2[k]-V22K*curlj1[k];
          K[ne2+i+(faceNodesj[k]-1)*ne3] +=
        V23K*curlj3[k]-V33K*curlj2[k];
          K[ne2+i+(np+faceNodesj[k]-1)*ne3] +=
        -V13K*curlj3[k]+V33K*curlj1[k];
          K[ne2+i+(np2+faceNodesj[k]-1)*ne3] +=
        V13K*curlj2[k]-V23K*curlj1[k];
        }

      // ***
      //  D
      // ***
      //
      for (k=0; k<3; k++)
        for (l=0; l<3; l++)
          {
        // D11
        D[faceNodesi[k]-1+(faceNodesj[l]-1)*np3] +=
          mu * (4*(curli3[k]*Vij*curlj3[l]+curli2[k]*Vij*curlj2[l])
            +curli1[k]*Vij*curlj1[l]);
        D[faceNodesi[k]-1+(faceNodesj[l]-1)*np3] -=
          4*mu2 * (curli3[k]*V22*curlj3[l]+curli2[k]*V33*curlj2[l]
               -curli3[k]*V23*curlj2[l]-curli2[k]*V23*curlj3[l]);
        // D12
        D[faceNodesi[k]-1+(np+faceNodesj[l]-1)*np3] +=
          mu * (-2*curli2[k]*Vij*curlj1[l]-curli1[k]*Vij*curlj2[l]);
        D[faceNodesi[k]-1+(np+faceNodesj[l]-1)*np3] -=
          4*mu2 * (-curli3[k]*V12*curlj3[l]+curli2[k]*V13*curlj3[l]
               +curli3[k]*V23*curlj1[l]-curli2[k]*V33*curlj1[l]);
        // D13
        D[faceNodesi[k]-1+(np2+faceNodesj[l]-1)*np3] -=
          mu * (2*curli3[k]*Vij*curlj1[l]+curli1[k]*Vij*curlj3[l]);
        D[faceNodesi[k]-1+(np2+faceNodesj[l]-1)*np3] -=
          4*mu2 * (-curli2[k]*V13*curlj2[l]+curli3[k]*V12*curlj2[l]
               -curli3[k]*V22*curlj1[l]+curli2[k]*V23*curlj1[l]);
        // D21
        D[np+faceNodesi[k]-1+(faceNodesj[l]-1)*np3] +=
          mu * (-2*curli1[k]*Vij*curlj2[l]-curli2[k]*Vij*curlj1[l]);
        D[np+faceNodesi[k]-1+(faceNodesj[l]-1)*np3] -=
          4*mu2 * (-curli3[k]*V12*curlj3[l]+curli3[k]*V13*curlj2[l]
               +curli1[k]*V23*curlj3[l]-curli1[k]*V33*curlj2[l]);
        // D22
        D[np+faceNodesi[k]-1+(np+faceNodesj[l]-1)*np3] +=
          mu * (4*(curli3[k]*Vij*curlj3[l]+curli1[k]*Vij*curlj1[l])
            +curli2[k]*Vij*curlj2[l]);
        D[np+faceNodesi[k]-1+(np+faceNodesj[l]-1)*np3] -=
          4*mu2 * (curli3[k]*V11*curlj3[l]+curli1[k]*V33*curlj1[l]
               -curli3[k]*V13*curlj1[l]-curli1[k]*V13*curlj3[l]);
        // D23
        D[np+faceNodesi[k]-1+(np2+faceNodesj[l]-1)*np3] +=
          mu * (-2*curli3[k]*Vij*curlj2[l]-curli2[k]*Vij*curlj3[l]);
        D[np+faceNodesi[k]-1+(np2+faceNodesj[l]-1)*np3] -=
          4*mu2 * (-curli1[k]*V23*curlj1[l]+curli1[k]*V13*curlj2[l]
               -curli3[k]*V11*curlj2[l]+curli3[k]*V12*curlj1[l]);
        // D31
        D[np2+faceNodesi[k]-1+(faceNodesj[l]-1)*np3] -=
          mu * (2*curli1[k]*Vij*curlj3[l]+curli3[k]*Vij*curlj1[l]);
        D[np2+faceNodesi[k]-1+(faceNodesj[l]-1)*np3] -=
          4*mu2 * (-curli2[k]*V13*curlj2[l]+curli2[k]*V12*curlj3[l]
               -curli1[k]*V22*curlj3[l]+curli1[k]*V23*curlj2[l]);
        // D32
        D[np2+faceNodesi[k]-1+(np+faceNodesj[l]-1)*np3] +=
          mu * (-2*curli2[k]*Vij*curlj3[l]-curli3[k]*Vij*curlj2[l]);
        D[np2+faceNodesi[k]-1+(np+faceNodesj[l]-1)*np3] -=
          4*mu2 * (-curli1[k]*V23*curlj1[l]+curli2[k]*V13*curlj1[l]
               -curli2[k]*V11*curlj3[l]+curli1[k]*V12*curlj3[l]);
        // D33
        D[np2+faceNodesi[k]-1+(np2+faceNodesj[l]-1)*np3] +=
          mu * (4*(curli2[k]*Vij*curlj2[l]+curli1[k]*Vij*curlj1[l])
            +curli3[k]*Vij*curlj3[l]);
        D[np2+faceNodesi[k]-1+(np2+faceNodesj[l]-1)*np3] -=
          4*mu2 * (curli2[k]*V11*curlj2[l]+curli1[k]*V22*curlj1[l]
               -curli2[k]*V12*curlj1[l]-curli1[k]*V12*curlj2[l]);
          }
    }
    }
}

void BEM3dElasticityVolume (esint nn, double *nodes,
                esint ne, esint *elemNodes,
                esint np, double *points, int order,
                //E = 1
                double nu,
                double *t, double *u, double *Vt_Wu)
{
  double lambda = /* E* */ nu/(1.0+nu)/(1.0-2.0*nu), mu = /* E* */ 0.5/(1+nu);
  esint i, j, k, faceNodesj[3];
  double ksi1, ksi2, wksi, A1, A2, A3, u1, u2, u3, v1, v2, v3, n1, n2, n3,Jacob;
  double V1, V2, V3, W1, W2, W3;
  double x1, x2, x3, xy1, xy2, xy3, xy12, xy13, xy22, xy23, xy32, xy33;
  double normxy, normxy3, normxy5;
  double t1j, t2j, t3j, u1j1, u1j2, u1j3, u2j1, u2j2, u2j3, u3j1, u3j2, u3j3;
  double A, B, C, D, E, F, G, H, I, J, K, L, M;
  double U11, U12, U13, U22, U23, U33;
  double U111, U112, U113, U121, U122, U123, U131, U132, U133;
  double U221, U222, U223, U231, U232, U233, U331, U332, U333;
  double u1ksi, u2ksi, u3ksi, un11, un22, un33, un, un12, un13, un23;

  quadTrin = TriangleQuadratureSizes[order];
  qTriPoints = TriangleQuadraturePoints[order];
  qTriWeights = TriangleQuadratureWeights[order];

  memset(Vt_Wu,0,3*np*sizeof(double));
  for (j=0; j<ne; j++)
    {
      t1j = t[j];
      t2j = t[ne+j];
      t3j = t[2*ne+j];
      for (k=0; k<3; k++)
    faceNodesj[k] = elemNodes[3*j+k] + 1;
      u1j1 = u[faceNodesj[0]-1];
      u1j2 = u[faceNodesj[1]-1];
      u1j3 = u[faceNodesj[2]-1];
      u2j1 = u[nn+faceNodesj[0]-1];
      u2j2 = u[nn+faceNodesj[1]-1];
      u2j3 = u[nn+faceNodesj[2]-1];
      u3j1 = u[2*nn+faceNodesj[0]-1];
      u3j2 = u[2*nn+faceNodesj[1]-1];
      u3j3 = u[2*nn+faceNodesj[2]-1];
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
      Jacob = sqrt(n1*n1+n2*n2+n3*n3);
      n1 /= Jacob; n2 /= Jacob; n3 /= Jacob;
      for (i=0; i<np; i++)
    {
      x1 = points[3*i];
      x2 = points[3*i+1];
      x3 = points[3*i+2];
      for (k=0, V1=0.0, V2=0.0, V3=0.0, W1=0.0, W2=0.0, W3=0.0;
           k<quadTrin; k++)
        {
          wksi = qTriWeights[k];
          ksi1 = qTriPoints[k].x;
          ksi2 = qTriPoints[k].y;
          u1ksi = u1j1+(u1j2-u1j1)*ksi1+(u1j3-u1j1)*ksi2;
          u2ksi = u2j1+(u2j2-u2j1)*ksi1+(u2j3-u2j1)*ksi2;
          u3ksi = u3j1+(u3j2-u3j1)*ksi1+(u3j3-u3j1)*ksi2;
          un11 = u1ksi*n1; un22 = u2ksi*n2; un33 = u3ksi*n3;
          un = un11+un22+un33;
          un12 = u1ksi*n2 + u2ksi*n1;
          un13 = u1ksi*n3 + u3ksi*n1;
          un23 = u2ksi*n3 + u3ksi*n2;
          xy1 = x1 - A1 - u1*ksi1 - v1*ksi2;
          xy12 = xy1*xy1; xy13 = xy12*xy1;
          xy2 = x2 - A2 - u2*ksi1 - v2*ksi2;
          xy22 = xy2*xy2; xy23 = xy22*xy2;
          xy3 = x3 - A3 - u3*ksi1 - v3*ksi2;
          xy32 = xy3*xy3; xy33 = xy32*xy3;
          normxy = sqrt(xy1*xy1+xy2*xy2+xy3*xy3);
          normxy3 = normxy*normxy*normxy;
          normxy5 = normxy3*normxy*normxy;
          U11 = (3.0-4.0*nu)/normxy + xy12/normxy3;
          U22 = (3.0-4.0*nu)/normxy + xy22/normxy3;
          U33 = (3.0-4.0*nu)/normxy + xy32/normxy3;
          U12 = xy1*xy2/normxy3; U13 = xy1*xy3/normxy3; U23=xy2*xy3/normxy3;
          A = xy1/normxy3; B = xy2/normxy3; C = xy3/normxy3;
          D = 3.0*xy13/normxy5; E = 3.0*xy23/normxy5; F = 3.0*xy33/normxy5;
          G = 3.0*xy12*xy2/normxy5; H = 3.0*xy12*xy3/normxy5;
          I = 3.0*xy22*xy1/normxy5; J = 3.0*xy22*xy3/normxy5;
          K = 3.0*xy32*xy1/normxy5; L = 3.0*xy32*xy2/normxy5;
          M = 3.0*xy1*xy2*xy3/normxy5;
          U111 = (1.0-4.0*nu)*A + D;
          U222 = (1.0-4.0*nu)*B + E;
          U333 = (1.0-4.0*nu)*C + F;
          U112 = (3.0-4.0*nu)*B + G;
          U113 = (3.0-4.0*nu)*C + H;
          U221 = (3.0-4.0*nu)*A + I;
          U223 = (3.0-4.0*nu)*C + J;
          U331 = (3.0-4.0*nu)*A + K;
          U332 = (3.0-4.0*nu)*B + L;
          U121 = G-B; U131 = H-C; U232 = J-C;
          U122 = I-A; U133 = K-A; U233 = L-B;
          U123 = U132 = U231 = M;
          V1 += wksi * (U11*t1j+U12*t2j+U13*t3j);
          V2 += wksi * (U12*t1j+U22*t2j+U23*t3j);
          V3 += wksi * (U13*t1j+U23*t2j+U33*t3j);
          W1 += wksi * ( lambda*(U111+U122+U133)*un +
                 2*mu*(U111*un11+U122*un22+U133*un33) +
                 mu*((U112+U121)*un12+(U113+U131)*un13+
                 (U123+U132)*un23) );
          W2 += wksi * ( lambda*(U121+U222+U233)*un +
                 2*mu*(U121*un11+U222*un22+U233*un33) +
                 mu*((U122+U221)*un12+(U123+U231)*un13+
                 (U223+U232)*un23) );
          W3 += wksi * ( lambda*(U131+U232+U333)*un +
                 2*mu*(U131*un11+U232*un22+U333*un33) +
                 mu*((U132+U231)*un12+(U133+U331)*un13+
                 (U233+U332)*un23) );
        }
      Vt_Wu[i] += Jacob * (V1-W1) * (1.0+nu) / 8.0 / M_PI / (1.0-nu);
      Vt_Wu[np+i] += Jacob * (V2-W2) * (1.0+nu) / 8.0 / M_PI / (1.0-nu);
      Vt_Wu[2*np+i] += Jacob * (V3-W3) * (1.0+nu) / 8.0 / M_PI / (1.0-nu);
    }
    }
}



//static void invertMatrix (esint n, double *M, double *invM)
//{
//  esint i, j, k;
//  double pivot, factor, entry;
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

//
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


//void deleteMyBEMData (MyBEMData* &bem)
//{
//  bem->nNodes = 0;
//  if (bem->S)
//    {
//      delete [] bem->S; bem->S = 0;
//    }
//}


void BEM3DElasticity(double *K, esint np, double *points, esint ne, esint *elements, double YoungModulus, double PoissonRatio)
{
  int i, j, idx, nElements3 = 3 * ne, nNodes3 = 3 * np;
  double *V, *KK, *D, *M;
  int order = 7;
  /*
  deleteMyBEMData(bem);
  bem->nNodes = nNodes;
  */
  V = new double[nElements3*nElements3];
  KK = new double[nElements3*nNodes3];
  D = new double[nNodes3*nNodes3];
  M = new double[nElements3*nNodes3];
  BEM3DElasticity(np, points, ne, elements, order, V, KK, D, M, PoissonRatio);
  /*
  mexPrintfMatrix(nElements3,nElements3,V,"V");
  mexPrintfMatrix(nElements3,nNodes3,KK,"K");
  mexPrintfMatrix(nElements3,nNodes3,M,"M");
  mexPrintfMatrix(nNodes3,nNodes3,D,"D");
  */
  //invV = new double[nElements3*nElements3];
  //invertMatrix(nElements3,V,invV);
  for (i=0, idx=0; i<nElements3; i++)
    for (j=0; j<nNodes3; j++, idx++)
      KK[idx] += 0.5*M[idx];
  //K = new double[nNodes3*nNodes3];
  //assembleSteklovPoincare(nElements3,nNodes3,invV,KK,D,K);
  ///*
  cholesky(nElements3,V);
  assembleSteklovPoincareCholesky(nElements3,nNodes3,V,KK,D,K);
  //*/
  for (i=0, idx=0; i<nNodes3; i++)
    for (j=0; j<nNodes3; j++, idx++)
      K[idx] *= YoungModulus;
  //bem->S = K;
  delete [] V;
  //delete [] invV;
  delete [] KK;
  delete [] D;
  delete [] M;
}

void BEM3DElasticityEval(double *results, esint np, double *points, esint ne, esint *elements, esint ni, double *inner, double YoungModulus, double PoissonRatio, double *dirichlet)
{
    esint i, j, idx, nElements3 = 3*ne, nNodes3 = 3*np;
    double *rhs, *neumann;
    double *V, *KK, *D, *M;
    int order = 7;
    /*
    deleteMyBEMData(bem);
    bem->nNodes = nNodes;
    */
    V = new double[nElements3*nElements3];
    KK = new double[nElements3*nNodes3];
    D = new double[nNodes3*nNodes3];
    M = new double[nElements3*nNodes3];
    BEM3DElasticity(np,points,ne,elements,order,V,KK,D,M,PoissonRatio);
    for (i=0, idx=0; i<nElements3; i++)
      for (j=0; j<nNodes3; j++, idx++)
        KK[idx] += 0.5*M[idx];

    rhs = new double[nElements3];
    memset(rhs,0,nElements3*sizeof(double));
    for (j=0, idx=0; j<nNodes3; j++)
      for (i=0; i<nElements3; i++, idx++)
        rhs[i] += KK[idx]*dirichlet[j];
    neumann = new double[nElements3];
    cholesky(nElements3,V);
    choleskySolve(nElements3,V,rhs,neumann);

    BEM3dElasticityVolume(np,points,ne,elements,ni,inner,order, PoissonRatio,neumann,dirichlet,results);
    delete [] neumann;
    delete [] rhs;
    delete [] V;
    delete [] KK;
    delete [] D;
    delete [] M;
}

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
//  double YoungModulus, PoissonRatio;
//
//  n = mxGetN(prhs[0]);
//  nodes = mxGetPr(prhs[0]);
//
//  m = mxGetN(prhs[1]);
//  elemNodes = mxGetPr(prhs[1]);
//
//  order = mxGetScalar(prhs[2]);
//
//  YoungModulus = mxGetScalar(prhs[3]);
//  PoissonRatio = mxGetScalar(prhs[4]);
//
//  /*
//  plhs[0] = mxCreateDoubleMatrix(3*m,3*m,mxREAL);
//  V = mxGetPr(plhs[0]);
//  plhs[1] = mxCreateDoubleMatrix(3*m,3*n,mxREAL);
//  K = mxGetPr(plhs[1]);
//  plhs[2] = mxCreateDoubleMatrix(3*n,3*n,mxREAL);
//  D = mxGetPr(plhs[2]);
//  plhs[3] = mxCreateDoubleMatrix(3*m,3*n,mxREAL);
//  M = mxGetPr(plhs[3]);
//  BEM3dElasticity(n,nodes,m,elemNodes,order,V,K,D,M,PoissonRatio);
//  */
//
//  ///*
//  plhs[0] = mxCreateDoubleMatrix(3*n,3*n,mxREAL);
//  K = mxGetPr(plhs[0]);
//  getElasticity(K,n,nodes,m,elemNodes,order,YoungModulus,PoissonRatio);
//  //*/
//
//  return;
//}
