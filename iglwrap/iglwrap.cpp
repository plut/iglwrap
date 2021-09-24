// TODO: write conversion functions Julia ↔ C++
#include <Eigen/Dense>
#include <stdlib.h>
#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include <igl/copyleft/cgal/mesh_boolean.h>
#include <igl/copyleft/cgal/piecewise_constant_winding_number.h>
#include "iglwrap.h"
#include <boost/bind/bind.hpp>
using namespace boost::placeholders;
using namespace Eigen;

typedef Matrix<double,Dynamic,3> VertexMatrix;
typedef Matrix<int,Dynamic,3> FaceMatrix;

// typedef Matrix<double, Dynamic, 3, RowMajor> vertices_t;
// typedef Matrix<int, Dynamic, 3, RowMajor> faces_t;

Eigen::IOFormat CommaInitFmt(32, DontAlignCols, ", ", ", ", "", "", " << ", ";");

struct Mesh {/*<<<*/
  Matrix<double,Dynamic,3> v;
  Matrix<int,Dynamic,3> f;
  Mesh(): v(), f() { }
  Mesh(int nv, int nf, const double *mv, const int *mf): v(nv, 3), f(nf, 3) {
    for(int i = 0; i < nv; i++) {
      for(int j = 0; j < 3; j++) {
	v(i, j) = mv[3*i+j];
      }
    }
    for(int i = 0; i < nf; i++) {
      for(int j = 0; j < 3; j++) {
	f(i, j) = mf[3*i+j]-1;
      }
    };
  };
};/*>>>*/

void to_jl(Mesh m, int *nv, int *nf, double **mv, int **mf) {/*<<<*/
  // sets mv to zero in case of failure
  *nv = m.v.rows();
  *nf = m.f.rows();
  *mv = (double *) malloc(3*(*nv)*sizeof(double));
  if(*mv == 0) {
    return;
  }
  *mf = (int *) malloc(3* (*nf)*sizeof(int));
  if(*mf == 0) {
    free(mv);
    *mv = 0;
    return;
  }
  for(int i = 0; i < *nv; i++) {
    for(int j = 0; j < 3; j++) {
      (*mv)[3*i+j] = m.v(i, j);
    }
  }
  for(int i = 0; i < *nf; i++) {
    for(int j = 0; j < 3; j++) {
      (*mf)[3*i+j] = 1 + m.f(i,j);
    }
  }
}/*>>>*/

int igl_mesh_boolean(
  int op,
  int nv1, int nf1, const double *mv1, const int *mf1,
  int nv2, int nf2, const double *mv2, const int *mf2,
  int *nv3, int *nf3, double **mv3, int **mf3, int **index) {


  Mesh m1(nv1, nf1, mv1, mf1),
       m2(nv2, nf2, mv2, mf2),
       m3;
  VectorXi j;
  igl::copyleft::cgal::mesh_boolean(m1.v,m1.f,m2.v,m2.f,
  igl::MeshBooleanType(op), m3.v,m3.f, j);

  *index = (int *) malloc(m3.f.rows()*sizeof(int));
  if (*index == 0) {
    return -1;
  }
  for(int i = 0; i < m3.f.rows(); i++) {
    (*index)[i] = j(i)+1;
  }
  to_jl(m3, nv3, nf3, mv3, mf3);
  if(*mv3 == 0) {
    free(*index);
    return -1;
  }
  return 0;
}

int igl_mesh_is_pwn(int nv, int nf, const double *mv, const int *mf) {
  Mesh m(nv, nf, mv, mf);
  return (int)igl::copyleft::cgal::piecewise_constant_winding_number(m.v, m.f);
}
