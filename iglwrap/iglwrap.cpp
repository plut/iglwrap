// TODO: write conversion functions Julia ↔ C++
#include <Eigen/Dense>
#include <stdlib.h>
#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include <igl/copyleft/cgal/mesh_boolean.h>
#include <igl/copyleft/cgal/minkowski_sum.h>
#include <igl/copyleft/cgal/piecewise_constant_winding_number.h>
#include <igl/offset_surface.h>
#include <igl/decimate.h>
#include <igl/loop.h>
#include <igl/copyleft/cgal/intersect_with_half_space.h>
#include "iglwrap.h"
#include <boost/bind/bind.hpp>

using namespace boost::placeholders;
using namespace Eigen;

typedef Matrix<double,Dynamic,3> VertexMatrix;
typedef Matrix<int,Dynamic,3> FaceMatrix;

// typedef Matrix<double, Dynamic, 3, RowMajor> vertices_t;
// typedef Matrix<int, Dynamic, 3, RowMajor> faces_t;

// Eigen::IOFormat CommaInitFmt(32, DontAlignCols, ", ", ", ", "", "", " << ", ";");

struct Mesh {/*««*/
	MatrixXd v;
	MatrixXi f;
//   Matrix<double,Dynamic,3> v;
//   Matrix<int,Dynamic,3> f;
  Mesh(): v(), f() { }
  Mesh(int nv, int nf, const double *mv, const int *mf, int d = 3):
			v(nv, d), f(nf, d) {
    for(int i = 0; i < nv; i++) {
      for(int j = 0; j < d; j++) {
				v(i, j) = mv[d*i+j];
      }
    }
    for(int i = 0; i < nf; i++) {
      for(int j = 0; j < d; j++) {
	f(i, j) = mf[d*i+j]-1;
      }
    };
  };
};/*»»*/

void to_jl(Mesh m, int *nv, int *nf, double **mv, int **mf) {/*««*/
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
}/*»»*/

int mesh_boolean(/*««*/
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
}/*»»*/
int minkowski_sum(/*««*/
  int nv1, int nf1, const double *mv1, const int *mf1,
  int nv2, int nf2, const double *mv2, const int *mf2, int dim2,
  int *nv3, int *nf3, double **mv3, int **mf3, int **index) {

  Mesh m1(nv1, nf1, mv1, mf1),
       m2(nv2, nf2, mv2, mf2, dim2),
       m3;
  MatrixXi j;
// 	std::cout << "mesh m1: "
// 		<< m1.v.rows() << " vertices:\n" << m1.v << "\n"
// 		<< m1.f.rows() << " faces:\n" << m1.f << "\n\n";
// 	std::cout << "mesh m2: "
// 		<< m2.v.rows() << " vertices:\n" << m2.v << "\n"
// 		<< m2.f.rows() << " faces:\n" << m2.f << "\n\n";
  igl::copyleft::cgal::minkowski_sum(m1.v,m1.f,m2.v,m2.f, true,
  	m3.v,m3.f, j);
// 	std::cout << "mesh m3: "
// 		<< m3.v.rows() << " vertices:\n" << m3.v << "\n"
// 		<< m3.f.rows() << " faces:\n" << m3.f << "\n\n";
// 	std::cout << "j has " << j.rows() << " rows:\n"  << j << "\n";

	int n3 = m3.f.rows();
  *index = (int *) malloc(2*n3*sizeof(int));
  if (*index == 0) {
    return -1;
  }
  for(int i = 0; i < n3; i++) {
    (*index)[i] = j(i,0)+1;
		(*index)[i+n3] = j(i,1)+1;
  }
  to_jl(m3, nv3, nf3, mv3, mf3);
  if(*mv3 == 0) {
    free(*index);
    return -1;
  }
  return 0;
}/*»»*/

int mesh_is_pwn(int nv, int nf, const double *mv, const int *mf) {
  Mesh m(nv, nf, mv, mf);
  return (int)igl::copyleft::cgal::piecewise_constant_winding_number(m.v, m.f);
}

int offset_surface(/*««*/
	int nv, int nf, const double *mv, const int *mf,
	double level, int grid,
	int *nvout, int *nfout, double **mvout, int **mfout) {

	Mesh m(nv, nf, mv, mf), mout;
	MatrixXd GV;
	MatrixXi side;
	MatrixXd S;

	igl::offset_surface(m.v, m.f, level, grid,
		igl::SignedDistanceType::SIGNED_DISTANCE_TYPE_DEFAULT,
		mout.v, mout.f, GV, side, S);

	to_jl(mout, nvout, nfout, mvout, mfout);
	if(*mvout == 0) {
		return -1;
	}
	return 0;
}/*»»*/
int decimate(int nv, int nf, const double *mv, const int *mf, int faces,/*««*/
	int *nvout, int *nfout, double **mvout, int **mfout, int **index) {
	Mesh m(nv, nf, mv, mf), mout;
	VectorXi fidx, vidx;
	igl::decimate(m.v, m.f, faces,
		mout.v, mout.f, fidx, vidx);

  *index = (int *) malloc(fidx.rows()*sizeof(int));
  if (*index == 0) {
    return -1;
  }
  for(int i = 0; i < fidx.rows(); i++) {
    (*index)[i] = fidx(i)+1;
  }
	to_jl(mout, nvout, nfout, mvout, mfout);
	if(*mvout == 0) {
		free(*index);
		return -1;
	}
	return 0;
}/*»»*/
int loop(int nv, int nf, const double *mv, const int *mf, int n,/*««*/
	int *nvout, int *nfout, double **mvout, int **mfout) {
	Mesh m(nv, nf, mv, mf), mout;
	VectorXi fidx, vidx;

	igl::loop(m.v, m.f, mout.v, mout.f, n);

	to_jl(mout, nvout, nfout, mvout, mfout);
	return 0;
}/*»»*/

int intersect_with_half_space(int nv, int nf, /*««*/
	const double *mv, const int *mf,
	const struct vec3d *p, const struct vec3d *n,
	int *nvout, int *nfout, double **mvout, int **mfout, int **index) {
	Mesh m(nv, nf, mv, mf), mout;
	Vector3d vp; vp << p->x, p->y, p->z;
	Vector3d vn; vn << n->x, n->y, n->z;
	VectorXi fidx;

	igl::copyleft::cgal::intersect_with_half_space(m.v, m.f, vp, vn,
		mout.v, mout.f, fidx);

  *index = (int *) malloc(fidx.rows()*sizeof(int));
  if (*index == 0) {
    return -1;
  }
  for(int i = 0; i < fidx.rows(); i++) {
    (*index)[i] = fidx(i)+1;
  }
	to_jl(mout, nvout, nfout, mvout, mfout);
	if(*mvout == 0) {
		free(*index);
		return -1;
	}
	return 0;
}/*»»*/

// vim: noet ts=2 sw=2 fmr=««,»»:
