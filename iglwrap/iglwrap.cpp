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
    std::cout << "\e[34mreading array of size " << (3*nv) << ":\e[m\n";
    for(int i = 0; i < 3*nv; i++) {
      std::cout << "  [" << i << "] = " << mv[i] << "\n";
    }
    std::cout << "(nv,nf)" << nv << " " << nf << "\n";
    for(int i = 0; i < nv; i++) {
      for(int j = 0; j < 3; j++) {
	std::cout << "v("<<i<<","<<j<<")=mv["<<(3*i+j)<<"]=";
	std::cout << mv[3*i+j] << "\n";
	v(i, j) = mv[3*i+j];
      }
    }
    for(int i = 0; i < nf; i++) {
      for(int j = 0; j < 3; j++) {
	f(i, j) = mf[3*i+j]-1;
	std::cout << "f("<<i<<","<<j<<")=mf["<<(3*i+j)<<"]=";
	std::cout << mf[3*i+j] << "\n";
      }
    };
    std::cout << "after init:\n";
    std::cout << "v=" << v << "\n";
    std::cout << "f=" << f << "\n";
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
}/*>>>*/

#define test_values \
std::cout << "before test_values\n"; \
Matrix<double, Dynamic, 3> tv1(26, 3), tv2(8,3); \
Matrix<int, Dynamic, 3> tf1(48, 3), tf2(12,3); \
tv1 << 3, 0, 0, 2.656368076959629842548338274355, 1.3941695161313054640572772768792, 0, 1.7041942401934677686625718706637, 2.4689515976809692787696803861763, 0, 0.3616100407659696425177742185042, 2.9781266222941620291919662122382, 0, -1.0638146611276062536433073546505, 2.8050487280562448333398606337141, 0, -2.2455322445133027997599128866568, 1.9893679747223862452898401897983, 0, -2.9128254522781560353905661031604, 0.71794699286267404936268121673493, 0, -2.9128254522781560353905661031604, -0.71794699286267404936268121673493, 0, -2.2455322445133027997599128866568, -1.9893679747223862452898401897983, 0, -1.0638146611276062536433073546505, -2.8050487280562448333398606337141, 0, 0.3616100407659696425177742185042, -2.9781266222941620291919662122382, 0, 1.7041942401934677686625718706637, -2.4689515976809692787696803861763, 0, 2.656368076959629842548338274355, -1.3941695161313054640572772768792, 0, 3, 0, 10, 2.656368076959629842548338274355, 1.3941695161313054640572772768792, 10, 1.7041942401934677686625718706637, 2.4689515976809692787696803861763, 10, 0.3616100407659696425177742185042, 2.9781266222941620291919662122382, 10, -1.0638146611276062536433073546505, 2.8050487280562448333398606337141, 10, -2.2455322445133027997599128866568, 1.9893679747223862452898401897983, 10, -2.9128254522781560353905661031604, 0.71794699286267404936268121673493, 10, -2.9128254522781560353905661031604, -0.71794699286267404936268121673493, 10, -2.2455322445133027997599128866568, -1.9893679747223862452898401897983, 10, -1.0638146611276062536433073546505, -2.8050487280562448333398606337141, 10, 0.3616100407659696425177742185042, -2.9781266222941620291919662122382, 10, 1.7041942401934677686625718706637, -2.4689515976809692787696803861763, 10, 2.656368076959629842548338274355, -1.3941695161313054640572772768792, 10; \
tv2 << 0, 0, 0, 0, 0, 3, 0, 3, 0, 0, 3, 3, 3, 0, 0, 3, 0, 3, 3, 3, 0, 3, 3, 3; \
tf1 << 8, 7, 5, 5, 0, 8, 6, 5, 7, 3, 5, 4, 10, 9, 8, 12, 11, 10, 0, 12, 10, 1, 3, 2, 3, 1, 0, 3, 0, 5, 10, 8, 0, 18, 20, 21, 21, 13, 18, 20, 18, 19, 17, 18, 16, 21, 22, 23, 23, 24, 25, 23, 25, 13, 15, 16, 14, 13, 14, 16, 18, 13, 16, 13, 21, 23, 0, 1, 13, 1, 2, 14, 2, 3, 15, 3, 4, 16, 4, 5, 17, 5, 6, 18, 6, 7, 19, 7, 8, 20, 8, 9, 21, 9, 10, 22, 10, 11, 23, 11, 12, 24, 12, 0, 25, 1, 14, 13, 2, 15, 14, 3, 16, 15, 4, 17, 16, 5, 18, 17, 6, 19, 18, 7, 20, 19, 8, 21, 20, 9, 22, 21, 10, 23, 22, 11, 24, 23, 12, 25, 24, 0, 13, 25; \
tf2 << 5, 4, 6, 6, 7, 5, 6, 2, 3, 3, 7, 6, 3, 1, 5, 5, 7, 3, 4, 0, 2, 2, 6, 4, 1, 0, 4, 4, 5, 1, 2, 0, 1, 1, 3, 2; \
std::cout << "after test_values\n"

    void test(void) {/*<<<*/
      // igl_mesh_boolean:
    Matrix<double, Dynamic, 3> tv1(26, 3), tv2(8,3);
    Matrix<int, Dynamic, 3> tf1(48, 3), tf2(12,3);
    tv1 << 3, 0, 0, 2.656368076959629842548338274355, 1.3941695161313054640572772768792, 0, 1.7041942401934677686625718706637, 2.4689515976809692787696803861763, 0, 0.3616100407659696425177742185042, 2.9781266222941620291919662122382, 0, -1.0638146611276062536433073546505, 2.8050487280562448333398606337141, 0, -2.2455322445133027997599128866568, 1.9893679747223862452898401897983, 0, -2.9128254522781560353905661031604, 0.71794699286267404936268121673493, 0, -2.9128254522781560353905661031604, -0.71794699286267404936268121673493, 0, -2.2455322445133027997599128866568, -1.9893679747223862452898401897983, 0, -1.0638146611276062536433073546505, -2.8050487280562448333398606337141, 0, 0.3616100407659696425177742185042, -2.9781266222941620291919662122382, 0, 1.7041942401934677686625718706637, -2.4689515976809692787696803861763, 0, 2.656368076959629842548338274355, -1.3941695161313054640572772768792, 0, 3, 0, 10, 2.656368076959629842548338274355, 1.3941695161313054640572772768792, 10, 1.7041942401934677686625718706637, 2.4689515976809692787696803861763, 10, 0.3616100407659696425177742185042, 2.9781266222941620291919662122382, 10, -1.0638146611276062536433073546505, 2.8050487280562448333398606337141, 10, -2.2455322445133027997599128866568, 1.9893679747223862452898401897983, 10, -2.9128254522781560353905661031604, 0.71794699286267404936268121673493, 10, -2.9128254522781560353905661031604, -0.71794699286267404936268121673493, 10, -2.2455322445133027997599128866568, -1.9893679747223862452898401897983, 10, -1.0638146611276062536433073546505, -2.8050487280562448333398606337141, 10, 0.3616100407659696425177742185042, -2.9781266222941620291919662122382, 10, 1.7041942401934677686625718706637, -2.4689515976809692787696803861763, 10, 2.656368076959629842548338274355, -1.3941695161313054640572772768792, 10;
    tv2 << 0, 0, 0, 0, 0, 3, 0, 3, 0, 0, 3, 3, 3, 0, 0, 3, 0, 3, 3, 3, 0, 3, 3, 3;
    tf1 << 8, 7, 5, 5, 0, 8, 6, 5, 7, 3, 5, 4, 10, 9, 8, 12, 11, 10, 0, 12, 10, 1, 3, 2, 3, 1, 0, 3, 0, 5, 10, 8, 0, 18, 20, 21, 21, 13, 18, 20, 18, 19, 17, 18, 16, 21, 22, 23, 23, 24, 25, 23, 25, 13, 15, 16, 14, 13, 14, 16, 18, 13, 16, 13, 21, 23, 0, 1, 13, 1, 2, 14, 2, 3, 15, 3, 4, 16, 4, 5, 17, 5, 6, 18, 6, 7, 19, 7, 8, 20, 8, 9, 21, 9, 10, 22, 10, 11, 23, 11, 12, 24, 12, 0, 25, 1, 14, 13, 2, 15, 14, 3, 16, 15, 4, 17, 16, 5, 18, 17, 6, 19, 18, 7, 20, 19, 8, 21, 20, 9, 22, 21, 10, 23, 22, 11, 24, 23, 12, 25, 24, 0, 13, 25;
    tf2 << 5, 4, 6, 6, 7, 5, 6, 2, 3, 3, 7, 6, 3, 1, 5, 5, 7, 3, 4, 0, 2, 2, 6, 4, 1, 0, 4, 4, 5, 1, 2, 0, 1, 1, 3, 2;
      Matrix<double,Dynamic,3> v3;
      Matrix<int,Dynamic,3> f3;
      VectorXi j;
      igl::copyleft::cgal::mesh_boolean(tv1,tf1,tv2,tf2,igl::MeshBooleanType(0),
    	  v3, f3, j);
      std::cout << "ok, computed " << v3.rows() << " vertices and " << f3.rows() << " faces\n";
    // std::cout << v3 << "\n" << f3 << "\n" << j << "\n";
    }/*>>>*/

int igl_mesh_boolean(
  int op,
  int nv1, int nf1, const double *mv1, const int *mf1,
  int nv2, int nf2, const double *mv2, const int *mf2,
  int *nv3, int *nf3, double **mv3, int **mf3, int **index) {

  test();
  // IDENTICAL TO test():
  test_values;
  std::cout << "before v3\n";
  return 0;
  Matrix<double,Dynamic,3> v3;
  std::cout << "before f3\n";
  Matrix<int,Dynamic,3> f3;
  std::cout << "before j\n";
  VectorXi j;
  std::cout << "before mesh_bool (test)\n";
  igl::copyleft::cgal::mesh_boolean(tv1,tf1,tv2,tf2,igl::MeshBooleanType(0),
	  v3, f3, j);
  std::cout << "ok, computed " << v3.rows() << " vertices and " << f3.rows() << " faces\n";
  return 0;

  //
  Mesh m1(nv1, nf1, mv1, mf1),
       m2(nv2, nf2, mv2, mf2),
       m3;
  igl::copyleft::cgal::mesh_boolean(tv1,tf1,tv2,tf2,igl::MeshBooleanType(0),
	  v3, f3, j);
  std::cout << "ok, computed " << v3.rows() << " vertices and " << f3.rows() << " faces\n";
  std::cout << "passed with flying colors!\n";

  std::cout << " // v1 == test: " << (m1.v == tv1) << "\n";
  std::cout << " // v2 == test: " << (m2.v == tv2) << "\n";
  std::cout << " // f1 == test: " << (m1.f == tf1) << "\n";
  std::cout << " // f2 == test: " << (m2.f == tf2) << "\n";
  std::cout << "// igl_mesh_boolean:\n";
  std::cout
    << "Matrix<double, Dynamic, 3> v1("<<nv1<<", 3), v2("<<nv2<<",3);\n"
    << "Matrix<int, Dynamic, 3> f1("<<nf1<<", 3), f2("<<nf2<<",3);\n"
    << "v1" << m1.v.format(CommaInitFmt) << "\n"
    << "v2" << m2.v.format(CommaInitFmt) << "\n"
    << "f1" << m1.f.format(CommaInitFmt) << "\n"
    << "f2" << m2.f.format(CommaInitFmt) << "\n\n"
    << "Matrix<double,Dynamic,3> v3;\n"
    << "Matrix<int,Dynamic,3> f3;\nVectorXi j;\n\n"
    << "igl::copyleft::cgal::mesh_boolean(v1,f1,v2,f2,"<<igl::MeshBooleanType(op)<<"v3, f3, j);\n";/*>>>*/
/*>>>*/
  std::cout << "matrix m1.v:\n" << m1.v << "\nmatrix m1.f:" << m1.f << "\n";
  std::cout << "matrix m2.v:\n" << m2.v << "\nmatrix m2.f:" << m2.f << "\n";
  std::cout << "calling mesh_boolean\n";
  igl::copyleft::cgal::mesh_boolean(tv1,tf1,tv2,tf2,igl::MeshBooleanType(op),
            v3, f3, j);
  std::cout << "more...\n";
  igl::copyleft::cgal::mesh_boolean(tv1,tf1, tv2, tf2,
    igl::MeshBooleanType(op), m3.v,m3.f, j);
  std::cout << "1st mesh_boolean OK\n";
  std::cout << std::flush;
  std::cout << "mesh_boolean (2nd)\n";
  return 0;
//   igl::copyleft::cgal::mesh_boolean(m1.v,m1.f,m2.v,m2.f,
//     igl::MeshBooleanType(op), m3.v,m3.f, j);
  std::cout << "done!\n";
//   std::cout << j << "\n";

  *index = (int *) malloc(m3.f.rows()*sizeof(int));
  if (*index == 0) {
    return -1;
  }
  for(int i = 0; i < m3.f.rows(); i++) {
    (*index)[i] = j(i)+1;
  }
  to_jl(m3, nv3, nf3, mv3, mf3);
  if(mv3 == 0) {
    free(*index);
    return -1;
  }
  std::cout << "copy operations done\n";
  return 0;
}

int igl_mesh_is_pwn(int nv, int nf, const double *mv, const int *mf) {
  Mesh m(nv, nf, mv, mf);
//   Matrix<double, Dynamic, 3> v = Map<vertices_t> (mv, nv, 3);
//   Matrix<int, Dynamic, 3> f = Map<faces_t> (mf, nf, 3);
  return (int)igl::copyleft::cgal::piecewise_constant_winding_number(m.v, m.f);
}
