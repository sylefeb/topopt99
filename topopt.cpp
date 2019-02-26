// --------------------------------------------------------------
// SL 2019-02-25
// A simple, easily hackable topopt code, directly inspired
// from "A 99 line topology optimization code written in MATLAB"
// MIT-license
// (c) Sylvain Lefebvre, https://github.com/sylefeb
// --------------------------------------------------------------
/*

This code is directly inspired from the fantastic paper

A 99 line topology optimization code written in MATLAB
Structural and Multidisciplinary Optimization 21(2), 2001, pp. 120-127
by Ole Sigmund.

http://www.topopt.mek.dtu.dk/Apps-and-software/A-99-line-topology-optimization-code-written-in-MATLAB

I simply re-implemented it in C/C++ using LibSL-small and Eigen.
And yes, it takes more than 99 lines of C/C++ ;-)

Notations
---------

- dens : optimized densities
- nelx : number of elements along X axis (nelx+1 corners) 
- nely : number of elements along Y axis (nely+1 corners)
- dc   : compliance gradient
- KE   : stiffness matrix for a single element

*/
// --------------------------------------------------------------

#include <LibSL/LibSL.h>
#include <Eigen/Sparse>

#include <iostream>
#include <ctime>
#include <cmath>
#include <set>
#include <limits>

// --------------------------------------------------------------

using namespace std;

// --------------------------------------------------------------

LIBSL_WIN32_FIX;

// --------------------------------------------------------------

// stiffness matrix for a unit square element
Array2D<double> KE;

// Computes the stiffness matrix of a single, unit square element
// (forward declaration, code at the end)
void lk(Array2D<double>& _KE);

// --------------------------------------------------------------

// matrix vector multiply
template <typename T>
void mv_mul(const Array2D<T>& A, const Array<T>& v, Array<T>& _res)
{
  _res.allocate(v.size());
  ForIndex(l, A.ysize()) {
    T al = 0;
    ForIndex(c, A.xsize()) {
      al = al + A.at(c, l) * v[c];
    }
    _res[l] = al;
  }
}

// vector-vector multiply
template <typename T>
T vv_mul(const Array<T>& a, const Array<T>& b)
{
  T res = 0;
  ForIndex(i, a.size()) {
    res = res + a[i] * b[i];
  }
  return res;
}

// --------------------------------------------------------------

// Applies optimality criterion (OC) to update densities (x) from compliance gradient (dc).
// Maintains the volume fraction
void OC(int nelx, int nely, Array2D<double>& _dens, double volfrac, const Array2D<double>& dc)
{
  double move = 0.2;
  Array2D<double> dens_new(_dens.xsize(), _dens.ysize());
  double vtot = 0.0;
  double l1 = 0, l2 = 100000;
  while (l2 - l1 > 1e-3) { // bisection search for volume preservation
    double lmid = 0.5*(l2 + l1);
    vtot = 0.0;
    ForArray2D(dens_new, c, l) {
      // OC term
      double Be = (-dc.at(c, l)) / lmid;
      dens_new.at(c, l) = max(0.02,
        max(_dens.at(c, l) - move,
          min(1.0,
            min(_dens.at(c, l) + move, _dens.at(c, l) * sqrt(Be))
          )
        )
      );
      vtot += dens_new.at(c, l);
    }
    if (vtot - volfrac * nelx * nely > 0.0) {
      l1 = lmid;
    } else {
      l2 = lmid;
    }
  }
  cerr << "Vol frac in result  " << vtot / (nelx* nely) << endl;
  _dens = dens_new;
}

// --------------------------------------------------------------

// Filters the result to prevent the 'checkerboard effect' due to the solver
// attempting to produce infinitely small 'bubbles' (composite).
void filter(int nelx, int nely, double rmin, const Array2D<double>& dens, Array2D<double>& _dc)
{
  Array2D<double> dcn(nelx, nely);
  dcn.fill(0);
  ForArray2D(dcn, i, j) {
    double sum = 0.0;
    ForRange(l, max(j - round(rmin), 0), min(j + round(rmin), nely - 1)) {
      ForRange(k, max(i - round(rmin), 0), min(i + round(rmin), nelx - 1)) {
        double fac = rmin - sqrt((double)(i - k)*(i - k) + (double)(j - l)*(j - l));
        sum += max(0.0f, fac);
        dcn.at(i, j) = dcn.at(i, j) + max(0.0f, fac) * dens.at(k, l) * _dc.at(k, l);
      }
    }
    dcn.at(i, j) = dcn.at(i, j) / (dens.at(i, j)*sum);
  }
  _dc = dcn;
}


// --------------------------------------------------------------

// Straightforward linear finite element solver
// The degrees of freedom (variables) are the corners of the elements, times two (x and y coordinates)
// Each corner at (cx,cy) is associated to two ids:
//   - grid id computed cx + (nelx + 1) * cy
//   - variable id which identifies it in the sparse system of equations
// The x/y coordinates are at respectively id*2+0 and id*2+1
void FE(int nelx, int nely, const Array2D<double>& dens, double penal, Array2D<v2f>& _U)
{
  // variables to lock: there are attachement points of the structure
  set<int> locked;
  // -> here we attach the left third of the bottom row (x in  [0,(nelx + 1) / 3 - 1] , y = nely
  ForIndex(x, (nelx + 1) / 3) {
    locked.insert((x + (nelx + 1) * nely) * 2 + 0); // x
    locked.insert((x + (nelx + 1) * nely) * 2 + 1); // y
  }
  // mapping between variable and ids
  Array<int>   grid2var(2 * (nelx + 1)*(nely + 1));
  Array<int>   var2grid(2 * (nelx + 1)*(nely + 1) - (int)locked.size());
  grid2var.fill(-1);
  var2grid.fill(-1);
  int varid = 0;
  ForIndex(gridid, 2 * (nelx + 1)*(nely + 1)) {
    if (locked.find(gridid) == locked.end()) {
      grid2var[gridid] = varid;
      var2grid[varid]  = gridid;
      varid++;
    }
  }
  sl_assert(varid == 2 * (nelx + 1)*(nely + 1) - (int)locked.size());
  // matrix A (FE equations)
  std::vector<Eigen::Triplet<double> > coefficients;
  ForIndex(x, nelx) {
    ForIndex(y, nely) {
      int p00 = (x + (y)* (nelx + 1));
      int p01 = (x + (y + 1) * (nelx + 1));
      int p11 = ((x + 1) + (y + 1) * (nelx + 1));
      int p10 = ((x + 1) + (y)* (nelx + 1));
      Array<int> corners(8);
      corners[0] = p00 * 2 + 0;
      corners[1] = p00 * 2 + 1;
      corners[2] = p10 * 2 + 0;
      corners[3] = p10 * 2 + 1;
      corners[4] = p11 * 2 + 0;
      corners[5] = p11 * 2 + 1;
      corners[6] = p01 * 2 + 0;
      corners[7] = p01 * 2 + 1;
      double pow_x = pow(dens.at(x, y), penal);
      // add coefficients only for non-locked variables
      ForIndex(k, 8) {
        if (grid2var[corners[k]] > -1) {
          ForIndex(l, 8) {
            if (grid2var[corners[l]] > -1) {
              coefficients.push_back(Eigen::Triplet<double>(grid2var[corners[k]], grid2var[corners[l]],
                KE.at(l, k) * pow_x
                ));
            }
          }
        }
      }
    }
  }
  Eigen::SparseMatrix<double> A(
    2 * (nelx + 1)*(nely + 1) - locked.size(),
    2 * (nelx + 1)*(nely + 1) - locked.size());
  A.setFromTriplets(coefficients.begin(), coefficients.end());
  // vector b, contains external forces
  Eigen::VectorXd b = Eigen::VectorXd(2 * (nelx + 1)*(nely + 1) - locked.size());
  ForIndex(i, b.size()) {
    b[i] = 0.0;
  }
  // -> here we apply external forces at two specific points (near the image top)
  b[grid2var[(nelx / 2     + 24 * (nelx + 1)) * 2 + 1]] = -1;
  b[grid2var[(nelx * 5 / 6 + 16 * (nelx + 1)) * 2 + 1]] = -1;
  //                                         ^^^^^^^ force in y direction
  // solver
  cerr << "solving ...";
  Eigen::SparseLU<Eigen::SparseMatrix<double> > solver(A);
  Eigen::VectorXd result = solver.solve(b);
  cerr << " done.\n";
  // store computed displacement
  _U.allocate(nelx + 1, nely + 1);
  ForIndex(varid, result.size()) {
    int gridid = var2grid[varid];
    int x = (gridid / 2) % (nelx + 1);
    int y = (gridid / 2) / (nelx + 1);
    int c = (gridid & 1);
    _U.at(x, y)[c] = (float)result[varid];
  }
}

// --------------------------------------------------------------

// Executes the global optimization loop until 'convergence' (change below threshold)
// volfrac is the target volume fraction in [0-1]
// penal   is the density penality (typically 3)
// rmin    controls the filter size
void topopt(int nelx, int nely, double volfrac, double penal, double rmin)
{
  ImageFloat1     img(nelx, nely);  // for image output
  Array2D<double> dens(nelx, nely); // optimized density grid
  Array2D<double> dens_old;         // result from previous iteration
  Array2D<double> dc(nelx, nely);   // compliance gradient

  // initialization with a random field of selected volume fraction (volfrac)
  double tot = 0;
  ForArray2D(dens, i, j) {
    dens.at(i, j) = rnd();
    tot += dens.at(i, j);
  }
  tot = tot / (nelx*nely);
  ForArray2D(dens, i, j) {
    dens.at(i, j) = volfrac * dens.at(i, j) / tot;
  }
  int loop = 0;
  double change = 1;
  while (change > 0.01) {

    loop++;
    dens_old = dens;

    // solve for displacement (finite element solution)
    Array2D<v2f> U;
    FE(nelx, nely, dens, penal, /*out*/ U);

    // option to output displacement magnitudes
    if (0) {
      ImageFloat1 img;
      img.pixels().allocate(nelx + 1, nely + 1);
      ForImage((&img), i, j) {
        img.pixel(i, j) = length(U.at(i, j));
      }
      img.remap(0.0f, 255.0f);
      saveImage(sprint("%03d_FE.tga", loop), img.cast<ImageL8>());
    }

    // compute compliance and gradient
    // -> for each element
    double c_comp = 0.0;
    ForIndex(ely, nely) {
      ForIndex(elx, nelx) {
        // get displacement for the element
        Array<double> corners(8);
        corners[0] = U.at(elx, ely)[0];
        corners[1] = U.at(elx, ely)[1];
        corners[2] = U.at(elx + 1, ely)[0];
        corners[3] = U.at(elx + 1, ely)[1];
        corners[4] = U.at(elx + 1, ely + 1)[0];
        corners[5] = U.at(elx + 1, ely + 1)[1];
        corners[6] = U.at(elx, ely + 1)[0];
        corners[7] = U.at(elx, ely + 1)[1];
        // compute compliance
        Array<double> KE_Ue;
        mv_mul(KE, corners, /*out*/KE_Ue);
        v2d Ue_KE_Ue = vv_mul(KE_Ue, corners);        
        double compliance = Ue_KE_Ue[0] + Ue_KE_Ue[1];
        double comp       = pow(dens.at(elx, ely), penal) * compliance;
        c_comp += comp;
        // gradient
        double dcomp    = -penal * pow(dens.at(elx, ely), penal - 1.0f) * compliance;
        dc.at(elx, ely) = dcomp;
      }
    }

    // filtering
    filter(nelx, nely, rmin, dens, dc/*out*/);

    // OC method
    OC(nelx, nely, dens/*out*/, volfrac, dc);

#if 0
    // For fun: kills a circle of density (obstacle)
    ForArray2D(dens, i, j) {
      if ( length(v2f((float)i, (float)j) - v2f((float)nelx/2, (float)nely/2)) < nelx/6.0f) {
        dens.at(i, j) = 0.01;
      }
    }
#endif

    // compute max change
    change = 0.0;
    ForArray2D(dens, i, j) {
      change = max(change, abs(dens.at(i, j) - dens_old.at(i, j)));
    }
    // output iteration stats
    cerr << sprint(" Loop %d, compliance = %f, change = %f\n", loop, c_comp, change);
    // output image
    ForImage((&img), i, j) {
      img.pixel(i,j) = 255.0f * (1.0f - (float)pow(dens.at(i,j),penal));
    }
    saveImage(sprint("%03d_struct.tga", loop), img.cast<ImageL8>());

  }

}

// --------------------------------------------------------------

// Program entry point.
int main(int argc, char **argv)
{
  try {

    // generate stiffness matrix for unit element
    lk(KE);

    // go ahead!
    topopt(256, 96, 0.3f, 3.0f, 1.2f);

  } catch (Fatal& e) {
    cerr << Console::red << e.message() << Console::gray << endl;
    return (-1);
  }

  return (0);
}

// --------------------------------------------------------------

// Computes the stiffness matrix of a single, unit square element
void lk(Array2D<double>& _KE)
{
  double E = 1.0f;
  double nu = 0.3f;
  double k[8];
  k[0] = 1.0f / 2.0 - nu / 6.0;
  k[1] = 1.0f / 8.0 + nu / 8.0;
  k[2] = -1.0f / 4.0 - nu / 12.0;
  k[3] = -1.0f / 8.0 + 3.0*nu / 8.0;
  k[4] = -1.0f / 4.0 + nu / 12.0;
  k[5] = -1.0f / 8.0 - nu / 8.0;
  k[6] = nu / 6.0;
  k[7] = 1.0f / 8.0 - 3.0*nu / 8.0;
  _KE.allocate(8, 8);

  _KE.at(0, 0) = E / (1.0 - nu * nu) * k[0];
  _KE.at(1, 0) = E / (1.0 - nu * nu) * k[1];
  _KE.at(2, 0) = E / (1.0 - nu * nu) * k[2];
  _KE.at(3, 0) = E / (1.0 - nu * nu) * k[3];
  _KE.at(4, 0) = E / (1.0 - nu * nu) * k[4];
  _KE.at(5, 0) = E / (1.0 - nu * nu) * k[5];
  _KE.at(6, 0) = E / (1.0 - nu * nu) * k[6];
  _KE.at(7, 0) = E / (1.0 - nu * nu) * k[7];

  _KE.at(0, 1) = E / (1.0 - nu * nu) * k[1];
  _KE.at(1, 1) = E / (1.0 - nu * nu) * k[0];
  _KE.at(2, 1) = E / (1.0 - nu * nu) * k[7];
  _KE.at(3, 1) = E / (1.0 - nu * nu) * k[6];
  _KE.at(4, 1) = E / (1.0 - nu * nu) * k[5];
  _KE.at(5, 1) = E / (1.0 - nu * nu) * k[4];
  _KE.at(6, 1) = E / (1.0 - nu * nu) * k[3];
  _KE.at(7, 1) = E / (1.0 - nu * nu) * k[2];

  _KE.at(0, 2) = E / (1.0 - nu * nu) * k[2];
  _KE.at(1, 2) = E / (1.0 - nu * nu) * k[7];
  _KE.at(2, 2) = E / (1.0 - nu * nu) * k[0];
  _KE.at(3, 2) = E / (1.0 - nu * nu) * k[5];
  _KE.at(4, 2) = E / (1.0 - nu * nu) * k[6];
  _KE.at(5, 2) = E / (1.0 - nu * nu) * k[3];
  _KE.at(6, 2) = E / (1.0 - nu * nu) * k[4];
  _KE.at(7, 2) = E / (1.0 - nu * nu) * k[1];

  _KE.at(0, 3) = E / (1.0 - nu * nu) * k[3];
  _KE.at(1, 3) = E / (1.0 - nu * nu) * k[6];
  _KE.at(2, 3) = E / (1.0 - nu * nu) * k[5];
  _KE.at(3, 3) = E / (1.0 - nu * nu) * k[0];
  _KE.at(4, 3) = E / (1.0 - nu * nu) * k[7];
  _KE.at(5, 3) = E / (1.0 - nu * nu) * k[2];
  _KE.at(6, 3) = E / (1.0 - nu * nu) * k[1];
  _KE.at(7, 3) = E / (1.0 - nu * nu) * k[4];

  _KE.at(0, 4) = E / (1.0 - nu * nu) * k[4];
  _KE.at(1, 4) = E / (1.0 - nu * nu) * k[5];
  _KE.at(2, 4) = E / (1.0 - nu * nu) * k[6];
  _KE.at(3, 4) = E / (1.0 - nu * nu) * k[7];
  _KE.at(4, 4) = E / (1.0 - nu * nu) * k[0];
  _KE.at(5, 4) = E / (1.0 - nu * nu) * k[1];
  _KE.at(6, 4) = E / (1.0 - nu * nu) * k[2];
  _KE.at(7, 4) = E / (1.0 - nu * nu) * k[3];

  _KE.at(0, 5) = E / (1.0 - nu * nu) * k[5];
  _KE.at(1, 5) = E / (1.0 - nu * nu) * k[4];
  _KE.at(2, 5) = E / (1.0 - nu * nu) * k[3];
  _KE.at(3, 5) = E / (1.0 - nu * nu) * k[2];
  _KE.at(4, 5) = E / (1.0 - nu * nu) * k[1];
  _KE.at(5, 5) = E / (1.0 - nu * nu) * k[0];
  _KE.at(6, 5) = E / (1.0 - nu * nu) * k[7];
  _KE.at(7, 5) = E / (1.0 - nu * nu) * k[6];

  _KE.at(0, 6) = E / (1.0 - nu * nu) * k[6];
  _KE.at(1, 6) = E / (1.0 - nu * nu) * k[3];
  _KE.at(2, 6) = E / (1.0 - nu * nu) * k[4];
  _KE.at(3, 6) = E / (1.0 - nu * nu) * k[1];
  _KE.at(4, 6) = E / (1.0 - nu * nu) * k[2];
  _KE.at(5, 6) = E / (1.0 - nu * nu) * k[7];
  _KE.at(6, 6) = E / (1.0 - nu * nu) * k[0];
  _KE.at(7, 6) = E / (1.0 - nu * nu) * k[5];

  _KE.at(0, 7) = E / (1.0 - nu * nu) * k[7];
  _KE.at(1, 7) = E / (1.0 - nu * nu) * k[2];
  _KE.at(2, 7) = E / (1.0 - nu * nu) * k[1];
  _KE.at(3, 7) = E / (1.0 - nu * nu) * k[4];
  _KE.at(4, 7) = E / (1.0 - nu * nu) * k[3];
  _KE.at(5, 7) = E / (1.0 - nu * nu) * k[6];
  _KE.at(6, 7) = E / (1.0 - nu * nu) * k[5];
  _KE.at(7, 7) = E / (1.0 - nu * nu) * k[0];
}

// --------------------------------------------------------------
