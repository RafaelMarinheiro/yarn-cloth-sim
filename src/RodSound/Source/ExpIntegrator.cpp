//
//  ExpIntegrator.cpp
//  Visualizer
//
//  Created by eschweic on 9/3/14.
//
//

#include "ExpIntegrator.h"

ExpIntegrator::ExpIntegrator(Rod& r, std::vector<RodEnergy*>& energies) : Integrator(r, energies) {
  // Stiffness coefficient
  alpha1 = 1.0e-9;
  // Mass coefficient
  alpha2 = 2.0;
  
  stiffness.resize(r.numDOF(), r.numDOF());
  setDamping();
}

/// Perform Arnoldi iteration to obtain a Krylov subspace. Results are written to Q and H; H is
/// assumed to have the correct dimensions already, while the dimensions of Q are set in this
/// method.
void ExpIntegrator::arnoldi(const Eigen::SparseMatrix<real>& A, const VecXe& v, MatXe& Q, MatXe& H) {
  // Assume H is (m x m), where m is the desired length of the subspace.
  int m = H.rows();
  H.setZero();
  // Set Q to be (n x m+1), where A is (n x n) and v is (n x 1).
  Q.resize(A.rows(), m+1);
  Q.col(0) = v.normalized();
  
  CHECK_NAN_VEC(v);
  CHECK_NAN_VEC(Q.col(0));
  
  for (int k=1; k<m+1; k++) {
    Q.col(k).noalias() = A * Q.col(k-1);
    for (int j=0; j<k; j++) {
      H(j, k-1) = Q.col(j).dot(Q.col(k));
      Q.col(k) -= H(j, k-1) * Q.col(j);
    }
    real norm = Q.col(k).norm();
    if (norm == 0.0) { // Check for convergence--if so, stop.
      for (int i=k; i<m; i++) {
        Q.col(i).setZero();
      }
      // Shrink Q to be (n x m)
      Q = Q.block(0, 0, A.rows(), m).eval();
      return;
    }
    if (k < m) {
      H(k, k-1) = norm;
      Q.col(k) /= norm;
    } else {
      // Shrink Q to be (n x m)
      Q = Q.block(0, 0, A.rows(), m).eval();
    }
  }
}

/// Performs a unitary function on a matrix and multiplies it by a given vector v, writing output
/// to z. This is done by diagonalizing the (m x m) Hessian given by Arnoldi iteration.
void ExpIntegrator::matFunc(const Eigen::SparseMatrix<real>& A, const VecXe& v, const std::size_t m,
                    const std::function<real(real)>& f, VecXe& z) {
  // First, check if input is zero. If it is, just return zero.
  if (v.norm() == 0.0) {
    z = VecXe::Zero(v.rows());
    return;
  }
  
  MatXe Q(v.rows(), m+1);
  MatXe H(m, m);
  arnoldi(A, v, Q, H); // sets Q and H
    
  Eigen::SelfAdjointEigenSolver<MatXe> saes(H);
  // z = T^T * f(D) * T * e1. z is now (m x 1).
  z.noalias() = saes.eigenvectors().transpose() *(saes.eigenvalues()
                                                  .unaryExpr(f)
                                                  .cwiseProduct(saes.eigenvectors().col(0)));
  // z = ||v|| * Q * z. z is now (n x 1).
  z = v.norm() * (Q * z);
  CHECK_NAN_VEC(z);
}

bool ExpIntegrator::integrate(Clock& c) {
  bool argFiltering = false; // TODO: what does this do?
  std::size_t krylovBasisSize = 5;
  
  real h = c.timestep();
  // Choose (psi(·), phi(·)) = (sinc^2(·), sinc(·)).
  std::function<real(real)> psi = [h] (real k) {
    assert(k >= 0.0);
    real x = h * sqrt(k);
    x = (x == 0.0) ? 1.0 : sin(x) / x;
    return x * x;
  };
  std::function<real(real)> phi = [h] (real k) {
    assert(k >= 0.0);
    real x = h * sqrt(k);
    return (x == 0.0) ? 1.0 : sin(x) / x;
  };
  std::function<real(real)> cosSqrt = [h] (real k) {
    assert(k >= 0.0);
    return cos(h * sqrt(k));
  };
  
  if (c.getTicks() % 1000 == 0) { // Periodically reset linearize (sets damping and stiffness)
    setDamping();
  }
  
  VecXe xPhi;
  matFunc(A, r.cur().pos, krylovBasisSize, phi, xPhi); // Sets xPhi
  
  Eigen::SparseMatrix<real> APhi;
  
  if (argFiltering) {
    VecXe del = xPhi - r.cur().pos;
    std::vector<Triplet> filteredStiffness;
    for (RodEnergy* e : energies) {
      if (e->energySource() == Internal) {
        e->eval(nullptr, &filteredStiffness, &del);
      }
    }
    Eigen::SparseMatrix<real> K(A.rows(), A.cols());
    K.setFromTriplets(filteredStiffness.begin(), filteredStiffness.end());
    APhi = -r.getInvMass().sparse * K;
    matFunc(APhi, r.cur().pos, krylovBasisSize, phi, xPhi);
    
    // TODO: Should we update the damping matrix as well?
  } else {
    APhi = A;
  }
  
  VecXe del = xPhi - r.cur().pos;
  VecXe extForcesPhi = VecXe::Zero(r.numDOF());
  for (RodEnergy* e : energies) {
    if (e->energySource() == External) {
      e->eval(&extForcesPhi, nullptr, &del);
    }
  }
  VecXe lambdaPhi = extForcesPhi + damping * r.cur().vel;
  VecXe xCos;
  matFunc(APhi, r.cur().pos, krylovBasisSize, cosSqrt, xCos); // sets xCos
  VecXe xPsi;
  matFunc(APhi, lambdaPhi, krylovBasisSize, psi, xPsi); // sets xPsi
  
  // Update rod configuration
  r.next().pos = 2.0 * xCos + (h * r.cur().vel - r.cur().pos) - h*h*(r.getInvMass().sparse * xPsi);
  r.next().vel = (r.next().pos - r.cur().pos) / h;
  r.next().dVel = r.next().vel - r.cur().vel;
  
  return true;
}

void ExpIntegrator::setStiffness() {
  std::vector<Triplet> triplets;
  for (RodEnergy* e : energies) {
    if (e->energySource() == Internal) {
      e->eval(nullptr, &triplets);
    }
  }
  stiffness.setFromTriplets(triplets.begin(), triplets.end());
  stiffness *= -1.0;
  
  A = r.getInvMass().sparse * stiffness;
}

void ExpIntegrator::setDamping() {
  setStiffness();
  
  // Rayleigh damping matrix
  damping = alpha1 * stiffness + alpha2 * r.getMass().sparse;
}

void ExpIntegrator::draw() {
  
}

void ExpIntegrator::testArnoldi() {
  Eigen::SparseMatrix<real> A(10, 10);
  std::vector<Triplet> tr;
  for (int i=0; i<10; i++) {
    for (int j=0; j<10; j++) {
      tr.push_back(Triplet(i, j, 1.0));
    }
    tr.push_back(Triplet(i, i, i+1.0));
  }
  A.setFromTriplets(tr.begin(), tr.end());
  
  VecXe v = VecXe::Zero(10);
  v(0) = 1.0;
  
  MatXe Q;
  MatXe H(5, 5);
  
  arnoldi(A, v, Q, H);
  
  std::cout << Q << "\n\n" << H << "\n\n\n";
}
