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
  alpha1 = 0.0;
  // Mass coefficient
  alpha2 = 0.0;
  
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
                    const std::function<complex(complex)>& f, VecXe& z) {
  // First, check if input is zero. If it is, just return zero.
  if (v.norm() == 0.0) {
    z = VecXe::Zero(v.rows());
    return;
  }
  
  // Otherwise, perform a the functions on a Krylov subspace.
  MatXe Q(v.rows(), m+1);
  MatXe H(m, m);
  arnoldi(A, v, Q, H); // sets Q and H. Q now has dimensions (n x m).
  
  if (H.isApprox(H.transpose())) {
    Eigen::SelfAdjointEigenSolver<MatXe> saes(H);
    // z = T^T * f(D) * T * e1. z is now (m x 1).
    z.noalias() = saes.eigenvectors() * (saes.eigenvalues()
                                         .unaryExpr(f)
                                         .cwiseProduct(saes.eigenvectors().row(0).transpose())).real();
  } else {
    Eigen::EigenSolver<MatXe> es(H);
    // zc = T^T * f(D) * T * e1. zc is now (m x 1).
    Eigen::Matrix<complex, Eigen::Dynamic, 1> zc;
    zc.noalias() = es.eigenvectors() * (es.eigenvalues()
                                        .unaryExpr(f)
                                        .cwiseProduct(es.eigenvectors().row(0).transpose()));
    // z is the real part of zc (hopefully zc has no imaginary component??). z is now (m x 1).
    z = zc.real();
  }
  // z = ||v|| * Q * z. z is now (n x 1).
  z = v.norm() * (Q * z);
  CHECK_NAN_VEC(z);
}

bool ExpIntegrator::integrate(Clock& c) {
  bool argFiltering = false; // TODO: what does this do?
  std::size_t krylovBasisSize = 15;
  
  real h = c.timestep();
  // Choose (psi(·), phi(·)) = (sinc^2(·), sinc(·)).
  std::function<complex(complex)> psi = [h] (complex k) {
    complex x = h * std::sqrt(k);
    x = (x == 0.0) ? 1.0 : sin(x) / x;
    return x * x;
  };
  std::function<complex(complex)> phi = [h] (complex k) {
    complex x = h * std::sqrt(k);
    return (x == 0.0) ? 1.0 : sin(x) / x;
  };
  std::function<complex(complex)> cosSqrt = [h] (complex k) {
    return cos(h * std::sqrt(k));
  };
  
  if (c.getTicks() % 100 == 0) { // Periodically linearize the system (sets damping and stiffness)
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
