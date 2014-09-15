//
//  ExpIntegrator.h
//  Visualizer
//
//  Created by eschweic on 9/3/14.
//
//

#ifndef __Visualizer__ExpIntegrator__
#define __Visualizer__ExpIntegrator__

#include "Integrator.h"

#endif /* defined(__Visualizer__ExpIntegrator__) */

class ExpIntegrator : public Integrator {
  Eigen::SparseMatrix<real> damping;
  Eigen::SparseMatrix<real> stiffness;
  Eigen::SparseMatrix<real> A; // cached InvMass * stiffness
  real alpha1;
  real alpha2;
  void setStiffness();
  void setDamping();
  
  void static testArnoldi();
  void static arnoldi(const Eigen::SparseMatrix<real>& A, const VecXe& v, MatXe& Q, MatXe& H);
  void static matFunc(const Eigen::SparseMatrix<real>& A, const VecXe& v, const std::size_t m,
                      const std::function<real(real)>& f, VecXe& z);
  
public:
  ExpIntegrator(Rod&, std::vector<RodEnergy*>&);
  bool integrate(Clock&);
  void draw();
};