//
//  Test.h
//  Visualizer
//
//  Created by eschweic on 9/19/14.
//
//

#ifndef Visualizer_Test_h
#define Visualizer_Test_h

#include "Defines.h"

#include "Rod.h"
#include "Energy.h"
#include "ExpIntegrator.h"

/// Builds a small test system.
static void test() {
  
  // Build a specific rod and test various quantities for verification purposes.
  VecXe pos(6);
  pos << 0.0, 0.0, 0.0, 0.0, 1.0, 0.0; // , 0.0, 2.0, 0.0;
  
  Vec3e u0(1.0, 0.0, 0.0);
  
  VecXe masses(2);
  masses << 1.0, 1.0; // , 1.0;
  
  Rod r(pos, u0, &masses, 1.0, 1.0);
  
  real h = 0.01;
  Clock c(h);
  
  std::vector<RodEnergy*> energies;
  
//  RodEnergy* bending = new Bending(r, Explicit);
//  energies.push_back(bending);
  
  RodEnergy* stretching = new Stretching(r, Explicit);
  energies.push_back(stretching);
  
  r.next().POS(1) += Vec3e(0.0, 0.01, 0.0);
//  r.next().vel = (r.next().pos - r.cur().pos) / h;
//  r.next().dVel = r.next().vel - r.cur().vel;
  
  r.next().updateReferenceFrames(r.cur());
  
  r.swapRods();
  
  VecXe f = VecXe::Zero(r.numDOF());
  std::vector<Triplet> triplets;
  Eigen::SparseMatrix<real> K(r.numDOF(), r.numDOF());
  
  stretching->eval(&f, &triplets);
  K.setFromTriplets(triplets.begin(), triplets.end());
  
  std::cout << "forces:\n" << f << "\n\n";
  std::cout << "stiffness:\n" << K.toDense() << "\n\n";
  
  ExpIntegrator integrator(r, energies);
  
  for (int i=0; i<100; i++) {
    integrator.integrate(c);
    r.next().updateReferenceFrames(r.cur());
    r.swapRods();
    c.increment();
    
    std::cout << i << ":\n" << r.cur().pos << "\n\n";
  }

  
  
}

/// Builds and tests a simple mass-spring system for 100 ticks.
static void testInt() {
  VecXe pos(3);
  pos << 0.0, 1e-15, 0.0;
  
  VecXe masses(1);
  masses << 0.01;
  
  Rod r(pos, Vec3e::Zero(), &masses, 1.0, 1.0);
  
  Clock c(0.01);
  
  std::vector<RodEnergy*> energies;

  Spring* spring = new Spring(r, Explicit, 0, 100);
  Vec3e clamp = Vec3e::Zero();
  spring->setClamp(clamp);
  energies.push_back(spring);
  
  ExpIntegrator integrator(r, energies);
  
  for (int i=0; i<100; i++) {
    integrator.integrate(c);
    r.swapRods();
    c.increment();
    
    std::cout << i << ":\n" << r.cur().POS(0) << "\n\n";
  }
}


#endif
