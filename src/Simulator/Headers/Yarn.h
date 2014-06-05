//
//  Yarn.h
//  Visualizer
//
//  Created by eschweickart on 2/12/14.
//
//

#ifndef __Visualizer__Yarn__
#define __Visualizer__Yarn__

#include "Segment.h"
#include <vector>

typedef Eigen::Vector2f Vec2f;

/// A Yarn at a specific point in time.
struct YarnStr
{
public:
  /// Control points (n+1 total) that define the yarn's position.
  std::vector<CtrlPoint> points;
  /// Segments (n total) that define the yarn.
  std::vector<Segment> segments;
};

/// A Yarn stepped though points in time.
class Yarn
{
private:
  float r = constants::radius;
  float shearModulus = constants::shearModulus;
  float youngsModulus = constants::youngsModulus;
  
  float xArea = constants::pi * r * r;
  float tCoeff = xArea * shearModulus * r * r / 2.0f;
  float sCoeff = xArea * youngsModulus;
  // Young's modulus is too strong for the following formula to be accurate for both bending and
  // stretching, so just use 1.0f instead.
  float bCoeff = 1.0f; // xArea * youngsModulus * r * r / 4.0f;
  
  /// The yarn at rest.
  YarnStr  restYS;
  /// The yarn at the current point in time.
  YarnStr* curYS;
  /// The yarn at the next point in time.
  YarnStr* nextYS;
  /// The rest Voronoi lengths defined at each internal control point of the yarn.
  std::vector<float> rvl;
  /// The rest curvature with respect to the previous edge (﻿﻿omega bar ^i _i)
  std::vector<Vec2f> rcp;
  /// The rest curvature with respect to the next edge (omega bar ^i _i+1)
  std::vector<Vec2f> rcn;
  /// The rest curvature defined at each internal control point of the yarn.
  std::vector<Vec2f> rc;
public:
  /// Constructs a yarn that is initially at rest given a vector of initial positions.
  /// The U vector for the first segment is then propagated down the yarn by parallel transport
  /// through space.
  Yarn(std::vector<Vec3f>& points, Vec3f u0) {
    for (Vec3f p : points) {
      CtrlPoint cp;
      cp.pos = p;
      cp.vel.setZero();
      cp.accel.setZero();
      restYS.points.push_back(cp);
    }
    curYS  = new YarnStr(restYS);
    nextYS = new YarnStr(restYS);
    for (int i=0; i<numCPs()-1; i++) {
      Segment restSeg(restYS.points[i], restYS.points[i+1], u0);
      restYS.segments.push_back(restSeg);
      
      Segment curSeg(curYS->points[i], curYS->points[i+1], u0);
      curYS->segments.push_back(curSeg);
      
      Segment nextSeg(nextYS->points[i], nextYS->points[i+1], u0);
      nextYS->segments.push_back(nextSeg);
    }
    for (int i=1; i<numSegs(); i++) {
      restYS.segments[i].parallelTransport(restYS.segments[i-1]);
      curYS->segments[i].setU(restYS.segments[i].getU());
      nextYS->segments[i].setU(restYS.segments[i].getU());
    }
    
    for (int i=1; i<numCPs()-1; i++) {
      const Segment& ePrev = restYS.segments[i-1];
      const Segment& eNext = restYS.segments[i];
      
      Vec3f curveBinorm = 2*ePrev.vec().cross(eNext.vec()) /
      (ePrev.length()*eNext.length() + ePrev.vec().dot(eNext.vec()));
      
      CHECK_NAN_VEC(curveBinorm);
      
      Vec2f restMatCurvePrev(curveBinorm.dot(ePrev.m2()), -(curveBinorm.dot(ePrev.m1())));
      Vec2f restMatCurveNext(curveBinorm.dot(eNext.m2()), -(curveBinorm.dot(eNext.m1())));
      Vec2f restMatCurve = 0.5f*(restMatCurvePrev + restMatCurveNext);
      
      rcp.push_back(restMatCurvePrev);
      rcn.push_back(restMatCurveNext);
      rc.push_back(restMatCurve);
      rvl.push_back(0.5f*(ePrev.length()+eNext.length()));
    }
  }
  
  // WARNING: The following methods are safe as long as curYS and nextYS are always allocated upon
  // construction and are never deleted until destruction.
  /// Get a const reference to the yarn at the current point in time.
  const YarnStr& cur() const  { return *curYS;  }
  /// Get a reference to the yarn at the next point in time.
  YarnStr& next() { return *nextYS; }
  /// Get a const reference to the yarn at the next point in time.
  const YarnStr& next() const { return *nextYS; }
  /// Get a const reference to the yarn at rest.
  const YarnStr& rest() const { return restYS; }
  
  /// Get the number of control points on the yarn.
  const size_t inline numCPs() const { return restYS.points.size(); }
  /// Get the number of control points associated with 2 edges.
  const size_t inline numIntCPs() const { return std::max((int)restYS.points.size()-2, 0); }
  /// Get the number of segments in the yarn.
  const size_t inline numSegs() const { return restYS.segments.size(); }
  
  /// Swaps the current and next yarns, e.g. at the end of a timestep.
  void inline swapYarns() {
    YarnStr* temp = curYS;
    curYS = nextYS;
    nextYS = temp;
  }
  
  // Get the rest Voronoi length for an internal control point.
  const float inline restVoronoiLength(size_t index) const {
    assert(index > 0 && "Voronoi length undifined at this control point.");
    return rvl[index-1];
  }
  
  // Get the rest curvature for an internal control point.
  const Vec2f inline restCurvePrev(size_t index) const {
    assert(index > 0 && "Curvature undefined at this control point.");
    return rcp[index-1];
  }
  
  // Get the rest curvature for an internal control point.
  const Vec2f inline restCurveNext(size_t index) const {
    assert(index > 0 && "Curvature undefined at this control point.");
    return rcn[index-1];
  }
  
  // Get the rest curvature for an internal control point.
  const Vec2f inline restCurve(size_t index) const {
    assert(index > 0 && "Curvature undefined at this control point.");
    return rc[index-1];
  }
  
  /// Get the yarn's radius
  const float inline radius() const { return r; }
  
  /// Get the twist coefficient
  const float inline twistCoeff() const { return tCoeff; }
  
  /// Get the stretch coefficient
  const float inline stretchCoeff() const { return sCoeff; }
  
  /// Get the bend coefficient (multiply by the identity matrix to get the bending matrix)
  const float inline bendCoeff() const { return bCoeff; }
  
  ~Yarn() {
    delete curYS;
    delete nextYS;
  }
};

class Spline {
  typedef Eigen::Vector4f Vec4f;
  Vec3f p[4];
  
public:
  
  Spline(const Vec3f p0, const Vec3f p1, const Vec3f p2, const Vec3f p3) {
    p[0] = p0; p[1] = p1; p[2] = p2; p[3] = p3;
  }
  
  Spline(const CtrlPoint& cp0, const CtrlPoint& cp1, const CtrlPoint& cp2, const  CtrlPoint& cp3) :
  Spline(cp0.pos, cp1.pos, cp2.pos, cp3.pos) { }
  
  Spline(const CtrlPoint& cp0, const CtrlPoint& cp1, const CtrlPoint&cp2, bool start) :
  Spline(start ? cp0.pos : (2*cp0.pos) - cp1.pos, start ? cp1.pos : cp0.pos,
         start ? cp2.pos : cp1.pos, start ? (2*cp2.pos) - cp1.pos : cp2.pos) { }

  Vec3f eval(float t) const {
    Vec4f u(t*t*t, t*t, t, 1);
    Vec4f uBasis(u.dot(constants::basis[0]),
                 u.dot(constants::basis[1]),
                 u.dot(constants::basis[2]),
                 u.dot(constants::basis[3]));
    Vec4f x(p[0].x(), p[1].x(), p[2].x(), p[3].x());
    Vec4f y(p[0].y(), p[1].y(), p[2].y(), p[3].y());
    Vec4f z(p[0].z(), p[1].z(), p[2].z(), p[3].z());
    return Vec3f(uBasis.dot(x), uBasis.dot(y), uBasis.dot(z));
  }
  
};

#endif /* defined(__Visualizer__Yarn__) */