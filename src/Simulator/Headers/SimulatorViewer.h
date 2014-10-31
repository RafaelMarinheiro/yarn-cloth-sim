#ifndef SIMULATOR_VIEWER_H
#define SIMULATOR_VIEWER_H

#include "Defines.h"

#include <QGLViewer/qglviewer.h>

#include "Resources.h"
#include "Util.h"
#include "Rod.h"
#include "Clock.h"
#include "Integrator.h"
#include "IMEXIntegrator.h"
#include "ConstraintIntegrator.h"
#include "YarnBuilder.h"

class SimulatorViewer : public QGLViewer{
public:
	virtual void init();
	virtual void animate();
	virtual void draw();

	virtual void keyPressEvent(QKeyEvent * event);
	virtual void mousePressEvent(QMouseEvent * event);
	virtual void mouseMoveEvent(QMouseEvent * event);
	virtual void mouseReleaseEvent(QMouseEvent * event);

private:

	void loadRodFile(std::string filename);
	void loadDefaultRod(int numPoints);
	void loadStdEnergies();
	void loadStdEnergiesAndConsts();

	Vec3e eyePos = Vec3e(50.0, 0.0, 0.0);
  	Vec3e targetPos = Vec3e(0.0, 0.0, 0.0);

	bool running = true;
	Rod * r;
	Clock c;
	Integrator * integrator = nullptr;
	ConstraintIntegrator * cIntegrator = nullptr;
	std::vector<RodEnergy*> energies;
	MouseSpring * mouseSpring;
	std::vector<RodConstraint*> constraints;

	real twist = 0.0;
	real rodTwist = 0.0;
	int numRodTwists = 0;

	//Interactive stuff;

	bool isMouseDown = false;
	Vec3e mousePosition;
	bool isRotate = false;
	bool isSetup = false;
};



#endif // SIMULATOR_VIEWER_H
