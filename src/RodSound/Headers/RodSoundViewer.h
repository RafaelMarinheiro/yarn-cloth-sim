#ifndef ROD_SOUND_VIEWER_H
#define ROD_SOUND_VIEWER_H

#include "Defines.h"
#include "RodSoundDefines.h"

#include <QGLViewer/qglviewer.h>

#include "Resources.h"
#include "Util.h"
#include "Rod.h"
#include "Clock.h"
#include "Integrator.h"
#include "ExIntegrator.h"
#include "ExpIntegrator.h"
#include "ConstraintIntegrator.h"
#include "YarnBuilder.h"
#include "FrameExporter.h"

class RodSoundViewer : public QGLViewer{
public:
	virtual void init();

	void update();
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

private:
	bool running = true;

	Vec3e eyePos = Vec3e(50.0, 0.0, 0.0);
	Vec3e targetPos = Vec3e(0.0, 0.0, 0.0);

	Rod * r = nullptr;
	Clock c;
	Integrator * integrator = nullptr;
	std::vector<RodEnergy*> energies;
	std::vector<RodConstraint*> constraints;
	MouseSpring * mouseSpring;

	real twist = 0.0;
	real rodTwist = 0.0;
	int numRodTwists = 0;

	//Interactive stuff;

	bool isMouseDown = false;
	Vec3e mousePosition;
	bool isRotate = false;
	bool isSetup = false;

	// Sound stuff
	#define SimulationLength 3.0 // in seconds
	const static std::size_t BufferSize = (std::size_t)(SampleRate * SimulationLength);
	double sampleBuffer[BufferSize];
	
	real tAtLastDraw = 0.0;
	bool stopNow = false;
	Vec3e ear2Pos = Vec3e(28.0, 10.0, 28.0);
	double sampleBuffer2[BufferSize];
	std::size_t curSample = 0;
	double sampleBuffer3[BufferSize];
	
	std::size_t multiSample = 21;
	FrameExporter fe;
};

#endif // ROD_SOUND_VIEWER_H
