/* 
* @Author: marinheiro
* @Date:   2014-10-08 20:25:13
* @Last Modified by:   Rafael Marinheiro
* @Last Modified time: 2014-10-30 23:37:05
*/

#include <iostream>
#include <fstream>

#include "glUtil.h"
#include "rodDraw.h"
#include "rodResource.h"

#include "SimulatorViewer.h"

#include <QKeyEvent>
#include <algorithm>

void SimulatorViewer::init(){
	this->camera()->setPosition(EtoQ(eyePos));
	this->camera()->setUpVector(qglviewer::Vec(0.0, 1.0, 0.0));
	this->camera()->lookAt(EtoQ(targetPos));
	this->camera()->setSceneCenter(EtoQ(targetPos));

	//Load the rod
	this->loadDefaultRod(42);

	#ifndef CONST_INTEGRATOR
		this->loadStdEnergies();
	#else
		this->loadStdEnergiesAndConsts();
	#endif //ifndef CONST_INTEGRATOR

	this->isSetup = true;
	camera()->setFlySpeed(1);
	setAnimationPeriod(16);
	// std::cout << animationPeriod() << std::endl;
	startAnimation();
}

void SimulatorViewer::keyPressEvent(QKeyEvent * event){
	if(event->key() == Qt::Key_R){
		if(event->modifiers() == Qt::ShiftModifier){
			this->twist = 0.0;
		} else{
			this->isRotate = !this->isRotate;
		}
	} else if(event->key() == Qt::Key_P){
		// std::cout << rod::resource::pathToResource("shader.glsl") << std::endl;
		this->running = !this->running;
		this->toggleAnimation();
	} else if(event->key() == Qt::Key_Y){
		this->running = false;
		this->loadRodFile("");
		this->loadStdEnergies();
	} else if(event->key() == Qt::Key_T){
		for(int i=1; i<r->numEdges(); i++) {
			std::cout << "edge " << i << " twist: " << (r->cur().rot(i) - r->cur().rot(i-1) + r->cur().refTwist(i))
			<< " (" << r->cur().rot(i) << " + " << r->cur().refTwist(i) <<
			" = " << r->cur().rot(i) + r->cur().refTwist(i) << ")\n";
		}
		std::cout << "numRodTwists: " << numRodTwists << "\n"; 
	} else{
		std::cout << this->camera()->flySpeed() << std::endl;
		QGLViewer::keyPressEvent(event);
	}
}

void SimulatorViewer::mousePressEvent(QMouseEvent * event){
	
	if(event->button() == Qt::LeftButton){
		if(!running) return;
		isMouseDown = true;
		qglviewer::Vec qorig;
		qglviewer::Vec qdirection;
		this->camera()->convertClickToLine(event->pos(), qorig, qdirection);
		Vec3e orig = QtoE(qorig);
		Vec3e direction = QtoE(qdirection);
		Vec3e plane = targetPos - orig;

		real rp = direction.dot(plane);
		if(rp != 0.0){
			real k = plane.dot(plane)/rp;
			mousePosition = orig + k * direction;
		} else{
			std::cerr << "Mouse ray did not intersect plane!\n";
		}
	} else{
		// // std::cout << "HIIII" << std::endl;
		// if(!r) return;
		// qglviewer::Vec qorig;
		// qglviewer::Vec qdirection;
		// this->camera()->convertClickToLine(event->pos(), qorig, qdirection);
		// Vec3e orig = QtoE(qorig);
		// Vec3e direction = QtoE(qdirection);

		// real tmin = INFINITY;
		// bool any = false;
		// for(int i = 0; i < r->numCPs(); i++){
		// 	bool hit = false;
		// 	real t = INFINITY;
		// 	std::cout << "HIIII" << std::endl;
		
		// 	{
		// 		Vec3e delta = r->cur().POS(i) - orig;
		// 		real k = delta.dot(direction);
		// 		real d = (delta - k*direction).norm();
		// 		real rr = constants::radius * 1.5;
				
		// 		if(d < rr){
		// 			real ddelta = sqrt(rr*rr - d*d);
		// 			if(d + ddelta > 0){
		// 				t = std::max(0.0, d-ddelta);
		// 				hit = true;
		// 			}
		// 		}
		// 	}

		// 	if(hit && t < tmin){
		// 		any = true;
		// 		tmin = t;
		// 	}
		// }

		// if(!any) return;
		// targetPos = orig + direction * tmin;

		// this->camera()->lookAt(EtoQ(targetPos));
		// this->camera()->setSceneCenter(EtoQ(targetPos));
		QGLViewer::mousePressEvent(event);
	}
}

void SimulatorViewer::mouseMoveEvent(QMouseEvent * event){
	// std::cout << "LEFT: " << (event->button() & Qt::LeftButton) << std::endl;
	// std::cout << "MID: " << (event->button() & Qt::MidButton) << std::endl;
	// std::cout << "RIGHT: " << (event->button() & Qt::RightButton) << std::endl;
	if(!running){
		QGLViewer::mouseMoveEvent(event);
	} else if(isMouseDown){
		qglviewer::Vec qorig;
		qglviewer::Vec qdirection;
		this->camera()->convertClickToLine(event->pos(), qorig, qdirection);
		Vec3e orig = QtoE(qorig);
		Vec3e direction = QtoE(qdirection);
		Vec3e plane = targetPos - orig;

		real rp = direction.dot(plane);
		if(rp != 0.0){
			real k = plane.dot(plane)/rp;
			mousePosition = orig + k * direction;
		} else{
			std::cerr << "Mouse ray did not intersect plane!\n";
		}
	} else{
		QGLViewer::mouseMoveEvent(event);	
	}
}

void SimulatorViewer::mouseReleaseEvent(QMouseEvent * event){
	if(!running) return;
	isMouseDown = false;
}

void SimulatorViewer::animate(){
	if (!running || !isSetup) return;

	Vec3e mp;
	if (isMouseDown) mp << mousePosition.x(), mousePosition.y(), mousePosition.z();
	mouseSpring->setMouse(mp, isMouseDown);
	
#ifndef CONST_INTEGRATOR
	while (!integrator->integrate(c)) {
		if (c.canDecreaseTimestep()) {
			c.suggestTimestep(c.timestep() / 2.0);
			std::cout << "Decreasing timestep: " << c.timestep() << "\n";
		} else {
			std::cout << "Simulation Failed!\n";
			running = false;
			return;
		}
	}
#else
	while (!cIntegrator->integrate(c)) {
		std::cout << "wat\n";
		throw;
	}
#endif // ifdef CONST_INTEGRATOR
	
	/// Update Bishop frame
	r->next().updateReferenceFrames(r->cur());
	
#ifndef CONST_INTEGRATOR
	/// Update material frame rotation
	if (isRotate) {
		twist += 2.0*constants::pi*c.timestep();
	}
	
	Vec3e uRef = parallelTransport(r->next().edge(0), r->next().edge(r->numEdges()-1), r->next().u[0]);
	real cosTwist = r->next().u[r->numEdges()-1].dot(uRef.normalized());
	real oldTwist = rodTwist;
	if (cosTwist >= 1.0) { // Avoid values like 1.0000000012 that introduce NaNs
		rodTwist = 0.0;
	} else if (cosTwist <= -1.0) {
		rodTwist = constants::pi;
	} else {
		rodTwist = acos(cosTwist);
	}
	// Flip the sign if necessary
	if (r->next().v(r->numEdges()-1).dot(uRef) > 0.0) {
		rodTwist = -rodTwist;
	}
	real diff = rodTwist - oldTwist;
	if (diff < -constants::pi) {
		numRodTwists += 1;
	} else if (diff > constants::pi) {
		numRodTwists -= 1;
	}
	r->next().rot(r->numEdges()-1) = twist - rodTwist;
	if (!static_cast<IMEXIntegrator*>(integrator)->setRotations()) {
		std::cout << "rotations failed";
	}
#endif // ifndef CONST_INTEGRATOR
	
	// Swap Rods
	r->swapRods();

	c.increment();	
}

void SimulatorViewer::draw(){
	if (!isSetup) {
		return;
	}

	// Clear out the window with grey
	rod::gl::clear(0.45, 0.45, 0.5);
	
	// Enable alpha blending and depth testing
	// gl::enableAlphaBlending();
	// gl::enableDepthRead(true);
	// gl::enableDepthWrite(true);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		
	// Draw framerate counter
	// gl::setMatricesWindow(getWindowSize());
	// std::stringstream ss;
	// ss << getAverageFps();
	// gl::drawStringRight(ss.str(),
	// 					Vec2c(getWindowWidth()-toPixels(10), getWindowHeight()-toPixels(20)),
	// 					Color(0.0, 0.0, 0.0),
	// 					Font("Arial", toPixels(12)));
	
	// // Set projection/modelview matrices
	// gl::setMatrices(cam);
	
	// Draw the rod and the normal of the bishop frame
	for(int i=0; i<r->numEdges(); i++) {
		Vec3e p0 = r->cur().POS(i);
		Vec3e p1 = r->cur().POS(i+1);
		rod::gl::drawLine(p0, p1);
		rod::gl::color(1.0, 1.0, 0.0);
		rod::gl::lineWidth(1.0);
		Vec3e u = r->cur().u[i];
		rod::gl::drawLine((p0+p1)/2.0, (p0+p1)/2.0+u*(p1-p0).norm()*2.0);
	}
	
	// m.apply();
	
	// l->setDiffuse(Color::white());
	// l->setAmbient(Color::white());
	// l->setPosition(Vec3c(0.0, 50.0, 0.0));
	// l->enable();
	
	// diffuseProg.bind();
	// for (int i=0; i<r->numCPs(); i++) {
	// 	gl::pushModelView();
	// 	gl::translate(EtoC(r->cur().POS(i)));
	// 	spheredl->draw();
	// 	gl::popModelView();
	// }
	// diffuseProg.unbind();
	
	// rodProg.bind();

	// floorTex.enableAndBind();
	// gl::draw(floor);
	// floorTex.disable();
	
	// rodProg.unbind();
	
	// // Draw rod edges
	// rodProg.bind();
	// rodTex.enableAndBind();
	// for (int i=0; i<r->numEdges(); i++) {
	// 	gl::pushModelView();
	// 	Vec3c v = EtoC(r->cur().edge(i).normalized());
		
	// 	gl::translate(EtoC(r->cur().POS(i)));
	// 	Quaternion<real> q(Vec3c(0.0, 1.0, 0.0), v);
	// 	real angle = acos(std::max((real)-1.0, std::min((real)1.0, (q*Vec3c(-1.0, 0.0, 0.0)).dot(EtoC(r->cur().u[i])))));
	// 	if ((q*Vec3c(-1.0, 0.0, 0.0)).dot(EtoC(r->cur().v(i))) > 0.0) angle = -angle;
	// 	gl::rotate(Quaternion<real>(v, angle));
	// 	gl::rotate(q);
	// 	gl::rotate(Vec3c(0.0, r->cur().rot(i)*180.0/constants::pi, 0.0));
	// 	gl::scale(1.0, r->cur().edgeLength(i), 1.0);
	// 	cylinderdl->draw();
	// 	gl::popModelView();
	// }
	// rodTex.unbind();
	// rodProg.unbind();
	
	for (RodEnergy* e : energies) {
		e->draw(c.timestep());
	}
#ifndef CONST_INTEGRATOR
	integrator->draw();
#endif // ifndef CONST_INTEGRATOR

}

void SimulatorViewer::loadRodFile(std::string filename) {
	// if (filename.empty()) filename = getOpenFilePath().string();
	if (filename.empty()) return;
	
	std::ifstream rodFile(filename);
	
	if (!rodFile.is_open()) {
		std::cerr << filename << " failed to open!\n";
		return;
	}
	
	std::string line;
	std::getline(rodFile, line);
	const std::size_t numPoints = std::stoi(line);
	
	VecXe rodPos(3*numPoints);
	
	for (int i=0; i<3*numPoints; i++) {
		std::string line;
		std::getline(rodFile, line);
		if(!line.empty()) {
			rodPos(i) = std::stof(line);
		}
	}
	Vec3e u;
	for (int i=0; i<3; i++) {
		std::string line;
		std::getline(rodFile, line);
		u(i) = std::stof(line);
	}
	assert((rodPos.segment<3>(3) - rodPos.segment<3>(0)).dot(u) < 5.0e-6);
	
	rodFile.close();
	if (r) delete r;
	r = new Rod(rodPos, u);
}

void SimulatorViewer::loadDefaultRod(int numPoints) {
	if (r) delete r;
	
	eyePos = Vec3e(40.0, 10.0, 0.0);
	targetPos = Vec3e(0.0, 10.0, 0.0);
	
	this->camera()->setPosition(EtoQ(eyePos));
	this->camera()->setUpVector(qglviewer::Vec(0.0, 1.0, 0.0));
	this->camera()->lookAt(EtoQ(targetPos));
	this->camera()->setSceneCenter(EtoQ(targetPos));
	
	Vec3e start = Vec3e(0.0, 1.6069, 0.0);
	Vec3e end   = Vec3e(0.0, 1.0, 0.0);
	
	Vec3e u     = (end-start).cross(Vec3e(0.0, 0.1, 0.0)).normalized();
	if (u.hasNaN() || u.norm() < 0.95) {
		u << 1.0, 0.0, 0.0;
	}
	
	VecXe rodPos(3*numPoints);
	for(int i=0; i < numPoints; i++) {
		real t = ((real) i) / (real) (numPoints -1);
		rodPos.segment<3>(3*i) = (1-t)*start + t*end;
	}
	
	real massPerPoint = 0.1;
	VecXe mass = VecXe::Constant(numPoints, massPerPoint);
	
	eyePos = Vec3e(5.0, 1.5, 0.0);
	targetPos = Vec3e(0.0, 1.5, 0.0);
	
	this->camera()->setPosition(EtoQ(eyePos));
	this->camera()->setUpVector(qglviewer::Vec(0.0, 1.0, 0.0));
	this->camera()->lookAt(EtoQ(targetPos));
	this->camera()->setSceneCenter(EtoQ(targetPos));
	
	r = new Rod(rodPos, Vec3e(0.0, 0.0, 1.0), &mass, 1e7, 80.0);
}

void SimulatorViewer::loadStdEnergies() {
	// Create Rod Energies - Add in the order they are most likely to fail during evaluation
	assert(r && "Tried to load energies on a null rod");
	for (RodEnergy* e : energies) {
		delete e;
	}
	energies.clear();
	
	
	RodEnergy* stretch = new Stretching(*r, Implicit);
	energies.push_back(stretch);
	
	RodEnergy* bending = new Bending(*r, Implicit);
	energies.push_back(bending);
	
	RodEnergy* twisting = new Twisting(*r, Explicit);
	energies.push_back(twisting);
	
	Vec3e dir(0.0, -9.8, 0.0);
	RodEnergy* gravity = new Gravity(*r, Explicit, dir);
	energies.push_back(gravity);
	
	mouseSpring = new MouseSpring(*r, Explicit, r->numCPs()-1, 100.0);
	energies.push_back(mouseSpring);
	
	RodEnergy* intContact = new IntContact(*r, Implicit);
	energies.push_back(intContact);
	
	Spring* clamp1 = new Spring(*r, Implicit, 0, 500.0);
	Vec3e pos0 = r->rest().POS(0);
	clamp1->setClamp(pos0);
	Spring* clamp2 = new Spring(*r, Implicit, 1, 1000.0);
	Vec3e pos1 = r->rest().POS(1);
	clamp2->setClamp(pos1);
	energies.push_back(clamp1);
	energies.push_back(clamp2);
	
	if (integrator) delete integrator;
	integrator = new IMEXIntegrator(energies, *r);
}

void SimulatorViewer::loadStdEnergiesAndConsts() {
	assert(r && "Tried to load energies and constraints on a null rod");
	for (RodEnergy* e : energies) {
		delete e;
	}
	for (RodConstraint* c : constraints) {
		delete c;
	}
	energies.clear();
	constraints.clear();
	
	Vec3e dir(0.0, -9.8, 0.0);
	RodEnergy* gravity = new Gravity(*r, Explicit, dir);
	energies.push_back(gravity);
	
	Spring* clamp1 = new Spring(*r, Implicit, 0, 500.0);
	Vec3e pos0 = r->rest().POS(0);
	clamp1->setClamp(pos0);
	Spring* clamp2 = new Spring(*r, Implicit, 1, 1000.0);
	Vec3e pos1 = r->rest().POS(1);
	clamp2->setClamp(pos1);
	energies.push_back(clamp1);
	energies.push_back(clamp2);
	
	mouseSpring = new MouseSpring(*r, Explicit, r->numCPs()-1, 10.0);
	energies.push_back(mouseSpring);
	
	RodConstraint* length = new Length(*r);
	constraints.push_back(length);
	
	if (cIntegrator) delete cIntegrator;
	cIntegrator = new ConstraintIntegrator(*r, energies, constraints);
}