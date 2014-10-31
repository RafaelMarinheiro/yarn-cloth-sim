#include <iostream>
#include <fstream>

#include "glUtil.h"
#include "rodDraw.h"
#include "rodResource.h"

#include "RodSoundViewer.h"

#include <QKeyEvent>

void RodSoundViewer::init(){
	this->camera()->setPosition(EtoQ(eyePos));
	this->camera()->setUpVector(qglviewer::Vec(0.0, 1.0, 0.0));
	this->camera()->lookAt(EtoQ(targetPos));
	this->camera()->setSceneCenter(EtoQ(targetPos));

	//Load the rod
	this->loadDefaultRod(50);

	this->loadStdEnergies();

	this->fe.viewer = this;
	this->setSnapshotFormat("PNG");
	PROFILER_START("Total");

	this->isSetup = true;
	camera()->setFlySpeed(1);
	setAnimationPeriod(16);
	// std::cout << animationPeriod() << std::endl;
	startAnimation();	
}

void RodSoundViewer::mousePressEvent(QMouseEvent * event){
	
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

void RodSoundViewer::mouseMoveEvent(QMouseEvent * event){
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

void RodSoundViewer::mouseReleaseEvent(QMouseEvent * event){
	if(!running) return;
	isMouseDown = false;
}

void RodSoundViewer::keyPressEvent(QKeyEvent * event){
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
	} else if(event->key() == Qt::Key_S){
		this->stopNow = true;
	} else{
		std::cout << this->camera()->flySpeed() << std::endl;
		QGLViewer::keyPressEvent(event);
	}
}

void RodSoundViewer::update(){
	if (!running || !isSetup) return;

	// std::cout << curSample << std::endl;

	if (curSample % 5000 == 0 && curSample != 0 && c.getTicks() % multiSample == 0) {
		std::cout << curSample << " / " << BufferSize << " (" << (curSample*100.0)/BufferSize << "%)\n";
		PROFILER_PRINT_ELAPSED();
		PROFILER_RESET_ALL();
		std::cout << "\n";
	}

	if (curSample >= BufferSize || stopNow) { // We're done!
		sampleBuffer[0] = 0.0; // To prevent the click of forces suddenly being applied
		double max = 0;
		for (int i=0; i<BufferSize; i++) {
			max = std::max(max, std::fabs(sampleBuffer[i]));
		}
		std::cout << "Max1: " << max << "\n";
		uint16_t buffer[BufferSize];
		for (int i=0; i<BufferSize; i++) {
			buffer[i] = toSample(sampleBuffer[i], max);
		}
		writeWAVData(rod::resource::pathToOutput("result.wav").data(), buffer,
								 curSample * sizeof(uint16_t), SampleRate, 1);
		
		sampleBuffer2[0] = 0.0;
		max = 0;
		for (int i=0; i<BufferSize; i++) {
			max = std::max(max, std::fabs(sampleBuffer2[i]));
		}
		std::cout << "Max2: " << max << "\n";
		for (int i=0; i<BufferSize; i++) {
			buffer[i] = toSample(sampleBuffer2[i], max);
		}
		writeWAVData(rod::resource::pathToOutput("result2.wav").data(), buffer,
								 curSample * sizeof(uint16_t), SampleRate, 1);
		
		sampleBuffer3[0] = 0.0;
		max = 0;
		for (int i=0; i<BufferSize; i++) {
			max = std::max(max, std::fabs(sampleBuffer3[i]));
		}
		std::cout << "Max3: " << max << "\n";
		for (int i=0; i<BufferSize; i++) {
			buffer[i] = toSample(sampleBuffer3[i], max);
		}
		writeWAVData(rod::resource::pathToOutput("result3.wav").data(), buffer,
								 curSample * sizeof(uint16_t), SampleRate, 1);
		
		fe.writeMPEG("result");
		// std::cout << "Total simulation time: " << app::getElapsedSeconds() << "\n"; // FIXME: This is inaccurate
		
		running = false;
		return;
	}

	PROFILER_START("Update");

	c.suggestTimestep(1.0 / (real) SampleRate / multiSample);
	// std::cout << "Timestep: " << c.timestep() << std::endl;
	// FIXME: Normally the frame exporter would suggest a timestep, but this interferes with the audio
	// recording, as it assumes all timesteps are 1/SampleRate. However, any error the frame exporter
	// experiences is small since 1/60 >> 1/SampleRate.
	// fe.suggestTimestep(c);
  

	Vec3e mp;
	if (isMouseDown) mp << mousePosition.x(), mousePosition.y(), mousePosition.z();
	mouseSpring->setMouse(mp, isMouseDown);
	
	if(!integrator->integrate(c)) throw;
	
	/// Update Bishop frame
	r->next().updateReferenceFrames(r->cur());

	// Sound Calculations
	if (c.getTicks() % multiSample == 0) {
		real sample = 0;
		real sample2 = 0;
		real avgX = 0;
		VecXe jerkVec = r->next().dVel - r->cur().dVel;
		for (int i=1; i<r->numCPs()-1; i++) {
			avgX += r->next().VEL(i).x();
			
			// Calculate jerk
			Vec3e jerk = jerkVec.segment<3>(3*i);
			// Project jerk to transverse plane
			Vec3e tPlaneNormal = (r->next().edge(i-1) + r->next().edge(i)).normalized();
			jerk = jerk - jerk.dot(tPlaneNormal) * tPlaneNormal; // Vector rejection of jerk from tPlaneNormal
			
			/*
			real m0 = r->restVoronoiLength(i)*constants::pi*r->radius()*r->radius()*constants::rhoAir;
			// Rotation to align system so that the cylinder is coaxial with the z-axis
			Eigen::Quaternion<real> q = Eigen::Quaternion<real>::FromTwoVectors(tPlaneNormal, Vec3e(0, 0, 1));
			Vec3e rotJerk = q * jerk;
			rotJerk = rotJerk.cwiseProduct(Vec3e(2.0*m0, 2.0*m0, m0));
			
			// Calculate sample contribution
			Vec3e earVec = CtoE(eyePos) - r->next().points[i].pos;
			sample +=  (q * earVec).dot(rotJerk) / (4.0 * constants::pi * constants::cAir * earVec.dot(earVec));
			
			earVec = ear2Pos - r->next().points[i].pos;
			sample2 +=  (q * earVec).dot(rotJerk) / (4.0 * constants::pi * constants::cAir * earVec.dot(earVec));
			*/
			 
			
			Vec3e earVec = eyePos - r->next().POS(i);
			// Calculate sample contribution
			sample += r->getCS()[i].calcSample(earVec, jerk);
		
			earVec = ear2Pos - r->next().POS(i);
			sample2 += r->getCS()[i].calcSample(earVec, jerk);
		}
		avgX = avgX/(r->numCPs()-2);
		sampleBuffer[curSample] = sample;
		sampleBuffer2[curSample] = sample2;
		
		sampleBuffer3[curSample] = r->next().VEL(r->numCPs()/2).x() - avgX;
		
		curSample++;
	}
	
	// Swap Rods
	r->swapRods();

	c.increment();
	PROFILER_STOP("Update");	
}

void RodSoundViewer::animate(){
	// std::cout << "Animarei" << std::endl;
	// std::cout << running << std::endl;
	// std::cout << fe.nextTimestep(c) << std::endl;
	// std::cout << 1.0 / (real) SampleRate << std::endl;
	while(running && 
			// app::getElapsedSeconds() - tAtLastDraw < 1.0/app::getFrameRate() &&
			fe.nextTimestep(c) > 0.0){
		this->update();
	}
	// std::cout << "Animando: " << fe.nextTimestep(c) << " " << c.time() << std::endl;
	fe.record(c);
	// std::cout << "Animei" << std::endl;
}

void RodSoundViewer::draw(){
	if (!isSetup) {
		return;
	}

	// tAtLastDraw = app::getElapsedSeconds();

	PROFILER_START("DRAW");
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

	integrator->draw();

	PROFILER_STOP("Draw");
}

void RodSoundViewer::loadRodFile(std::string filename) {
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

void RodSoundViewer::loadDefaultRod(int numPoints) {
	if (r) delete r;
	
	eyePos = Vec3e(5.0, 1.5, 0.0);
	targetPos = Vec3e(0.0, 1.5, 0.0);
	
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
	
	real totalMass = constants::radius * constants::radius * constants::pi * (start - end).norm() * constants::rhoRod;
	real massPerPoint = totalMass / numPoints;
	VecXe mass = VecXe::Constant(numPoints, massPerPoint);
	
	r = new Rod(rodPos, Vec3e(0.0, 0.0, 1.0), &mass, 1e7, 80.0);

	real l = (start - end).norm();
	real kappa = sqrt(r->youngsModulus * r->getCS()[0].areaMoment()(0, 0) /
                    (constants::rhoRod * r->getCS()[0].area() * l * l * l * l));
	std::cout << "kappa: " << kappa << "\n";
	real h = l / (numPoints-1);
	real k = 1.0 / (44100*multiSample);
	real mu = kappa * k / (h * h);
	std::cout << "mu: " << mu << "\n";
	
	real fmax = asin(2.0 * mu) / (constants::pi * k);
	std::cout << "fmax: " << fmax << "\n";
}

void RodSoundViewer::loadStdEnergies() {
	// Create Rod Energies - Add in the order they are most likely to fail during evaluation
	assert(r && "Tried to load energies on a null rod");
	for (RodEnergy* e : energies) {
		delete e;
	}
	energies.clear();
	
	RodEnergy* stretch = new Stretching(*r, Explicit);
//  energies.push_back(stretch);
	// OR
	//RodConstraint* length = new Length(*r);
	//constraints.push_back(length);
	
	RodEnergy* bending = new Bending(*r, Explicit);
	energies.push_back(bending);
	
	RodEnergy* fembending = new FEMBending(*r, Explicit);
//  energies.push_back(fembending);
	
	RodEnergy* twisting = new Twisting(*r, Explicit);
//  energies.push_back(twisting);
	
	Vec3e dir(0.0, -9.8, 0.0);
	RodEnergy* gravity = new Gravity(*r, Explicit, dir);
//  energies.push_back(gravity);
	
	mouseSpring = new MouseSpring(*r, Explicit, r->numCPs()-1, 100.0);
	energies.push_back(mouseSpring);
	
	Vec3e floorNormal(0.0, 1.0, 0.0);
	Vec3e floorOrigin = Vec3e::Zero();
	RodEnergy* floor = new PlaneContact(*r, Explicit, floorNormal, floorOrigin, 5000.0);
//  energies.push_back(floor);
	
	Vec3e imp1dir(1.0e-10, 0.0, 0.0);
	Vec3e imp2dir(-1.0e-5, 0.0, 0.0);
	Vec3e imp3dir(0.0, 0.0, 1.0e-10);
	RodEnergy* imp1 = new Impulse(*r, Explicit, c, 0.2, 0.201, imp1dir, 3);
	RodEnergy* imp2 = new Impulse(*r, Explicit, c, 0.2, 0.201, imp2dir, r->numCPs()-2);
	RodEnergy* imp3 = new Impulse(*r, Explicit, c, 0.2, 0.201, imp3dir, r->numCPs()-1);
	energies.push_back(imp1);  // energies.push_back(imp2); // energies.push_back(imp3);
	
	RodEnergy* spr1 = new Spring(*r, Explicit, 0, 1e5);
	Vec3e spr1clamp = r->rest().POS(0) + Vec3e(0.0, 0.1, 0.0);
	static_cast<Spring*>(spr1)->setClamp(spr1clamp);
	RodEnergy* spr2 = new Spring(*r, Explicit, r->numCPs()-1, 1e5);
	Vec3e spr2clamp = r->rest().POS(r->numCPs()-1);
	static_cast<Spring*>(spr2)->setClamp(spr2clamp);
//  energies.push_back(spr1); energies.push_back(spr2);
	
	/*
	
	RodEnergy* intContact = new IntContact(*r, Explicit);
	energies.push_back(intContact);
	
	 */
	
	if (integrator) delete integrator;
	integrator = new ExIntegrator(*r, energies, &constraints);
	// integrator = new ExpIntegrator(*r, energies);
}