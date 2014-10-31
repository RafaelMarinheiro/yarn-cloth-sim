//
//  FrameExporter.h
//  Visualizer
//
//  Created by eschweickart on 6/12/14.
//
//

#ifndef Visualizer_FrameExporter_h
#define Visualizer_FrameExporter_h

#include "Defines.h"

#include "Constants.h"
#include "Clock.h"

#include <QGLViewer/qglviewer.h>
#include "rodResource.h"

#include <iostream>

class FrameExporter {
  
  real fr;
  real lastFrame;
  std::size_t frameCount = 0;
public:
  FrameExporter(real framerate = 1.0/60.0) : fr(framerate), lastFrame(-fr) { }
  QGLViewer * viewer;
  void const inline suggestTimestep(Clock& c) const {
    // NOTE: this may be negative, but clock accounts for this.
    c.suggestTimestep(nextTimestep(c));
  }
  
  real const inline nextTimestep(const Clock& c) const {
    return fr - (c.time() - lastFrame);
  }
  
  void record(const Clock& c) {
    if (c.time() - lastFrame < fr){
      // std::cout << "WHYYY" << std::endl;
      return;
    }
    if (c.time() - lastFrame > fr * 1.01) {
      std::cout << "Warning: jumpy frames in Framewriter\n";
    }
    
    // std::cout << "Frame count: " << frameCount << std::endl;
    frameCount++;
    lastFrame = c.time();
    // std::cout << "lastFrame: " << lastFrame << std::endl;
    // std::cout << "c.time(): " << c.time() << std::endl;
    std::stringstream path;
    path << "png/" << frameCount << ".png";
    // path << constants::ResultPath << "png/" << frameCount << ".png";
    viewer->saveSnapshot(QString(rod::resource::pathToOutput(path.str()).data()), true);
    // writeImage(path.str(), ci::app::copyWindowSurface());
  }
  
  void writeMPEG(char const* fname) {
    std::string mm = rod::resource::pathToResource("MovieMaker.sh") + " " + rod::resource::pathToOutput("") + " " + fname;
    system(mm.data());
  }
};

#endif
