#ifndef GL_UTIL_H
#define GL_UTIL_H

#ifdef __APPLE__
	#include <OpenGL/gl.h>
	#include <OpenGL/glu.h>
#else
	#ifdef _WIN32
	  #include <windows.h>
	#endif
	
	#include <GL/gl.h>
	#include <GL/glu.h>
#endif

#endif // GL_UTIL_H
