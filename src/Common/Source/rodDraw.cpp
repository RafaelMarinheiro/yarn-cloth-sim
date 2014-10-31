/* 
* @Author: marinheiro
* @Date:   2014-10-08 18:29:16
* @Last Modified by:   marinheiro
* @Last Modified time: 2014-10-08 22:11:46
*/

#include "rodDraw.h"

#include "glUtil.h"

namespace rod{
	namespace gl{

		void clear(float r, float g, float b){
			glClearColor(r, g, b, 1.0);
		}

		void clear(double r, double g, double b){
			glClearColor(r, g, b, 1.0);
		}

		void color(float r, float g, float b){
			glColor3f(r, g, b);
		}

		void color(double r, double g, double b){
			glColor3d(r, g, b);
		}

		void lineWidth(float w){
			glLineWidth(w);
		}

		void lineWidth(double w){
			glLineWidth(w);
		}

		void drawLineFloat(const Eigen::Matrix<float, 3, 1> & a, const Eigen::Matrix<float, 3, 1> & b){
			glBegin(GL_LINES);
				glVertex3f(a.x(), a.y(), a.z());
				glVertex3f(a.x(), b.y(), b.z());
			glEnd();
		}

		void drawLineDouble(const Eigen::Matrix<double, 3, 1> & a, const Eigen::Matrix<double, 3, 1> & b){
			glBegin(GL_LINES);
				glVertex3d(a.x(), a.y(), a.z());
				glVertex3d(a.x(), b.y(), b.z());
			glEnd();
		}

		void drawLine(const Vec3e & a, const Vec3e & b){ drawLineDouble(a, b); }
	}
}