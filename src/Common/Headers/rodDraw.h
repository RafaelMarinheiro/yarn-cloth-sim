#ifndef ROD_DRAW_H
#define ROD_DRAW_H

#include "Defines.h"

namespace rod{
	namespace gl{
		void clear(float r, float g, float b);
		void clear(double r, double g, double b);

		void color(float r, float g, float b);
		void color(double r, double g, double b);

		void lineWidth(float w);
		void lineWidth(double w);

		void drawLineFloat(const Eigen::Matrix<float, 3, 1> & a, const Eigen::Matrix<float, 3, 1> & b);
		void drawLineDouble(const Eigen::Matrix<double, 3, 1> & a, const Eigen::Matrix<double, 3, 1> & b);
		void drawLine(const Vec3e & a, const Vec3e & b);
		
	}
}

#endif // ROD_DRAW_H
