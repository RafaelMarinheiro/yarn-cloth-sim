/* 
* @Author: marinheiro
* @Date:   2014-10-08 21:14:51
* @Last Modified by:   marinheiro
* @Last Modified time: 2014-10-08 21:52:42
*/

#include "SimulatorViewer.h"
#include <qapplication.h>

int main(int argc, char** argv)
{
	// Read command lines arguments.
	QApplication application(argc, argv);
	
	// Instantiate the viewer.
	SimulatorViewer viewer;

	viewer.setWindowTitle("Simulator");

	// Make the viewer window visible on screen.
	viewer.show();

	// Run main loop.
	return application.exec();
}