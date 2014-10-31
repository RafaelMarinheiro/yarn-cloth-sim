#ifndef ROD_APP_H
#define ROD_APP_H

#include <string>

namespace rod{
	namespace app{
		std::string getAppBinPath();
		std::string getAppResourcesPath();
		std::string getAppOutputPath();
	}
}

#endif // ROD_APP_H
