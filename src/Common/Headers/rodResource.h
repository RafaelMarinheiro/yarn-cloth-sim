#ifndef ROD_RESOURCE_H
#define ROD_RESOURCE_H

#include <string>

	#define ROD_RESOURCE( LOCALPREFIX, PATH ) #PATH

	namespace rod{
		namespace resource{
			std::string pathToResource(const std::string & resource);
			std::string pathToOutput(const std::string & output);
		}
	}

#endif // ROD_RESOURCE_H
