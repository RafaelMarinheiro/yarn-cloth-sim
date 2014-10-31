/* 
* @Author: marinheiro
* @Date:   2014-10-08 19:48:51
* @Last Modified by:   marinheiro
* @Last Modified time: 2014-10-08 19:52:24
*/

#include "rodResource.h"
#include "rodApp.h"

namespace rod{
	namespace resource{
		std::string pathToResource(const std::string & resource){
			return rod::app::getAppResourcesPath() + resource;
		}

		std::string pathToOutput(const std::string & output){
			return rod::app::getAppOutputPath() + output;
		}
	}
}