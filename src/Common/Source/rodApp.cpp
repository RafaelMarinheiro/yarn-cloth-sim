#include "rodApp.h"

#include <QCoreApplication>

#include <boost/filesystem.hpp>

#include <iostream>

namespace rod{
	namespace app{
		std::string getAppExecutablePath(){
			return QCoreApplication::applicationFilePath().toStdString();
		}

		std::string getAppBinPath(){
			boost::filesystem::path p(QCoreApplication::applicationFilePath().toStdString());
			while(p.filename().string() != "bin"){
				p = p.parent_path();
			}

			// std::cout << p << std::endl;
			return p.string();
		}

		std::string getAppResourcesPath(){
			boost::filesystem::path p(getAppBinPath());
			p = p.parent_path().parent_path();
			p += "/resources/";

			// std::cout << p << std::endl;
			return p.string();
		}

		std::string getAppOutputPath(){
			boost::filesystem::path p(getAppBinPath());
			p = p.parent_path().parent_path();
			p += "/result/";
			// std::cout << p << std::endl;
			return p.string();
		}
	}
}