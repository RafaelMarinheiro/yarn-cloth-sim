//
//  Shader.h
//  SPH
//
//  Created by Rafael Farias Marinheiro on 2/18/14.
//  Copyright (c) 2014 Rafael Farias Marinheiro. All rights reserved.
//

#ifndef __SPH__Shader__
#define __SPH__Shader__

#include "glUtil.h"
#include <map>
#include <string>

namespace rod{
    namespace gl{
        using namespace std;
        
        class Shader
        {
        public:
            Shader(void);
            ~Shader(void);
            void loadFromString(GLenum whichShader, const string& source);
            void loadFromFile(GLenum whichShader, const string& filename);
            void createAndLinkProgram();
            void use();
            void unUse();
            void addAttribute(const string& attribute);
            void addUniform(const string& uniform);
            GLuint getProgram() const;
            //An indexer that returns the location of the attribute/uniform
            GLuint operator[](const string& attribute);
            GLuint operator()(const string& uniform);
            //Program deletion
            void deleteProgram() {glDeleteProgram(_program);_program=-1;}
        private:
            enum ShaderType {VERTEX_SHADER, FRAGMENT_SHADER, GEOMETRY_SHADER};
            GLuint	_program;
            int _totalShaders;
            GLuint _shaders[3];//0->vertexshader, 1->fragmentshader, 2->geometryshader
            map<string,GLuint> _attributeList;
            map<string,GLuint> _uniformLocationList;
        };	
        
    }
}

#endif /* defined(__SPH__Shader__) */
