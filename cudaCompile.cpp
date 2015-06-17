#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cstdlib>
#include <sstream>

int main(int argc, char * argv []){

std::string codeFile;


if (argc != 2){
std::cout << "USAGE: cudaCompile <YourCodeFileHere (do not add .cu extension!)> " << std::endl;
std::cout << "Example: I want to compile cosSphere.cu so I run: cudaCompile cosSphere" << std::endl;
}

else{
	codeFile  = std::string( argv[1] ); //name of text file with listed image filenames
	std::string cmd = " \" \"C:\\Program Files (x86)\\Microsoft Visual Studio 10.0\\VC\\bin\\vcvars32.bat\" && \"C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA\\v5.5\\bin\\nvcc.exe\" --use-local-env --cl-version 2010 -ccbin \"C:\\Program Files (x86)\\Microsoft Visual Studio 10.0\\VC\\bin\" --machine 32 -arch sm_20 -I \"C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA\\v5.5\\include\" -I .\\acquisition -I .\\functors -I .\\kernels -I .\\params -I .\\primitives -I .\\substrate -I .\\util -I .\\blochdiff " ;
	cmd += "-o " + codeFile + ".exe " + codeFile + ".cu \""; 
	// std::cout << cmd << std::endl;
	system(cmd.c_str());
}

}