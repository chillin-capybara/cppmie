//
// Created by Marcell Pigniczki on 24.10.21.
//

#include <iostream>
#include "cppmie/cppmie.h"

#define CONST_PI (3.14159265358979323846)

int main(){
	std::cout << "Hello World!" << std::endl;
	double diameter = 1.0;  // Diameter in microns
	double wavelength = 0.6328; // Wavelength in microns
	double x = CONST_PI * diameter / wavelength;
	std::complex<double> m{1.5, 0.0};

	cppmie::mie(x, m.real());

	return 0;
}