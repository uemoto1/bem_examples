//      def.hpp
//
//      Copyright 2014 uemoto <uemoto@yoshida-lab>
//
//      This program is free software; you can redistribute it and/or modify
//      it under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//
//      This program is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//
//      You should have received a copy of the GNU General Public License
//      along with this program; if not, write to the Free Software
//      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//      MA 02110-1301, USA.
//
//
#ifndef __def_hpp__
#define __def_hpp__

#include <iostream>
#include <cstdio>
#include <string>
#include <sstream>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/IterativeLinearSolvers>

// Principal constants
const double pi = 3.14159265358979;
const double eular = 0.577215664901;
const double hbar = 6.58211889e-16;  // eV/s
const double cspeed = 2.99792458e+8; // m/s

// Aliases for convenience.
typedef std::complex<double> dcmplx;
typedef Eigen::Vector2d vec2d;
typedef Eigen::Matrix2d mat2d;
typedef Eigen::Vector2cd vec2z;
typedef Eigen::Matrix2cd mat2z;
typedef Eigen::MatrixXcd matx;
typedef Eigen::VectorXcd vecx;

using std::cout;
using std::cin;
using std::endl;
using std::printf;
using std::stod;
using std::log;
using std::sqrt;
using std::sin;
using std::cos;
using std::vector;
using std::string;
using std::getline;
using std::stringstream;

#endif
