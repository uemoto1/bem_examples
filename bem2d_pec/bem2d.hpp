//      bem2d.hpp
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
#ifndef __bem2d_hpp__
#define __bem2d_hpp__

#include "def.hpp"
#include "bessel.hpp"

class CIncident {
public:
  vec2d ivec;
  dcmplx a;
  double k;
  double w;

  CIncident();
  CIncident(const vec2d &, double, const dcmplx &);
  dcmplx field(const vec2d &);
};

class CObject {
public:
  vector<vec2d> node;
  Eigen::VectorXcd data;

  CObject();
  double length(int);
  bool inside(const vec2d &);
};

dcmplx h0(double);
dcmplx green(double, const vec2d &);
void bem2d(CObject &, CIncident &);
dcmplx field(const vec2d &, CObject &, CIncident &);

#endif
