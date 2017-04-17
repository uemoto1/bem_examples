//      bessel.cpp
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

#include "bessel.hpp"

double besj0(double z) {
  double j0;
  if (abs(z) < bessel_lim) {
    double d = z * 0.5;
    double a = 1.0;
    double s = a;
    for (int k = 1; k < bessel_order; k++) {
      a = -a * (d / k) * (d / k);
      s += a;
    }
    j0 = s;
  } else {
    j0 = sqrt(2.0 / (pi * z)) * cos(z - 0.25 * pi);
  }
  return j0;
}

double besj1(double z) {
  double j1;
  if (abs(z) < bessel_lim) {
    double d = 0.5 * z;
    double a = d;
    double s = a;
    for (int k = 1; k <= bessel_order; k++) {
      a = -a * (d * d) / (k * (k + 1));
      s += a;
    }
    j1 = s;
  } else {
    j1 = sqrt(2.0 / (pi * z)) * cos(z - 0.75 * pi);
  }
  return j1;
}

double besy0(double z) {
  double y0;
  if (abs(z) < bessel_lim) {
    double d = 0.5 * z;
    double a = 2.0 * besj0(z) * (eular + log(d));
    double b = 1.0;
    double c = 0.0;
    double s = a;
    for (int k = 1; k < bessel_order; k++) {
      b = -b * (d / k) * (d / k);
      c += (2.0 / k);
      s -= b * c;
    }
    y0 = s / pi;
  } else {
    y0 = sqrt(2.0 / (pi * z)) * sin(z - 0.25 * pi);
  }
  return y0;
}

double besy1(double z) {
  double y1;
  if (abs(z) < bessel_lim) {
    double d = z * 0.5;
    double a = 2.0 * besj1(z) * (eular * log(d));
    double b = d;
    double c = 1.0;
    double s = a - b * c - 1.0 / d;
    for (int k = 1; k < bessel_order; k++) {
      b = -b * (d * d) / (k * (k + 1));
      c += 1.0 / k + 1.0 / (k + 1.0);
      s -= b * c;
    }
    y1 = s / pi;
  } else {
    y1 = sqrt(2.0 / (pi * z)) * sin(z - 0.75 * pi);
  }
  return y1;
}
