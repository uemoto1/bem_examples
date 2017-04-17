//      bem2d.cpp
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
#include "bem2d.hpp"

dcmplx h0(double x) { return dcmplx(besj0(x), besy0(x)); }

dcmplx green(double k, const vec2d &r) {
  return dcmplx(0.0, 0.25) * h0(k * r.norm());
}

CIncident::CIncident() {
  a = 1.0;
  k = 1.0;
  w = 1.0;
  ivec = vec2d(1.0, 0.0);
}

CIncident::CIncident(const vec2d &kvec, double omega, const dcmplx &amp) {
  a = amp;
  k = kvec.norm();
  w = omega;
  ivec = kvec / k;
}

dcmplx CIncident::field(const vec2d &r) {
  return a * exp(dcmplx(0.0, k * ivec.dot(r)));
}

CObject::CObject() {}

double CObject::length(int i) {
  int j = (i + 1) % node.size();
  return (node[i] - node[j]).norm();
}

bool CObject::inside(const vec2d &r) {
  int n = node.size();
  double x = r(0), y = r(1);

  int count = 0;
  for (int i = 0; i < n; i++) {
    int j = (i + 1) % n;
    vec2d ri = node[i];
    vec2d rj = node[j];
    double xi = ri(0), yi = ri(1);
    double xj = rj(0), yj = rj(1);

    if (((yi <= y) && (y < yj)) || ((yj <= y) && (y < yi))) {
      double dx = xj - xi;
      double dy = yj - yi;
      double fy = y - yi;
      double xc = xi + dx / dy * fy;
      if (x < xc) {
        count++;
      }
    }
  }
  return ((0 < count) && (count % 2 != 0));
}

void bem2d(CObject &obj, CIncident &inc) {
  int N = obj.node.size();
  double k = inc.k;
  double w = inc.w;

  Eigen::MatrixXcd G(N, N);
  Eigen::MatrixXcd A(N, N);
  Eigen::VectorXcd b(N);

  dcmplx fa(1.0, 0.073804);
  dcmplx fb(0.0, 0.636619);

  // Generating Green's function
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      if (i == j) {
        G(i, j) = 0.0;
      } else {
        vec2d ri = obj.node[i];
        vec2d rj = obj.node[j];
        G(i, j) = green(k, ri - rj);
      }
    }
  }

  // Compose BEM matrix
  for (int i = 0; i < N; i++) {
    vec2d ri = obj.node[i];
    b[i] = inc.field(ri);

    for (int j = 0; j < N; j++) {
      int jp = (j + 1) % N;
      int jm = (j + N - 1) % N;
      double ljm = obj.length(jm);
      double lj = obj.length(j);

      if (i == j) {
        dcmplx sjm = 0.5 * fa - 0.75 * fb + 0.5 * fb * log(k * ljm);
        dcmplx sj = 0.5 * fa - 0.75 * fb + 0.5 * fb * log(k * lj);
        A(i, i) = dcmplx(0.0, 0.25) * (ljm * sjm + lj * sj);
      } else {
        dcmplx sjm = (ljm / 6) * G(i, jm);
        dcmplx sj = (ljm / 3 + lj / 3) * G(i, j);
        dcmplx sjp = (lj / 6) * G(i, jp);
        A(i, j) = sjm + sj + sjp;
      }
    }
  }
  cout << "# Solving linear equation..." << endl;
  // Solve DDA Matrix
  Eigen::BiCGSTAB<Eigen::MatrixXcd> solver(A);

  // int step = 0;
  obj.data = solver.solveWithGuess(b, obj.data); // Initial state
  // solver.setMaxIterations(iter_step);
  // do {
  //  grid.data = solver.solveWithGuess(b, grid.data);
  //  step += iter_step;
  //  cerr << "# Step:" << step << " Error:" << solver.error() <<
  // endl;
  //} while (solver.info() != Eigen::Success && step < max_step);

  // obj.data = A.householderQr().solve(b);
}

dcmplx field(const vec2d &r, CObject &obj, CIncident &inc) {
  int n = obj.node.size();
  double k = inc.k;
  dcmplx ef = inc.field(r);

  for (int i = 0; i < n; i++) {
    int j = (i + 1) % n;
    double l = obj.length(i);
    dcmplx fi = obj.data(i);
    dcmplx fj = obj.data(j);
    dcmplx fij = (fi + fj) * 0.5;
    vec2d ri = obj.node[i];
    vec2d rj = obj.node[j];
    vec2d rij = (ri + rj) * 0.5;
    dcmplx g = green(k, r - rij);
    ef -= l * g * fij;
  }
  return ef;
}

int main() {
  CIncident inc;
  CObject obj;

  string line;
  while (getline(cin, line)) {
    stringstream ss(line);
    string temp;
    vector<double> buff;
    while (ss >> temp)
      buff.push_back(stod(temp));

    if (2 <= buff.size()) {
      vec2d r;
      r(0) = buff[0];
      r(1) = buff[1];
      obj.node.push_back(r);
    }
  }

  obj.data = vecx::Zero(obj.node.size());

  bem2d(obj, inc);

  for (int i = -150; i <= 150; i++) {
    //int j = 0;
    for (int j = -150; j <= 150; j++) {
      dcmplx ef(0.0, 0.0);
      vec2d r(i * 0.1, j * 0.1);
      if (obj.inside(r)) {
        ef = 0.0;
      } else {
        ef = field(r, obj, inc);
      }
      printf("%f\t%f\t%e\t%e\n", r(0), r(1), ef.real(), ef.imag());
    }
  }

  return 0;
}
