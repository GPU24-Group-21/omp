#include "Models.h"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <omp.h>

using namespace std;

struct BlockResult {
  float vSum[2];
  float vvSum;
};

struct PropertiesData {
  float vSum[2];
  float vvSum;
  float ke;
  float keSum;
  float keSum2;
  float totalEnergy;
  float totalEnergy2;
  float pressure;
  float pressure2;
  int cycleCount; // 添加周期计数器
};

// Constants
Config config;

bool debug = false;
bool outfile = false;

uint32_t RAND_SEED_P = 17;

// Statistical variables
// velocity sum
double vSum[2] = {0, 0};
// kinetic energy (Ek)
double keSum = 0;
double keSum2 = 0;
// total energy (E)
double totalEnergy = 0;
double totalEnergy2 = 0;
// pressure(P)
double pressure = 0;
double pressure2 = 0;

float rCut = 0;
float region[2] = {0, 0};
float velMag = 0;

// random
double random_r() {
  RAND_SEED_P = (RAND_SEED_P * IMUL + IADD) & MASK;
  return SCALE * RAND_SEED_P;
}

void random_velocity(float &v1, float &v2) {
  const double s = 2.0 * M_PI * random_r();
  v1 = cos(s);
  v2 = sin(s);
}

// Read input file
void readToken(std::ifstream &file, const string &token) {
  string str;
  file >> str;
  if (str != token) {
    std::cerr << "Error: token not found: " << token << std::endl;
    exit(1);
  }
}

void readConfig(const string &filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error: file not found" << std::endl;
    exit(1);
  }

  // Read config
  readToken(file, "deltaT");
  file >> config.deltaT;
  readToken(file, "density");
  file >> config.density;
  readToken(file, "initUcell_x");
  file >> config.initUcell_x;
  readToken(file, "initUcell_y");
  file >> config.initUcell_y;
  readToken(file, "stepAvg");
  file >> config.stepAvg;
  readToken(file, "stepEquil");
  file >> config.stepEquil;
  readToken(file, "stepLimit");
  file >> config.stepLimit;
  readToken(file, "temperature");
  file >> config.temperature;

  if (debug) {
    cout << "=========== Config ===========" << endl;
    cout << "  deltaT: " << config.deltaT << endl;
    cout << "  density: " << config.density << endl;
    cout << "  initUcell_x: " << config.initUcell_x << endl;
    cout << "  initUcell_y: " << config.initUcell_y << endl;
    cout << "  stepAvg: " << config.stepAvg << endl;
    cout << "  stepEquil: " << config.stepEquil << endl;
    cout << "  stepLimit: " << config.stepLimit << endl;
    cout << "  temperature: " << config.temperature << endl;
  }
}

void readMoo(const string &filename, long N, Molecule *molecules) {
  std::ifstream file(filename);

  if (!file.is_open()) {
    std::cerr << "Error: file not found" << std::endl;
    exit(1);
  }

  readToken(file, "rCut");
  file >> rCut;
  readToken(file, "region");
  file >> region[0] >> region[1];
  readToken(file, "velMag");
  file >> velMag;

  if (debug) {
    cout << "=========== Pre Defined Props ===========" << endl;
    cout << "  rCut: " << rCut << endl;
    cout << "  region: " << region[0] << " " << region[1] << endl;
    cout << "  velMag: " << velMag << endl;
  }

  // rest til the end of the file is the molecule
  for (int i = 0; i < N; i++) {
    Molecule m;
    m.id = i;
    file >> m.pos[0] >> m.pos[1] >> m.vel[0] >> m.vel[1] >> m.acc[0] >>
        m.acc[1];
    molecules[i] = m;
  }

  if (debug) {
    cout << "=========== Molecules ===========" << endl;
    for (int i = 0; i < N; i++) {
      cout << "  id: " << molecules[i].id << " pos: " << molecules[i].pos[0]
           << " " << molecules[i].pos[1] << " vel: " << molecules[i].vel[0]
           << " " << molecules[i].vel[1] << endl;
    }
  }
}

// Output result
void outputResult(const string &filename, const int n,
                  const Molecule *molecules, const int step,
                  const double dTime) {
  ofstream file;
  file.open(filename);

  if (!file.is_open()) {
    std::cerr << "Error: file not found" << std::endl;
    exit(1);
  }

  file << "step: " << to_string(step) << endl;
  file << "ts: " << dTime << endl;
  file << "====================" << endl;
  file << setprecision(5) << fixed;
  const int mark1 = n / 2 + n / 8;
  const int mark2 = n / 2 + n / 8 + 1;
  for (int i = 0; i < n; i++) {
    if (i == mark1 || i == mark2) {
      file << "m-" << molecules[i].id << " ";
    } else {
      file << "o-" << molecules[i].id << " ";
    }
    file << molecules[i].pos[0] << " " << molecules[i].pos[1] << endl;
  }
}

void outputMolInitData(const int n, const Molecule *molecules, const float rCut,
                       float region[2], const float velMag) {
  ofstream file;
  file.open("mols.in");
  file << setprecision(5) << fixed;
  file << "rCut " << rCut << endl;
  file << "region " << region[0] << " " << region[1] << endl;
  file << "velMag " << velMag << endl;
  file << "vSum " << vSum[0] << " " << vSum[1] << endl;
  for (int i = 0; i < n; i++) {
    file << molecules[i].pos[0] << " " << molecules[i].pos[1] << " "
         << molecules[i].vel[0] << " " << molecules[i].vel[1] << " "
         << molecules[i].acc[0] << " " << molecules[i].acc[1] << endl;
  }
}

// ============= Core function =============
// Toroidal functions
void toroidal(float &x, float &y, const float region[2]) {
  if (x < -0.5 * region[0])
    x += region[0];
  if (x >= 0.5 * region[0])
    x -= region[0];
  if (y < -0.5 * region[1])
    y += region[1];
  if (y >= 0.5 * region[1])
    y -= region[1];
}

void leapfrog(const int n, Molecule *mols, const bool pre, const float deltaT) {
  for (int i = 0; i < n; i++) {
    // v(t + Δt/2) = v(t) + (Δt/2)a(t)
    mols[i].vel[0] += 0.5 * deltaT * mols[i].acc[0];
    mols[i].vel[1] += 0.5 * deltaT * mols[i].acc[1];

    if (pre) {
      // r(t + Δt) = r(t) + Δt v(t + Δt/2)
      mols[i].pos[0] += deltaT * mols[i].vel[0];
      mols[i].pos[1] += deltaT * mols[i].vel[1];

      toroidal(mols[i].pos[0], mols[i].pos[1], region);

      mols[i].acc[0] = 0;
      mols[i].acc[1] = 0;
    }
  }
}

void evaluateForce(const int n, Molecule *mols, double &uSum, double &virSum) {
  for (size_t i = 0; i < n - 1; i++) {
    for (size_t j = i + 1; j < n; j++) {
      // Make DeltaRij: (sum of squared RJ1-RJ2)
      float dr[2] = {mols[i].pos[0] - mols[j].pos[0],
                     mols[i].pos[1] - mols[j].pos[1]};
      toroidal(dr[0], dr[1], region);
      const double rr = dr[0] * dr[0] + dr[1] * dr[1];

      // case dr2 < Rc^2
      if (rr < rCut * rCut) {
        const double r = sqrt(rr);
        const double fcVal = 48.0 * EPSILON * pow(SIGMA, 12) / pow(r, 13) -
                             24.0 * EPSILON * pow(SIGMA, 6) / pow(r, 7);
        // update the acc
        mols[i].acc[0] += fcVal * dr[0];
        mols[i].acc[1] += fcVal * dr[1];
        mols[j].acc[0] -= fcVal * dr[0];
        mols[j].acc[1] -= fcVal * dr[1];

        // The completed Lennard-Jones.
        uSum += 4.0 * EPSILON * pow(SIGMA / r, 12) / r - pow(SIGMA / r, 6);
        virSum += fcVal * rr;
      }
    }
  }
}

void evaluateProperties(const int n, const Molecule *mols, const double &uSum,
                        const double &virSum) {

  vSum[0] = 0;
  vSum[1] = 0;

  double vvSum = 0;

  for (int i = 0; i < n; i++) {
    vSum[0] += mols[i].vel[0];
    vSum[1] += mols[i].vel[1];
    vvSum += mols[i].vel[0] * mols[i].vel[0] + mols[i].vel[1] * mols[i].vel[1];
  }

  const double ke = 0.5 * vvSum / n;

  const double energy = ke + uSum / n;
  const double p = config.density * (vvSum + virSum) / (n * 2);

  keSum += ke;
  totalEnergy += energy;
  pressure += p;

  keSum2 += ke * ke;
  totalEnergy2 += energy * energy;
  pressure2 += p * p;
}

void stepSummary(const int n, const int step, const double dTime) {
  // average and standard deviation of kinetic energy, total energy, and
  // pressure
  double keAvg = keSum / config.stepAvg;
  double totalAvg = totalEnergy / config.stepAvg;
  double pressureAvg = pressure / config.stepAvg;

  double keStd = sqrt(max(0.0, keSum2 / config.stepAvg - keAvg * keAvg));
  double totalStd =
      sqrt(max(0.0, totalEnergy2 / config.stepAvg - totalAvg * totalAvg));
  double pressureStd =
      sqrt(max(0.0, pressure2 / config.stepAvg - pressureAvg * pressureAvg));

  cout << fixed << setprecision(8) << step << "\t" << dTime << "\t"
       << vSum[0] / n << "\t" << totalAvg << "\t" << totalStd << "\t" << keAvg
       << "\t" << keStd << "\t" << pressureAvg << "\t" << pressureStd << endl;

  // reset the sum
  keSum = 0;
  keSum2 = 0;
  totalEnergy = 0;
  totalEnergy2 = 0;
  pressure = 0;
  pressure2 = 0;
}

void launchSequentail(int N, Molecule *mols) {
  auto start_time = chrono::high_resolution_clock::now();

  int step = 0;
  while (step < config.stepLimit) {
    step++;
    const double deltaT = static_cast<double>(step) * config.deltaT;
    double uSum = 0;
    double virSum = 0;

    leapfrog(N, mols, true, config.deltaT);
    boundaryCondition(N, mols);
    evaluateForce(N, mols, uSum, virSum);
    leapfrog(N, mols, false, config.deltaT);
    evaluateProperties(N, mols, uSum, virSum);
    if (config.stepAvg > 0 && step % config.stepAvg == 0) {
      stepSummary(N, step, deltaT);
    }
    // output the result
    if (outfile) {
      outputResult("output/" + to_string(step - 1) + ".out", N, mols, step - 1,
                   config.deltaT);
    }
  }

  auto end_time = chrono::high_resolution_clock::now();
  auto duration =
      chrono::duration_cast<chrono::microseconds>(end_time - start_time);
  cout << "[CPU Time] " << duration.count() << "ms - " << fixed
       << setprecision(4) << duration.count() / 1000000.0 << "s" << endl;
}

void launchOMP(int N, Molecule *mols) {
  auto start_time = chrono::high_resolution_clock::now();
  int step = 0;
  while (step < config.stepLimit) {
    step++;
    const double deltaT = static_cast<double>(step) * config.deltaT;
    double uSum = 0;
    double virSum = 0;

    leapfrog(N, mols, true, config.deltaT);
    evaluateForce(N, mols, uSum, virSum);
    leapfrog(N, mols, false, config.deltaT);
    evaluateProperties(N, mols, uSum, virSum);
    if (config.stepAvg > 0 && step % config.stepAvg == 0) {
      stepSummary(N, step, deltaT);
    }
    if (outfile) {
      outputResult("output/" + to_string(step - 1) + ".out", N, mols, step - 1,
                   config.deltaT);
    }
  }

  auto end_time = chrono::high_resolution_clock::now();
  auto duration =
      chrono::duration_cast<chrono::microseconds>(end_time - start_time);
  cout << "[CPU Time] " << duration.count() << "ms - " << fixed
       << setprecision(4) << duration.count() / 1000000.0 << "s" << endl;
}

/*
  Validate the result
*/
void validate(const int N, const Molecule *mols, const Molecule *mols2) {
  bool error = false;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < 2; j++) {
      if (fabs(mols[i].pos[j] - mols2[i].pos[j]) >= 0.009) {
        std::cerr << "Error: pos mismatch: mols[" << i << "].pos[" << j
                  << "] = " << mols[i].pos[j] << " != " << mols2[i].pos[j]
                  << endl;
        error = error || true;
      }
      if (fabs(mols[i].vel[j] - mols2[i].vel[j]) >= 0.009) {
        std::cerr << "Error: vel mismatch: mols[" << i << "].vel[" << j
                  << "] = " << mols[i].vel[j] << " != " << mols2[i].vel[j]
                  << endl;
        error = error || true;
      }
      if (fabs(mols[i].acc[j] - mols2[i].acc[j]) >= 0.009) {
        std::cerr << "Error: acc mismatch: mols[" << i << "].acc[" << j
                  << "] = " << mols[i].acc[j] << " != " << mols2[i].acc[j]
                  << endl;
        error = error || true;
      }
    }
  }

  if (!error) {
    cout << "Validation passed" << endl;
  }
}

// Main function
int main(const int argc, char *argv[]) {
  // Parse arguments
  if (argc < 4) {
    std::cerr << "Usage: " << argv[0]
              << " <config file> <0: no file output, 1:output step file>"
              << std::endl;
    return 1;
  }

  const string filename = argv[1];
  readConfig(filename);
  system("mkdir -p output");

  const int mSize = config.initUcell_x * config.initUcell_y;
  Molecule molecules_c[mSize], molecules_o[mSize];
  cout << "Size: " << config.initUcell_x << "x" << config.initUcell_y << "("
       << mSize << " mols)" << endl;
  outfile = atoi(argv[2]);
  const int mode = atoi(argv[3]);
  rCut = pow(2.0, 1.0 / 6.0 * SIGMA);
  
  // Region size
  region[0] = 1.0 / sqrt(config.density) * config.initUcell_x;
  region[1] = 1.0 / sqrt(config.density) * config.initUcell_y;

  // Velocity magnitude
  velMag = sqrt(NDIM * (1.0 - 1.0 / mSize) * config.temperature);

  if (debug) {
    cout << "=========== Random Init ===========" << endl;
    cout << "  rCut: " << rCut << endl;
    cout << "  region: " << region[0] << " " << region[1] << endl;
    cout << "  velMag: " << velMag << endl;
  }

  const double gap[2] = {region[0] / config.initUcell_x,
                         region[1] / config.initUcell_y};

  for (int y = 0; y < config.initUcell_x; y++) {
    for (int x = 0; x < config.initUcell_y; x++) {
      Molecule m;
      // assign molecule id
      m.id = y * config.initUcell_y + x;

      // assign position
      m.pos[0] = (x + 0.5) * gap[0] + region[0] * -0.5;
      m.pos[1] = (y + 0.5) * gap[1] + region[1] * -0.5;

      // assign velocity
      random_velocity(m.vel[0], m.vel[1]);
      m.multiple_vel(velMag);

      // update the vsum
      vSum[0] += m.vel[0];
      vSum[1] += m.vel[1];

      // assign acceleration
      m.multiple_acc(0);

      // add to list
      molecules_c[y * config.initUcell_y + x] = m;
      molecules_o[y * config.initUcell_y + x] = m;
    }
  }

  for (int i = 0; i < mSize; i++) {
    molecules_c[i].vel[0] -= vSum[0] / mSize;
    molecules_c[i].vel[1] -= vSum[1] / mSize;

    molecules_o[i].vel[0] -= vSum[0] / mSize;
    molecules_o[i].vel[1] -= vSum[1] / mSize;
  }

  if (outfile) {
    outputMolInitData(mSize, molecules_c, rCut, region, velMag);
  }

  cout << "=========== Sequential Version ===========" << endl;
  // print the header
  cout
      << "Step\tTime\t\tvSum\t\tE.Avg\t\tE.Std\t\tK.Avg\t\tK.Std\t\tP.Avg\t\tP."
         "Std"
      << endl;
  launchSequentail(mSize, molecules_c);
  cout << "=========== OMP Version ===========" << endl;
  cout
      << "Step\tTime\t\tvSum\t\tE.Avg\t\tE.Std\t\tK.Avg\t\tK.Std\t\tP.Avg\t\tP."
         "Std"
      << endl;
  launchKernel(mSize, molecules_o);

  // validate the result
  cout << "=========== Validation ===========" << endl;
  validate(mSize, molecules_c, molecules_o);
  cout << "=========== Done ===========" << endl;
  return 0;
}