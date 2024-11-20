#include "Models.h"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <string>

using namespace std;

float fastPow(float base, int exp) {
  float result = 1.0f;
  while (exp > 0) {
    if (exp & 1)
      result *= base;
    base *= base;
    exp >>= 1;
  }
  return result;
}

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
bool verbose = false;
uint32_t RAND_SEED_P = 17;

// Statistical variables
// velocity sum
float vSum[2] = {0, 0};
// kinetic energy (Ek)
float keSum = 0;
float keSum2 = 0;
// total energy (E)
float totalEnergy = 0;
float totalEnergy2 = 0;
// pressure(P)
float pressure = 0;
float pressure2 = 0;

float rCut = 0;
float region[2] = {0, 0};
float velMag = 0;

// random
float random_r() {
  RAND_SEED_P = (RAND_SEED_P * IMUL + IADD) & MASK;
  return SCALE * RAND_SEED_P;
}

void random_velocity(float &v1, float &v2) {
  const float s = 2.0 * M_PI * random_r();
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
  readToken(file, "stepAvg");
  file >> config.stepAvg;
  readToken(file, "stepLimit");
  file >> config.stepLimit;
  readToken(file, "temperature");
  file >> config.temperature;
}

// Output result
void outputResult(const string &folder, const int n, const Molecule *molecules,
                  const int step, const float dTime) {
  if (!filesystem::exists(folder))
    filesystem::create_directories(folder);
  ofstream file;
  file.open(folder + "/" + to_string(step) + ".out");
  if (!file.is_open()) {
    std::cerr << "Error: file not found" << std::endl;
    exit(1);
  }
  file << "step " << to_string(step) << endl;
  file << "ts " << dTime << endl;
  file << setprecision(5) << fixed;
  file << "E " << totalEnergy << endl;
  file << "KE " << keSum << endl;
  file << "P " << pressure << endl;
  file << "sE " << totalEnergy2 << endl;
  file << "sKE " << keSum2 << endl;
  file << "sP " << pressure2 << endl;
  file << "====================" << endl;
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
                       float region[2], const float velMag, int size,
                       bool omp) {
  // create the folder if not exist
  if (!filesystem::exists("output"))
    filesystem::create_directories(string("output/") + (omp ? "omp" : "cpu") +
                                   "/" + to_string(size));

  ofstream file;
  file.open(string("output/") + (omp ? "omp/" : "cpu/") + to_string(size) +
            "/init");
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

void leapfrog_omp(const int n, Molecule *mols, const bool pre,
                  const float deltaT) {
#pragma omp parallel for schedule(dynamic)
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

void evaluateForce(const int n, Molecule *mols, float &uSum, float &virSum) {
  for (size_t i = 0; i < n - 1; i++) {
    for (size_t j = i + 1; j < n; j++) {
      // Make DeltaRij: (sum of squared RJ1-RJ2)
      float dr[2] = {mols[i].pos[0] - mols[j].pos[0],
                     mols[i].pos[1] - mols[j].pos[1]};
      toroidal(dr[0], dr[1], region);
      const float rr = dr[0] * dr[0] + dr[1] * dr[1];

      // case dr2 < Rc^2
      if (rr < rCut * rCut) {
        const float r = sqrt(rr);
        const float fcVal =
            48.0 * EPSILON * fastPow(SIGMA, 12) / fastPow(r, 13) -
            24.0 * EPSILON * fastPow(SIGMA, 6) / fastPow(r, 7);

        // update the acc
        mols[i].acc[0] += fcVal * dr[0];
        mols[i].acc[1] += fcVal * dr[1];
        mols[j].acc[0] -= fcVal * dr[0];
        mols[j].acc[1] -= fcVal * dr[1];

        // The completed Lennard-Jones.
        uSum +=
            4.0 * EPSILON * fastPow(SIGMA / r, 12) / r - fastPow(SIGMA / r, 6);

        virSum += fcVal * rr;
      }
    }
  }
}

void evaluateForce_omp(const int n, Molecule *mols, float &uSum,
                       float &virSum) {
#pragma omp parallel for reduction(+ : uSum, virSum) schedule(dynamic)
  for (size_t i = 0; i < n - 1; i++) {
    for (size_t j = i + 1; j < n; j++) {
      // Make DeltaRij: (sum of squared RJ1-RJ2)
      float dr[2] = {mols[i].pos[0] - mols[j].pos[0],
                     mols[i].pos[1] - mols[j].pos[1]};
      toroidal(dr[0], dr[1], region);
      const float rr = dr[0] * dr[0] + dr[1] * dr[1];

      // case dr2 < Rc^2
      if (rr < rCut * rCut) {
        const float r = sqrt(rr);
        const float fcVal =
            48.0 * EPSILON * fastPow(SIGMA, 12) / fastPow(r, 13) -
            24.0 * EPSILON * fastPow(SIGMA, 6) / fastPow(r, 7);
        // update the acc
#pragma omp atomic
        mols[i].acc[0] += fcVal * dr[0];
#pragma omp atomic
        mols[i].acc[1] += fcVal * dr[1];
#pragma omp atomic
        mols[j].acc[0] -= fcVal * dr[0];
#pragma omp atomic
        mols[j].acc[1] -= fcVal * dr[1];

        // The completed Lennard-Jones.
        uSum +=
            4.0 * EPSILON * fastPow(SIGMA / r, 12) / r - fastPow(SIGMA / r, 6);
        virSum += fcVal * rr;
      }
    }
  }
}

void stepSummary(const int n, const int step, const float dTime) {
  // cal avg and std of kinetic energy, total energy, and pressure
  float keAvg = keSum / config.stepAvg;
  float totalAvg = totalEnergy / config.stepAvg;
  float pressureAvg = pressure / config.stepAvg;

  float keStd = sqrt(max(0.0f, keSum2 / config.stepAvg - keAvg * keAvg));
  float totalStd =
      sqrt(max(0.0f, totalEnergy2 / config.stepAvg - totalAvg * totalAvg));
  float pressureStd =
      sqrt(max(0.0f, pressure2 / config.stepAvg - pressureAvg * pressureAvg));

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

void evaluateProperties(const int n, const Molecule *mols, const float &uSum,
                        const float &virSum, int step, float dTime,
                        bool summary) {

  vSum[0] = 0;
  vSum[1] = 0;

  float vvSum = 0;

  for (int i = 0; i < n; i++) {
    vSum[0] += mols[i].vel[0];
    vSum[1] += mols[i].vel[1];
    vvSum += mols[i].vel[0] * mols[i].vel[0] + mols[i].vel[1] * mols[i].vel[1];
  }

  const float ke = 0.5 * vvSum / n;
  const float energy = ke + uSum / n;
  const float p = config.density * (vvSum + virSum) / (n * 2);

  keSum += ke;
  totalEnergy += energy;
  pressure += p;

  keSum2 += ke * ke;
  totalEnergy2 += energy * energy;
  pressure2 += p * p;

  if (summary) {
    stepSummary(n, step, dTime);
  }
}

void launchSequentail(int N, Molecule *mols, const int size) {
  auto start_time = chrono::high_resolution_clock::now();

  int step = 0;
  while (step < config.stepLimit) {
    step++;
    const float deltaT = config.deltaT * step;
    float uSum = 0;
    float virSum = 0;
    leapfrog(N, mols, true, config.deltaT);
    evaluateForce(N, mols, uSum, virSum);
    leapfrog(N, mols, false, config.deltaT);
    evaluateProperties(N, mols, uSum, virSum, step, deltaT,
                       config.stepAvg > 0 && step % config.stepAvg == 0);
    if (verbose)
      outputResult("output/cpu/" + to_string(size), N, mols, step - 1,
                   config.deltaT);
  }

  auto end_time = chrono::high_resolution_clock::now();
  auto duration =
      chrono::duration_cast<chrono::microseconds>(end_time - start_time);
  cout << "[Seq Time] " << duration.count() << "ms - " << fixed
       << setprecision(4) << duration.count() / 1000000.0 << "s" << endl;
}

void launchOMP(int N, Molecule *mols, const int size) {
  auto start_time = chrono::high_resolution_clock::now();
  int step = 0;
  while (step < config.stepLimit) {
    step++;
    const float deltaT = config.deltaT * step;
    float uSum = 0;
    float virSum = 0;
    leapfrog_omp(N, mols, true, config.deltaT);
    evaluateForce_omp(N, mols, uSum, virSum);
    leapfrog_omp(N, mols, false, config.deltaT);
    evaluateProperties(N, mols, uSum, virSum, step, deltaT,
                       config.stepAvg > 0 && step % config.stepAvg == 0);

    if (verbose)
      outputResult("output/omp/" + to_string(size), N, mols, step - 1,
                   config.deltaT);
  }

  auto end_time = chrono::high_resolution_clock::now();
  auto duration =
      chrono::duration_cast<chrono::microseconds>(end_time - start_time);
  cout << "[OMP Time] " << duration.count() << "ms - " << fixed
       << setprecision(4) << duration.count() / 1000000.0 << "s" << endl;
}

// Main function
int main(const int argc, char *argv[]) {
  // Parse arguments
  if (argc < 4) {
    std::cerr << "Usage: " << argv[0] << " <config file> <size> <mode>"
              << std::endl;
    std::cerr << "  mode: 0 - Sequential, 1 - OMP" << std::endl;
    return 1;
  }

  const string filename = argv[1];
  const int size = atoi(argv[2]);
  const int mode = atoi(argv[3]);

  if (argc == 5) {
    verbose = atoi(argv[4]);
  }

  readConfig(filename);
  const int mSize = size * size;
  Molecule molecules[mSize];
  rCut = fastPow(2.0, 1.0 / 6.0 * SIGMA);

  // Region size
  region[0] = 1.0 / sqrt(config.density) * size;
  region[1] = 1.0 / sqrt(config.density) * size;

  // Velocity magnitude
  velMag = sqrt(NDIM * (1.0 - 1.0 / mSize) * config.temperature);
  const float gap[2] = {region[0] / size, region[1] / size};

  for (int y = 0; y < size; y++) {
    for (int x = 0; x < size; x++) {
      Molecule m;
      // assign molecule id
      m.id = y * size + x;

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
      molecules[y * size + x] = m;
    }
  }

  for (int i = 0; i < mSize; i++) {
    molecules[i].vel[0] -= vSum[0] / mSize;
    molecules[i].vel[1] -= vSum[1] / mSize;
  }

  if (verbose)
    outputMolInitData(mSize, molecules, rCut, region, velMag, size, mode == 1);
  if (mode == 0) {
    cout << "=========== Sequential Version ===========" << endl;
    cout << "Step\tTime\t\tvSum\t\tE.Avg\t\tE.Std\t\tK.Avg\t\tK.Std\t\tP."
            "Avg\t\tP.\t\tStd"
         << endl;
    launchSequentail(mSize, molecules, size);
  } else {
    cout << "=========== OMP Version ===========" << endl;
    cout << "Step\tTime\t\tvSum\t\tE.Avg\t\tE.Std\t\tK.Avg\t\tK.Std\t\tP."
            "Avg\t\tP.\t\tStd"
         << endl;
    launchOMP(mSize, molecules, size);
  }
  return 0;
}