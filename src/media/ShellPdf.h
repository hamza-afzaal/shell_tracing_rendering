#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_SHELL_PDF_H
#define PBRT_SHELL_PDF_H

#include <algorithm>
#include <numeric>
#include <vector>
#include "geometry.h"

using namespace std;
using std::acos;
// float toRadians = 3.1415 / 180.0;

namespace pbrt {

struct PdfBin {
    int alpha = 0;
    int theta = 0;
    int phi = 0;
};

class ShellPdf {
  public:
    ShellPdf(){};
    ShellPdf(int binSize);

    void populatePdf(float alpha, float theta, float phi);
    int searchBinIndex(const vector<float> &bin, float val);
    const vector<float> generatedPdf1D() const;
    const vector<float> generatedPdf1D(int alphaIndex) const;
    const vector<float> generatedPdf1D(int alphaIndex, int thetaIndex) const;

	const void PrecomputeAlphaPdf();
	const void PrecomputeThetaPdf();
	const void PrecomputePhiPdf();
    
    float AlphaBinValue(int index) const;
    float ThetaBinValue(int index) const;
    float PhiBinValue(int index) const;


	int SampleCount() const { return countRegistered; }
    int BinSize() const { return binSize; }
    const vector<vector<vector<int>>> *getPdf() const { return &pdf; }

    const vector<float> *GetAlphaPdf() const { return &alphaPdf; }
    const vector<float> *GetThetaPdf(int alphaIndex) const { return &(thetaPdf[alphaIndex]); }
    const vector<float> *GetPhiPdf(int alphaIndex, int thetaIndex) const { return &(phiPdf[alphaIndex][thetaIndex]); }

  private:
    vector<vector<vector<int>>> pdf;

	vector<float> alphaPdf;
    vector<vector<float>> thetaPdf;
    vector<vector<vector<float>>> phiPdf;

    vector<float> alphaBin;
    vector<float> thetaBin;
    vector<float> phiBin;
    int binSize;

    int countRegistered = 0;
};

struct SpherePdf {
    float radius;
    ShellPdf pdf;
};

// free functions
float computeRaySphereIntersect(Point3f center, Float radius, Ray &r);

}  // namespace pbrt

#endif
