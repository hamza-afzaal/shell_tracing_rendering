#include "ShellPdf.h"

pbrt::ShellPdf::ShellPdf(int binSize) : binSize(binSize) {
    float acosDiscSplit = 1.f / float(binSize);
    float discSplit = 2.f * 3.1415 / float(binSize);  // check this for -1 too

    // PdfBin bin;
    vector<vector<vector<int>>> pdf(
        binSize, vector<vector<int>>(binSize, vector<int>(binSize, 0)));

    this->pdf = pdf;

    for (int i = binSize - 1; i >= 0; --i) {
        alphaBin.push_back(acos(2.f * (i * acosDiscSplit) -
                                1.f));  // see if you need to convert in degrees
        thetaBin.push_back(acos(2.f * (i * acosDiscSplit) - 1.f));
        phiBin.push_back((binSize - i) * discSplit);
    }
}

void pbrt::ShellPdf::populatePdf(float alpha, float theta, float phi) {
    int alphaIndex = searchBinIndex(alphaBin, alpha);
    int thetaIndex = searchBinIndex(thetaBin, theta);
    int phiIndex = searchBinIndex(phiBin, phi);

    if (alphaIndex != -1 || thetaIndex != -1 || phiIndex != -1) {
        pdf[alphaIndex][thetaIndex][phiIndex] += 1;
        countRegistered++;
    }
}

// test this
int pbrt::ShellPdf::searchBinIndex(vector<float> const &bin, float val) {
    auto const it = std::upper_bound(bin.begin(), bin.end(), val);

    if (it >= bin.end()) return bin.size() - 1;

    return it - bin.begin();
}

const vector<float> pbrt::ShellPdf::generatedPdf1D() const {
    vector<float> pdf1d(binSize, 0.f);
    vector<vector<float>> pdf2d(binSize, vector<float>(binSize, 0));
    for (int i = 0; i < pdf.size(); i++) {  // need to invert this coz push_back
        for (int j = 0; j < pdf[i].size(); j++) {
            pdf2d[i][j] =
                std::accumulate(pdf[i][j].begin(), pdf[i][j].end(), 0.f);
        }
    }
    for (int i = 0; i < pdf2d.size(); i++) {
        pdf1d[i] = std::accumulate(pdf2d[i].begin(), pdf2d[i].end(), 0.f) /
                   countRegistered;
    }

    float accum = std::accumulate(pdf1d.begin(), pdf1d.end(), 0.f);
    return pdf1d;
}

const vector<float> pbrt::ShellPdf::generatedPdf1D(int alphaIndex) const {
    vector<float> pdf1d(binSize, 0.f);
    vector<vector<int>> pdf2d = pdf[alphaIndex];

    for (int i = 0; i < pdf2d.size(); i++) {
        pdf1d[i] = std::accumulate(pdf2d[i].begin(), pdf2d[i].end(), 0.f);
    }

    float sumPdf = std::accumulate(pdf1d.begin(), pdf1d.end(), 0.f);
    for (int i = 0; i < pdf1d.size(); i++) {
        pdf1d[i] = pdf1d[i] / sumPdf;
    }

    float accum = std::accumulate(pdf1d.begin(), pdf1d.end(), 0.f);

    return pdf1d;
}

const vector<float> pbrt::ShellPdf::generatedPdf1D(int alphaIndex,
                                                   int thetaIndex) const {
    vector<int> trimPdf = pdf[alphaIndex][thetaIndex];
    vector<float> pdf1d(binSize, 0.f);

    float sumPdf = std::accumulate(trimPdf.begin(), trimPdf.end(), 0.f);
    for (int i = 0; i < trimPdf.size(); i++) {
        pdf1d[i] = trimPdf[i] / sumPdf;
    }

    float accum = std::accumulate(pdf1d.begin(), pdf1d.end(), 0.f);

    return pdf1d;
}

const void pbrt::ShellPdf::PrecomputeAlphaPdf() {
    alphaPdf = vector<float>(binSize, 0.f);
    vector<vector<float>> pdf2d(binSize, vector<float>(binSize, 0));
    for (int i = 0; i < pdf.size(); i++) {  // need to invert this coz push_back
        for (int j = 0; j < pdf[i].size(); j++) {
            pdf2d[i][j] =
                std::accumulate(pdf[i][j].begin(), pdf[i][j].end(), 0.f);
        }
    }
    for (int i = 0; i < pdf2d.size(); i++) {
        alphaPdf[i] = std::accumulate(pdf2d[i].begin(), pdf2d[i].end(), 0.f) /
                      countRegistered;
    }
}

const void pbrt::ShellPdf::PrecomputeThetaPdf() {
    thetaPdf = vector<vector<float>>(binSize, vector<float>(binSize, 0.f));

    for (int i = 0; i < pdf.size(); i++) {  // need to invert this coz push_back
        for (int j = 0; j < pdf[i].size(); j++) {
            thetaPdf[i][j] =
                std::accumulate(pdf[i][j].begin(), pdf[i][j].end(), 0.f);
        }
    }
    for (int i = 0; i < thetaPdf.size(); i++) {
        float sum =
            std::accumulate(thetaPdf[i].begin(), thetaPdf[i].end(), 0.f);
        for (int j = 0; j < thetaPdf[i].size(); j++) {
            thetaPdf[i][j] = thetaPdf[i][j] / sum;
        }
    }
}

const void pbrt::ShellPdf::PrecomputePhiPdf() {
    phiPdf = vector<vector<vector<float>>>(
        binSize, vector<vector<float>>(binSize, vector<float>(binSize, 0)));

    for (int i = 0; i < pdf.size(); i++) {  // need to invert this coz push_back
        for (int j = 0; j < pdf[i].size(); j++) {
            float sum =
                std::accumulate(pdf[i][j].begin(), pdf[i][j].end(), 0.f);
			for (int k = 0; k < pdf[j].size(); k++)
			{
                phiPdf[i][j][k] = pdf[i][j][k] / sum;
			}
        }
    }
}

float pbrt::ShellPdf::AlphaBinValue(int index) const {
    if (index < 0 && index >= alphaBin.size()) return 0.f;
    return alphaBin[index];
}

float pbrt::ShellPdf::ThetaBinValue(int index) const {
    if (index < 0 && index >= thetaBin.size()) return 0.f;
    return thetaBin[index];
}

float pbrt::ShellPdf::PhiBinValue(int index) const {
    if (index < 0 && index >= phiBin.size()) return 0.f;
    return phiBin[index];
}

float pbrt::computeRaySphereIntersect(Point3f center, Float radius, Ray &r) {
    Vector3f D = r.o - center;  // direction of ray is r.d
    float delta = (pow(Dot(r.d, D), 2) - D.LengthSquared() + pow(radius, 2)) /
                  r.d.LengthSquared();
    if (delta > 0.f) {
        return -Dot(r.d, D) + sqrt(delta);
    } else if (delta == 0)
        return 0.f;
    else
        return -1.f;
}
