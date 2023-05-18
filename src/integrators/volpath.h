
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_VOLPATH_H
#define PBRT_INTEGRATORS_VOLPATH_H

// integrators/volpath.h*
#include <numeric>
#include "integrator.h"
#include "lightdistrib.h"
#include "media/ShellPdf.h"
#include "pbrt.h"
#include "transform.h"


namespace pbrt {

// VolPathIntegrator Declarations
class VolPathIntegrator : public SamplerIntegrator {
  public:
    // VolPathIntegrator Public Methods
    VolPathIntegrator(int maxDepth, std::shared_ptr<const Camera> camera,
                      std::shared_ptr<Sampler> sampler,
                      const Bounds2i &pixelBounds, Float rrThreshold = 1,
                      const std::string &lightSampleStrategy = "spatial")
        : SamplerIntegrator(camera, sampler, pixelBounds),
          maxDepth(maxDepth),
          rrThreshold(rrThreshold),
          lightSampleStrategy(lightSampleStrategy) {}
    void Preprocess(const Scene &scene, Sampler &sampler);
    Spectrum Li(const RayDifferential &ray, const Scene &scene,
                Sampler &sampler, MemoryArena &arena, int depth) const;

    float ComputeNearestDistance(const RayDifferential &ray,
                                 const Scene &scene) const;
    Transform ComputeShellRayTransform(Vector3f wo) const;
    int SampleFromDistribution(const vector<float> *alphaPdf, Float *prob,
                               Sampler &sampler) const;

    const vector<float> generatedPdf1D(const vector<vector<int>> *pdf,
                                       int alphaIndex,
                                       int binSize = 64,
                                       int countRegistered = 40000) const {
        vector<float> pdf1d(binSize, 0.f);
        vector<vector<int>> pdf2d = pdf[alphaIndex];

        for (int i = 0; i < pdf2d.size(); i++) {
            pdf1d[i] = std::accumulate(pdf2d[i].begin(), pdf2d[i].end(), 0.f);
        }

        float sumPdf = std::accumulate(pdf1d.begin(), pdf1d.end(), 0.f);
        for (int i = 0; i < pdf1d.size(); i++) {
            pdf1d[i] = pdf1d[i] / sumPdf;
        }
        return pdf1d;
    }

  private:
    // VolPathIntegrator Private Data
    const int maxDepth;
    const Float rrThreshold;
    const std::string lightSampleStrategy;
    std::unique_ptr<LightDistribution> lightDistribution;

    std::vector<Point3f> platonicSolidVertx = {
        Point3f(0.f, 1.f, 1.6180f),  Point3f(0.f, -1.f, 1.6180f),
        Point3f(0.f, 1.f, -1.6180f), Point3f(0.f, -1.f, -1.6180f),

        Point3f(1.f, 1.6180f, 0.f),  Point3f(-1.f, 1.6180f, 0.f),
        Point3f(1.f, -1.6180f, 0.f), Point3f(-1.f, -1.6180f, 0.f),

        Point3f(1.6180f, 0.f, 1.f),  Point3f(-1.6180f, 0.f, 1.f),
        Point3f(1.6180f, 0.f, -1.f), Point3f(-1.6180f, 0.f, -1.f)};
};

VolPathIntegrator *CreateVolPathIntegrator(
    const ParamSet &params, std::shared_ptr<Sampler> sampler,
    std::shared_ptr<const Camera> camera);

}  // namespace pbrt

#endif  // PBRT_INTEGRATORS_VOLPATH_H
