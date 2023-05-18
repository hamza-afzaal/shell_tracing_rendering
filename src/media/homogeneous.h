
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

#ifndef PBRT_MEDIA_HOMOGENEOUS_H
#define PBRT_MEDIA_HOMOGENEOUS_H

// media/homogeneous.h*
#include <random>
#include "sampling.h"
#include "ShellPdf.h"
#include "medium.h"
#include "parallel.h"

namespace pbrt {

// HomogeneousMedium Declarations
class HomogeneousMedium : public Medium {
  public:
    // HomogeneousMedium Public Methods
    HomogeneousMedium(const Spectrum &sigma_a, const Spectrum &sigma_s, Float g,
                      Sampler &sampler)
        : sigma_a(sigma_a), sigma_s(sigma_s), sigma_t(sigma_s + sigma_a), g(g) {
	        //Preprocess(64, 4.f, 5.f, 1.f, 10000000, sampler);
    }
    Spectrum Tr(const Ray &ray, Sampler &sampler) const;
    Spectrum Sample(const Ray &ray, Sampler &sampler, MemoryArena &arena,
                    MediumInteraction *mi, Vector3f *rayOffset = NULL) const;
	const SpherePdf * getShellPdf(float radius) const;
    const vector<SpherePdf> *getShellPdf() const { return &sphere_pdf; }

    void Preprocess(int binSize, float minRadius, float maxRadius,
                    float stepSize, int sampleCount);

    int PopulatePdf(float minRadius, float maxRadius, float stepSize,
                     int binSize);

  private:
    enum RAY_STATE { RAY_ABSORBED, RAY_SCATTERED, RAY_EXITTED };

    // HomogeneousMedium Private Data
    const Spectrum sigma_a, sigma_s, sigma_t;
    const Float g;
    vector<SpherePdf> sphere_pdf;
};

}  // namespace pbrt

#endif  // PBRT_MEDIA_HOMOGENEOUS_H
