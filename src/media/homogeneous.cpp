
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

// media/homogeneous.cpp*
#include "media/homogeneous.h"
#include "interaction.h"
#include "paramset.h"
#include "sampler.h"
#include "stats.h"

namespace pbrt {

std::uniform_real_distribution<float> scopeRandomFloats(0.f, 1.f);
std::default_random_engine randGen;

// HomogeneousMedium Method Definitions
Spectrum HomogeneousMedium::Tr(const Ray &ray, Sampler &sampler) const {
    ProfilePhase _(Prof::MediumTr);
    return Exp(-sigma_t * std::min(ray.tMax * ray.d.Length(), MaxFloat));
}

// sample the homogeneous medium -> produces the direction for the ray going
// out??
Spectrum HomogeneousMedium::Sample(const Ray &ray, Sampler &sampler,
                                   MemoryArena &arena, MediumInteraction *mi,
                                   Vector3f *rayOffset) const {
    ProfilePhase _(Prof::MediumSample);
    // Sample a channel and distance along the ray
    int channel = std::min((int)(sampler.Get1D() * Spectrum::nSamples),
                           Spectrum::nSamples - 1);
    Float dist = -std::log(1 - sampler.Get1D()) / sigma_t[channel];
    Float t =
        std::min(rayOffset ? dist / ray.d.Length() : ray.d.Length(), ray.tMax);
    bool sampledMedium = rayOffset
                             ? Distance(ray.o, ray.o + *rayOffset) < ray.tMax
                             : t < ray.tMax;
    //if (sampledMedium)
    //    printf("The tMax is %f and bool %f\n", ray.tMax,
    //           rayOffset ? Distance(ray.o, ray.o + *rayOffset)
    //                     : 0.f);
    if (sampledMedium) //rayOffset ? ray.o + *rayOffset : 
        *mi = MediumInteraction(ray(t), -ray.d,
                                ray.time, this,
                                ARENA_ALLOC(arena, HenyeyGreenstein)(g));

    // Compute the transmittance and sampling density
    Spectrum Tr = Exp(-sigma_t * std::min(t, MaxFloat) * ray.d.Length());

    // Return weighting factor for scattering from homogeneous medium
    Spectrum density = sampledMedium ? (sigma_t * Tr) : Tr;
    Float pdf = 0;
    for (int i = 0; i < Spectrum::nSamples; ++i) pdf += density[i];
    pdf *= 1 / (Float)Spectrum::nSamples;
    if (pdf == 0) {
        CHECK(Tr.IsBlack());
        pdf = 1;
    }
    return sampledMedium ? (Tr * sigma_s / pdf) : (Tr / pdf);
}

const SpherePdf * HomogeneousMedium::getShellPdf(float radius) const {
    if (sphere_pdf.size() == 0) return NULL;

    int samp = 0;
    for (int i = 0; i < sphere_pdf.size(); i++) {
        if (radius > sphere_pdf[i].radius)
            samp = i;
        else
            break;
    }

    return &sphere_pdf[samp];

    // auto const it = std::upper_bound(
    //    sphere_pdf.begin(), sphere_pdf.end(), radius,
    //    [](double value, const SpherePdf &pdf) { return value < pdf.radius;
    //    });

    // if (it >= sphere_pdf.end()) {
    //    return sphere_pdf[sphere_pdf.size() - 1];
    //}

    // return sphere_pdf[it - sphere_pdf.begin()];
}

void HomogeneousMedium::Preprocess(int binSize, float minRadius,
                                   float maxRadius, float stepSize,
                                   int sampleCount) {
    // Do not use the sampler here, just use a random number generator to ease
    // your life a bit....

    float currRadius = minRadius;
    PhaseFunction *phase = new HenyeyGreenstein(g);
    int count = PopulatePdf(minRadius, maxRadius, stepSize, binSize);

    //(maxRadius - minRadius) / stepSize;

    // populate first then use the datastructure to update values
    // for (float i = 0; i < count; i++) {
    //    Point2i radius(i, 0);

    ParallelFor2D(
        [&](Point2i radius) {
            for (int i = 0; i < sampleCount; i++) {
                Ray r(Point3f(0.f, 0.f, 0.f),
                      Vector3f(0.f, 0.f,
                               1.f));  // initialized a ray pointing upwards
                RAY_STATE state = RAY_SCATTERED;
                bool firstScatter = true;
                while (state == RAY_SCATTERED) {
                    float t_exit = computeRaySphereIntersect(
                        Point3f(0.f, 0.f, 0.f), currRadius, r);

                    int channel = std::min(
                        (int)(scopeRandomFloats(randGen) *
                              Spectrum::nSamples),  // sampler is not working
                                                    // -> check this
                        Spectrum::nSamples - 1);
                    float t_scatter = -std::log(scopeRandomFloats(randGen)) /
                                      sigma_s[channel];  // 1 - sampler.Get1D()
                    float t_absorb = -std::log(scopeRandomFloats(randGen)) /
                                     sigma_a[channel];  // 1 - sampler.Get1D()

                    if (min(t_absorb, t_scatter) >=
                        t_exit) {  //  || Distance(Point3f(0.f, 0.f, 0.f), r.o)
                                   //  < currRadius
                        r = Ray(r(t_exit), r.d);
                        state = RAY_EXITTED;
                    } else {
                        firstScatter = false;
                        if (t_absorb < t_scatter)
                            state = RAY_ABSORBED;
                        else {
                            Vector3f wo = -r.d, wi;
                            phase->Sample_p(
                                wo, &wi,
                                Point2f(scopeRandomFloats(randGen),
                                        scopeRandomFloats(randGen)));
                            Point3f p = r.o + (r.d * t_scatter);
                            r = Ray(p, wi);
                        }
                    }
                }

                if (state == RAY_EXITTED && !firstScatter) {
                    Vector3f posVec = r.o - Point3f(0.f, 0.f, 0.f);
                    Float alpha = SphericalTheta(Normalize(posVec));
                    Float theta = SphericalTheta(r.d);
                    Float phi = SphericalPhi(r.d);

                    sphere_pdf[radius.x].pdf.populatePdf(alpha, theta, phi);
                }
            }
            //}
        },
        Point2i(count, 1));
    ParallelFor2D(
        [&](Point2i radius) {

		sphere_pdf[radius.x].pdf.PrecomputeAlphaPdf();
		sphere_pdf[radius.x].pdf.PrecomputeThetaPdf();
		sphere_pdf[radius.x].pdf.PrecomputePhiPdf();

	},
        Point2i(count, 1));
    /*vector<float> val = sphere_pdf[0].pdf.generatedPdf1D();
    vector<float> val2 = sphere_pdf[0].pdf.generatedPdf1D(34);
    vector<float> val3 = sphere_pdf[0].pdf.generatedPdf1D(33, 43);*/
}

int HomogeneousMedium::PopulatePdf(float minRadius, float maxRadius,
                                    float stepSize, int binSize) {
    int entries = 0;

	float step = stepSize;
    float currRadius = minRadius;
    while (currRadius < maxRadius) {
        ShellPdf shellPdf(binSize);
        SpherePdf spherePdf = {currRadius, shellPdf};
        sphere_pdf.push_back(spherePdf);

        currRadius += std::exp(step);
        step += stepSize;
        entries++;
    }
    return entries;
}
}  // namespace pbrt

//
// while (currRadius < maxRadius) {  // we need to incorporate the sample size
//                                  // here as well ....
//                                  // variable initializations
//    ShellPdf shellPdf(binSize);
//    float *stratifiedSamples = new float[sampleCount]();
//    // StratifiedSample1D(stratifiedSamples, sampleCount, RNG(1), true);
//    sampler.StartPixel(Point2i(int(scopeRandomFloats(randGen) * 1024.f),
//                               int(scopeRandomFloats(randGen) * 1024.f)));
//
//    for (int i = 0; i < sampleCount; i++) {
//        Ray r(Point3f(0.f, 0.f, 0.f),
//              Vector3f(0.f, 1.f, 0.f));  // initialized a ray pointing upwards
//        RAY_STATE state = RAY_SCATTERED;
//        bool firstScatter = true;
//        while (state == RAY_SCATTERED) {
//            float t_exit = computeRaySphereIntersect(Point3f(0.f, 0.f, 0.f),
//                                                     currRadius, r);
//
//            int channel =
//                std::min((int)(scopeRandomFloats(randGen) *
//                               Spectrum::nSamples),  // sampler is not working
//                                                     // -> check this
//                         Spectrum::nSamples - 1);
//            float t_scatter = -std::log(scopeRandomFloats(randGen)) /
//                              sigma_s[channel];  // 1 - sampler.Get1D()
//            float t_absorb = -std::log(scopeRandomFloats(randGen)) /
//                             sigma_a[channel];  // 1 - sampler.Get1D()
//
//            if (min(t_absorb, t_scatter) >=
//                t_exit) {  //  || Distance(Point3f(0.f, 0.f, 0.f), r.o) <
//                           //  currRadius
//                r = Ray(r(t_exit), r.d);
//                state = RAY_EXITTED;
//            } else {
//                firstScatter = false;
//                if (t_absorb < t_scatter)
//                    state = RAY_ABSORBED;
//                else {
//                    Vector3f wo = -r.d, wi;
//                    phase->Sample_p(wo, &wi,
//                                    Point2f(scopeRandomFloats(randGen),
//                                            scopeRandomFloats(randGen)));
//                    Point3f p = r.o + (r.d * t_scatter);
//                    r = Ray(p, wi);
//                }
//            }
//        }
//
//        if (state == RAY_EXITTED && !firstScatter) {
//            Vector3f posVec = r.o - Point3f(0.f, 0.f, 0.f);
//            Float alpha = SphericalTheta(Normalize(posVec));
//            Float theta = SphericalTheta(r.d);
//            Float phi = SphericalPhi(r.d);
//
//            shellPdf.populatePdf(alpha, theta, phi);
//        }
//    }
//
//    // don't know what arena alloc is, but check it if this does not work
//    // Vector3f wo = -r.d,
//    //         wi;  // check this why the direction of the ray is negative
//    //         ...
//    // phase->Sample_p(wo, &wi, sampler.Get2D());
//
//    // implement the part where you test the time for diffuse, scatter or
//    // exit populate the pdf move on to the next steps add check to make
//    // sure the distance from the origin and the point does not exceed the
//    // radius
//    // r(t_exit)
//
//    // increment the radius size
//
//    SpherePdf spherePdf = {currRadius, shellPdf};
//    sphere_pdf.push_back(spherePdf);
//    currRadius += stepSize;
//}