
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

// integrators/volpath.cpp*
#include "integrators/volpath.h"
#include "bssrdf.h"
#include "camera.h"
#include "film.h"
#include "interaction.h"
#include "paramset.h"
#include "scene.h"
#include "stats.h"

namespace pbrt {

STAT_INT_DISTRIBUTION("Integrator/Path length", pathLength);
STAT_COUNTER("Integrator/Volume interactions", volumeInteractions);
STAT_COUNTER("Integrator/Surface interactions", surfaceInteractions);

// VolPathIntegrator Method Definitions
void VolPathIntegrator::Preprocess(const Scene &scene, Sampler &sampler) {
    lightDistribution =
        CreateLightSampleDistribution(lightSampleStrategy, scene);
}

Spectrum VolPathIntegrator::Li(const RayDifferential &r, const Scene &scene,
                               Sampler &sampler, MemoryArena &arena,
                               int depth) const {
    ProfilePhase p(Prof::SamplerIntegratorLi);
    Spectrum L(0.f), beta(1.f);
    RayDifferential ray(r);
    bool specularBounce = false;
    int bounces;
    // Added after book publication: etaScale tracks the accumulated effect
    // of radiance scaling due to rays passing through refractive
    // boundaries (see the derivation on p. 527 of the third edition). We
    // track this value in order to remove it from beta when we apply
    // Russian roulette; this is worthwhile, since it lets us sometimes
    // avoid terminating refracted rays that are about to be refracted back
    // out of a medium and thus have their beta value increased.
    Float etaScale = 1;

    for (bounces = 0;; ++bounces) {
        // Intersect _ray_ with scene and store intersection in _isect_
        SurfaceInteraction isect;
        bool foundIntersection = scene.Intersect(ray, &isect);

        // shift this inside the medium check code
        // compile the pdf to do the 1D inverse method
        //
        // Code for looping through the SHELL PDF
        
        Vector3f rayOriginOffset;
        bool isRayOffset = false;
        
        float boundaryDistance = 0.f;
        const SpherePdf *spherePdf =
            ray.medium ? ray.medium->getShellPdf(boundaryDistance) : NULL;

        if (ray.medium != 0 && spherePdf &&
                    spherePdf->pdf.SampleCount()) {  // 
            
			// compute the distance from boundary once 
			boundaryDistance = ComputeNearestDistance(ray, scene);
            spherePdf = ray.medium->getShellPdf(boundaryDistance);

            // check this value, do not use hard coded one here
            if (boundaryDistance > 3.f) {
                Vector3f wo = -ray.d;
                auto shellPdf = ray.medium->getShellPdf(boundaryDistance);
                Transform rayTransform = ComputeShellRayTransform(wo);
                Float prob;

                // the bin index of the angles
                int alphaBinIndex = SampleFromDistribution(
                    shellPdf->pdf.GetAlphaPdf(), &prob, sampler);

                // the angles for all the values
                float alphaAngle = shellPdf->pdf.AlphaBinValue(alphaBinIndex);
                float betaAngle = 2 * Pi * sampler.Get1D();

                // compute the movement direction of the ray...
                rayOriginOffset = SphericalDirection(
                    std::sin(alphaAngle), std::cos(alphaAngle), betaAngle);
                rayOriginOffset = shellPdf->radius * Normalize(rayTransform(rayOriginOffset));
                //printf("The radius value is %f\n", shellPdf->radius);
                isRayOffset = true;
			}
		}


        // Sample the participating medium, if present
        MediumInteraction mi; // check if the tMax is the problem in the medium sample func
        if (ray.medium) beta *= ray.medium->Sample(ray, sampler, arena, &mi, isRayOffset ? &rayOriginOffset : NULL);
        if (beta.IsBlack()) break;

        // Handle an interaction with a medium or a surface
        if (mi.IsValid()) {
            // Terminate path if ray escaped or _maxDepth_ was reached
            if (bounces >= maxDepth) break;

            ++volumeInteractions;
            // Handle scattering at point in medium for volumetric path tracer
            const Distribution1D *lightDistrib =
                lightDistribution->Lookup(mi.p);
            L += beta * UniformSampleOneLight(mi, scene, arena, sampler, true,
                                              lightDistrib);

            // Code for looping through the SHELL PDF
            if (spherePdf && spherePdf
                    ->pdf.SampleCount() != 0) {
                //float boundaryDistance = ComputeNearestDistance(ray, scene);
                // check this value, do not use hard coded one here
                if (boundaryDistance > 3.f) {
                    Vector3f wo = -ray.d;
                    auto shellPdf =
                        mi.GetMedium()->getShellPdf(boundaryDistance);
                    Transform rayTransform = ComputeShellRayTransform(wo);
                    Float prob;

                    // the bin index of the angles
                    int alphaBinIndex = SampleFromDistribution(
                        shellPdf->pdf.GetAlphaPdf(), &prob, sampler);
                    //generatedPdf1D(shellPdf->pdf.getPdf());
                    
                    int thetaBinIndex = SampleFromDistribution(
                        shellPdf->pdf.GetThetaPdf(alphaBinIndex), &prob,
                        sampler);
                    int phiBinIndex = SampleFromDistribution(
                        shellPdf->pdf.GetPhiPdf(alphaBinIndex,
                                                    thetaBinIndex),
                        &prob, sampler);

                    // the angles for all the values
                    /*float alphaAngle =
                        shellPdf.pdf.AlphaBinValue(alphaBinIndex);
                    float betaAngle = 2 * Pi * sampler.Get1D();*/
                    float thetaAngle =
                        shellPdf->pdf.ThetaBinValue(thetaBinIndex);
                    float phiAngle = shellPdf->pdf.PhiBinValue(phiBinIndex);

                    // this is in the original // remember the transformation
                    // application
                    Vector3f rayDir = SphericalDirection(
                        std::sin(thetaAngle), std::cos(thetaAngle), phiAngle);
                    rayDir = boundaryDistance * rayTransform(rayDir);

                     //compute the movement direction of the ray...
                    /*Vector3f rayOriginOffset = SphericalDirection(
                        std::sin(alphaAngle), std::cos(alphaAngle), betaAngle);
                    rayOriginOffset = rayTransform(rayOriginOffset);
                    Point3f rayOrigin = ray.o + rayOriginOffset;*/

					ray = mi.SpawnRay(rayDir);

					//Vector3f wi;
     //               mi.phase->Sample_p(wo, &wi, sampler.Get2D());
     //               ray = mi.SpawnRay(wi);

                    // think about the ray spawning direction... (this is really
                    // important)

                    /*Distribution1D alphaDistrib(&alphaPdf[0],
                    alphaPdf.size()); Float alphaVal; int alphaNum =
                    alphaDistrib.SampleDiscrete(sampler.Get1D(), &alphaVal);

                                        vector<float> thetaPdf =
                    shellPdf.pdf.generatedPdf1D(alphaNum); Float alphaVal;*/

                    // shellPdf.pdf.generatedPdf1D();
                } else {
                    Vector3f wo = -ray.d, wi;
                    mi.phase->Sample_p(wo, &wi, sampler.Get2D());
                    ray = mi.SpawnRay(wi);
                }
            } else {
                Vector3f wo = -ray.d, wi;
                mi.phase->Sample_p(wo, &wi, sampler.Get2D());
                ray = mi.SpawnRay(wi);
            }

			           /*Vector3f wo = -ray.d, wi;
                       mi.phase->Sample_p(wo, &wi, sampler.Get2D());
                       ray = mi.SpawnRay(wi);
*/
            specularBounce = false;
        } else {
            ++surfaceInteractions;
            // Handle scattering at point on surface for volumetric path tracer

            // Possibly add emitted light at intersection
            if (bounces == 0 || specularBounce) {
                // Add emitted light at path vertex or from the environment
                if (foundIntersection)
                    L += beta * isect.Le(-ray.d);
                else
                    for (const auto &light : scene.infiniteLights)
                        L += beta * light->Le(ray);
            }

            // Terminate path if ray escaped or _maxDepth_ was reached
            if (!foundIntersection || bounces >= maxDepth) break;

            // Compute scattering functions and skip over medium boundaries
            isect.ComputeScatteringFunctions(ray, arena, true);
            if (!isect.bsdf) {
                ray = isect.SpawnRay(ray.d);
                bounces--;
                continue;
            }

            // Sample illumination from lights to find attenuated path
            // contribution
            const Distribution1D *lightDistrib =
                lightDistribution->Lookup(isect.p);
            L += beta * UniformSampleOneLight(isect, scene, arena, sampler,
                                              true, lightDistrib);

            // Sample BSDF to get new path direction
            Vector3f wo = -ray.d, wi;
            Float pdf;
            BxDFType flags;
            Spectrum f = isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdf,
                                              BSDF_ALL, &flags);
            if (f.IsBlack() || pdf == 0.f) break;
            beta *= f * AbsDot(wi, isect.shading.n) / pdf;
            DCHECK(std::isinf(beta.y()) == false);
            specularBounce = (flags & BSDF_SPECULAR) != 0;
            if ((flags & BSDF_SPECULAR) && (flags & BSDF_TRANSMISSION)) {
                Float eta = isect.bsdf->eta;
                // Update the term that tracks radiance scaling for refraction
                // depending on whether the ray is entering or leaving the
                // medium.
                etaScale *=
                    (Dot(wo, isect.n) > 0) ? (eta * eta) : 1 / (eta * eta);
            }
            ray = isect.SpawnRay(wi);

            // Account for attenuated subsurface scattering, if applicable
            if (isect.bssrdf && (flags & BSDF_TRANSMISSION)) {
                // Importance sample the BSSRDF
                SurfaceInteraction pi;
                Spectrum S = isect.bssrdf->Sample_S(
                    scene, sampler.Get1D(), sampler.Get2D(), arena, &pi, &pdf);
                DCHECK(std::isinf(beta.y()) == false);
                if (S.IsBlack() || pdf == 0) break;
                beta *= S / pdf;

                // Account for the attenuated direct subsurface scattering
                // component
                L += beta *
                     UniformSampleOneLight(pi, scene, arena, sampler, true,
                                           lightDistribution->Lookup(pi.p));

                // Account for the indirect subsurface scattering component
                Spectrum f = pi.bsdf->Sample_f(pi.wo, &wi, sampler.Get2D(),
                                               &pdf, BSDF_ALL, &flags);
                if (f.IsBlack() || pdf == 0) break;
                beta *= f * AbsDot(wi, pi.shading.n) / pdf;
                DCHECK(std::isinf(beta.y()) == false);
                specularBounce = (flags & BSDF_SPECULAR) != 0;
                ray = pi.SpawnRay(wi);
            }
        }

        // Possibly terminate the path with Russian roulette
        // Factor out radiance scaling due to refraction in rrBeta.
        Spectrum rrBeta = beta * etaScale;
        if (rrBeta.MaxComponentValue() < rrThreshold && bounces > 3) {
            Float q = std::max((Float).05, 1 - rrBeta.MaxComponentValue());
            if (sampler.Get1D() < q) break;
            beta /= 1 - q;
            DCHECK(std::isinf(beta.y()) == false);
        }
    }
    ReportValue(pathLength, bounces);
    return L;
}

float VolPathIntegrator::ComputeNearestDistance(const RayDifferential &ray,
                                                const Scene &scene) const {
    if (ray.medium) {
        float minRadius = 10000000.f;
        Point3f origin = ray.o;

        for (Point3f vertx : platonicSolidVertx) {
            SurfaceInteraction ssMeshIsect;
            Vector3f dir = Normalize(vertx - Point3f(0.f, 0.f, 0.f));
            RayDifferential rayClone(origin, dir, Infinity, 0.f, ray.medium);

            scene.Intersect(rayClone, &ssMeshIsect);

            float dist = Distance(origin, ssMeshIsect.p);
            if (minRadius > dist) minRadius = dist;
        }
        return minRadius;
    }
    return -1.f;
}

Transform VolPathIntegrator::ComputeShellRayTransform(Vector3f wo) const {
    Float theta = Degrees(SphericalTheta(wo));
    Float phi = Degrees(SphericalPhi(wo));

    Transform thetaTransform(Rotate(theta, Vector3f(0.f, 1.f, 0.f))
                                 .GetMatrix());  // 1, 0, 0 // 0, 1, 0
    Transform phiTransform(Rotate(phi, Vector3f(0.f, 0.f, 1.f))
                               .GetMatrix());  // 0, 0, 1 // 0, 0, 1

    // Transform ve = tTransform * pTransform;
    Transform phiThetaTransform(
        (phiTransform * thetaTransform));  // .GetInverseMatrix()
    /*Transform s(ve.GetInverseMatrix());
    Transform s1(ev.GetInverseMatrix());*/

    /*
    Vector3f qrs = ve(wo);
    Vector3f qrt = ev(wo);
    Vector3f rrs = s(wo);*/
    // Vector3f topVec = phiThetaTransform(wo);  // this one
    return phiThetaTransform;
}

int VolPathIntegrator::SampleFromDistribution(const vector<float> *pdf, Float *prob,
                                              Sampler &sampler) const {
    Distribution1D distrib(&(*pdf)[0], pdf->size());
    return distrib.SampleDiscrete(sampler.Get1D(), prob);
}

VolPathIntegrator *CreateVolPathIntegrator(
    const ParamSet &params, std::shared_ptr<Sampler> sampler,
    std::shared_ptr<const Camera> camera) {
    int maxDepth = params.FindOneInt("maxdepth", 5);
    int np;
    const int *pb = params.FindInt("pixelbounds", &np);
    Bounds2i pixelBounds = camera->film->GetSampleBounds();
    if (pb) {
        if (np != 4)
            Error("Expected four values for \"pixelbounds\" parameter. Got %d.",
                  np);
        else {
            pixelBounds = Intersect(pixelBounds,
                                    Bounds2i{{pb[0], pb[2]}, {pb[1], pb[3]}});
            if (pixelBounds.Area() == 0)
                Error("Degenerate \"pixelbounds\" specified.");
        }
    }
    Float rrThreshold = params.FindOneFloat("rrthreshold", 1.);
    std::string lightStrategy =
        params.FindOneString("lightsamplestrategy", "spatial");
    return new VolPathIntegrator(maxDepth, camera, sampler, pixelBounds,
                                 rrThreshold, lightStrategy);
}

}  // namespace pbrt
