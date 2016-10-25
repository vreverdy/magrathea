/* ******************************* RAYTRACER ******************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        RAYTRACING
// TITLE :          Raytracer
// DESCRIPTION :    Main function of the raytracer
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           raytracer.cpp
/// \brief          Main function of the raytracer
/// \author         Vincent Reverdy (vince.rev@gmail.com)
/// \date           2012-2013
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/



// ------------------------------ PREPROCESSOR ------------------------------ //
// Include C++
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <utility>
#include <random>
#include <string>
#include <vector>
#include <deque>
#include <array>
#include <tuple>
#include <atomic>
#include <mutex>
// Include libs
#include <mpi.h>
// Include project
#include "../magrathea/simplehyperoctree.h"
#include "../magrathea/simplehyperoctreeindex.h"
#include "../magrathea/hypersphere.h"
#include "../magrathea/hypercube.h"
#include "../magrathea/constants.h"
#include "../magrathea/evolution.h"
#include "../magrathea/timer.h"
#include "cone.h"
#include "utility.h"
#include "gravity.h"
#include "input.h"
#include "output.h"
#include "integrator.h"
// Misc
using namespace magrathea;
// -------------------------------------------------------------------------- //



// ---------------------------------- MAIN ---------------------------------- //
// Main function
/// \brief          Main function.
/// \details        Main function of the program. It takes a namelist as a 
///                 parameter and run the raytracing using MPI parallelization.
///                 To compile it, uses: 
///                 mpic++ -Wall -Wextra -pedantic -std=c++11 -O3 raytracer.cpp 
///                 -o raytracer -lpthread
///                 To run it, uses:
///                 ./raytracer raytracer.txt
/// \param[in]      argc Number of arguments.
/// \param[in]      argv List of arguments.
/// \return         Zero on success, error code otherwise.
int main(int argc, char* argv[])
{
    // Constants
    /*using integer = int;
    using uint = unsigned int;
    using real = double;
    using floating = float;
    using point = std::array<real, 3>;
    using position = std::ratio<0>;
    using extent = std::ratio<1>;
    using indexing = unsigned long long int;//__uint128_t;
    using element = std::pair<SimpleHyperOctreeIndex<indexing, 3>, Gravity<floating, 3> >;
    using evolution = Evolution<Photon<real, 3> >;
    static constexpr uint zero = 0;
    static constexpr uint one = 1;
    static constexpr uint two = 2;
    static constexpr uint three = 3;
    static constexpr uint dimension = 3;
    static constexpr real pi = Constants<real>::pi();
    static constexpr real msun = Constants<real>::msun();
    static constexpr real g = Constants<real>::g();
    static constexpr real rposition = static_cast<real>(position::num)/static_cast<real>(position::den);
    static constexpr point center({{rposition, rposition, rposition}});
    static constexpr real diameter = static_cast<real>(extent::num)/static_cast<real>(extent::den);
    static constexpr uint digits = std::numeric_limits<real>::max_digits10;
    static constexpr char dot = '.';
    static constexpr char dotc = 'd';
    static const std::string all = "all";
    static const std::string namelist = argc > 1 ? std::string(argv[1]) : std::string("raytracer.txt");
    
    // Parameters
    std::map<std::string, std::string> parameter = Input::parse(namelist);
    const uint propagation = std::stoul(parameter["propagation"]);
    const uint preparation = std::stoul(parameter["preparation"]);
    const uint visualization = std::stoul(parameter["visualization"]);
    const uint homogeneous = std::stoul(parameter["homogeneous"]);
    const uint schwarzschild = std::stoul(parameter["schwarzschild"]);
    const uint tdependent = std::stoul(parameter["tdependent"]);
    const uint test = std::stoul(parameter["test"]);
    const uint seed = std::stoul(parameter["seed"]);
    const uint allocation = std::stoul(parameter["allocation"]);
    const uint alphacoeff = std::stoul(parameter["alphacoeff"]);
    const uint microcoeff = std::stoul(parameter["microcoeff"]);
    const real mpc = std::stod(parameter["mpc"]);
    const real rhoch2 = std::stod(parameter["rhoch2"]);
    const std::string paramfile = parameter["paramfile"];
    const std::string evolfile = parameter["evolfile"];
    const std::string cubedir = parameter["cubedir"];
    const std::string cubefmt = parameter["cubefmt"];
    const std::string conedir = parameter["conedir"];
    const std::string conefmt = parameter["conefmt"];
    const std::string outputdir = parameter["outputdir"];
    const std::string outputsep = parameter["outputsep"];
    const std::string outputint = parameter["outputint"];
    const std::string outputopening = parameter["outputopening"];
    const std::string outputprefix = parameter["outputprefix"];
    const std::string outputsuffix = parameter["outputsuffix"];
    const std::string outputvisu = parameter["outputvisu"];
    const std::string outputhomo = parameter["outputhomo"];
    const std::string outputschw = parameter["outputschw"];
    const std::string outputtdep = parameter["outputtdep"];
    const std::string outputtest = parameter["outputtest"];
    const std::string outputrk4 = parameter["outputrk4"];
    const std::string outputeul = parameter["outputeul"];
    const std::string outputngp = parameter["outputngp"];
    const std::string outputcic = parameter["outputcic"];
    const std::string outputref = parameter["outputref"];
    const std::string outputtree = parameter["outputtree"];
    const std::string outputstat = parameter["outputstat"];
    const uint correction = std::stoul(parameter["correction"]);
    const uint coarsecorrection = std::stoul(parameter["coarsecorrection"]);
    const uint acorrection = std::stoul(parameter["acorrection"]);
    const uint nreference = std::stoul(parameter["nreference"]);
    const uint ncoarse = std::stoul(parameter["ncoarse"]);
    const uint ncones = std::stoul(parameter["ncones"]);
    const uint ntrajectories = std::stoul(parameter["ntrajectories"]);
    const uint nsteps = std::stoul(parameter["nsteps"]);
    const uint nstat = std::stoul(parameter["nstat"]);
    const uint nbundlemin = std::stoul(parameter["nbundlemin"]);
    const uint nbundlecnt = std::stoul(parameter["nbundlecnt"]);
    const real openingmin = std::stod(parameter["openingmin"]);
    const uint openingcnt = std::stoul(parameter["openingcnt"]);
    const std::string interpolation = parameter["interpolation"];
    const std::string statistic = parameter["statistic"];
    const uint makestat = std::stoul(parameter["makestat"]);
    const integer savemode = std::stol(parameter["savemode"]);
    const uint savetree = std::stoul(parameter["savetree"]);
    const real lboxmpch0 = std::stod(parameter["lboxmpch0"]);
    const real lboxmpc0 = std::stod(parameter["lboxmpc0"]);
    const real x0 = std::stod(parameter["x0"]);
    const real y0 = std::stod(parameter["y0"]);
    const real z0 = std::stod(parameter["z0"]);
    const real x = std::stod(parameter["x"]);
    const real y = std::stod(parameter["y"]);
    const real z = std::stod(parameter["z"]);
    const real rvismin = std::stod(parameter["rvismin"]);
    const real rvismax = std::stod(parameter["rvismax"]);
    const real vislim = std::stod(parameter["vislim"]);
    const real massmsun = std::stod(parameter["massmsun"]);
    const uint levelmin = std::stoul(parameter["levelmin"]);
    const uint levelmax = std::stoul(parameter["levelmax"]);
    const real criterion = std::stod(parameter["criterion"]);
    const uint stepmin = std::stoul(parameter["stepmin"]);
    const uint stepmul = std::stoul(parameter["stepmul"]);
    const uint stepinc = std::stoul(parameter["stepinc"]);
    const uint stepcnt = std::stoul(parameter["stepcnt"]);
    const real yshiftmin = std::stod(parameter["yshiftmin"]);
    const real yshiftmul = std::stod(parameter["yshiftmul"]);
    const real yshiftinc = std::stod(parameter["yshiftinc"]);
    const real yshiftcnt = std::stoul(parameter["yshiftcnt"]);
    
    // Initialization
    std::vector<std::string> interpolations = (interpolation == all) ? (std::vector<std::string>({"redshift", "t", "a", "r"})) : (std::vector<std::string>({interpolation}));
    std::vector<std::string> statistics = (statistic == all) ? (std::vector<std::string>({"distance", "distance2", "homogeneous", "inhomogeneous"})) : (std::vector<std::string>({statistic}));
    uint interpcnt = interpolations.size();
    uint statcnt = statistics.size();
    uint interpcase = zero;
    uint statcase = zero;
    std::string interp;
    std::string stat;
    Timer<real> timer;
    FileList conefile(conefmt, zero, ncones, zero, conedir);*/
    //SimpleHyperOctree<real, SimpleHyperOctreeIndex<uint, dimension>, std::string, dimension, position, extent> filetree;
    //SimpleHyperOctree<real, SimpleHyperOctreeIndex<indexing, dimension>, Gravity<floating, dimension>, dimension, position, extent> octree;
    //SimpleHyperOctree<real, SimpleHyperOctreeIndex<indexing, dimension>, Gravity<floating, dimension>, dimension, position, extent> homotree;
    /*std::mt19937 engine(seed > zero ? seed : std::random_device()());
    std::uniform_real_distribution<real> distribution(zero, one);
    HyperSphere<dimension, point> sphere(center, diameter/two);
    HyperCube<dimension, point> cube(center, diameter);
    HyperSphere<dimension, point> microsphere(center, diameter/microcoeff);
    std::vector<std::string> filelist;
    std::vector<point> tiling(ncones);
    std::vector<Cone<point> > cone(ncones);
    std::vector<evolution> trajectory(ntrajectories);
    std::vector<Photon<real, dimension> > photons(ntrajectories);
    std::vector<real> random(ntrajectories);
    std::vector<SimpleHyperOctreeIndex<indexing, dimension> > index;
    std::deque<std::atomic<uint> > count;
    std::array<std::vector<real>, two+two> cosmology;
    Cone<point> conic;
    Photon<real, dimension> photon;
    Evolution<Photon<real, dimension> > reference;
    real alpha = (two*pi)/alphacoeff;
    real h = zero;
    real omegam = zero;
    real lboxmpch = zero;
    real amin = zero;
    uint nfiles = zero;
    uint nbundle = zero;
    uint size = zero;
    real opening = zero;
    integer nthreads = std::thread::hardware_concurrency();
    integer ntasks = nthreads*zero;
    integer rank = zero;
    std::ofstream stream;
    std::string filename;
    std::mutex mutex;
    uint statmod = zero;
    uint statlength = zero;
    uint statsize = zero;
    uint statgsize = zero;
    std::vector<real> statrefx;
    std::vector<real> statrefy;
    std::vector<real> statrefz;
    std::vector<std::vector<real> > statx;
    std::vector<std::vector<real> > staty;
    std::vector<real> statmean;
    std::vector<real> statstd;
    std::vector<real> statgmean;
    std::vector<real> statgstd;
    
    // Message passing interface
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Generate cones
    alpha = std::sqrt(pi)*std::asin(sphere.uniform<dimension-one>(std::begin(tiling), std::end(tiling)).first/sphere.diameter());
    Utility::parallelize(ncones, [=, &tiling, &cone](const uint i){cone[i].assign(sphere.position(), tiling[i], alpha);});
    
    // Read cosmology
    cosmology = Input::acquire(cubedir, paramfile, evolfile, h, omegam, lboxmpch, mpc);
    Input::homogenize(octree.assign(nreference, zero));
    reference.append(Integrator::launch(center[zero], center[one], center[two], center[zero]+diameter/two, center[one], center[two]));
    Integrator::integrate<-1>(reference, cosmology, octree, lboxmpch*mpc/h, std::pow(two, static_cast<uint>(std::log2(std::get<0>(cosmology).size()/std::pow(two, nreference)+one)+one)+one));
    cosmology = Input::correct(cosmology, reference);
    octree.clear();
    
    // Construct octree
    if (propagation || visualization || test) {
        Input::load(octree, conefile[rank]);
        Input::correct(octree, correction, coarsecorrection, acorrection, amin);
        octree.shrink();
        if (propagation || visualization || test) {
            Input::sistemize(octree, h, omegam, lboxmpch, mpc, rhoch2);
        } else if (homogeneous) {
            Input::homogenize(octree);
        }
    } else if (preparation) {
        octree.reserve(allocation);
        Input::filetree(filetree, cubedir, cubefmt);
        for (uint icone = zero; icone < ncones; ++icone) {
            if (icone%static_cast<uint>(ntasks) == static_cast<uint>(rank)) {
                octree.clear();
                filelist.clear();
                conic = cone[icone];
                Input::prepare(filelist, filetree, microsphere, conic);
                std::shuffle(std::begin(filelist), std::end(filelist), engine);
                nfiles = filelist.size();
                for (uint ifile = zero; ifile < nfiles; ++ifile) {
                    Input::import(octree, filelist[ifile], ncoarse, [=, &octree, &microsphere, &conic](const element& e){return Input::collide(octree, std::get<0>(e), microsphere, conic);});
                }
                Input::save(octree, conefile[icone]);
            }
        }
    } else if (homogeneous) {
        octree.assign(levelmin, zero);
        Input::homogenize(octree);
    } else if (schwarzschild) {
        octree.assign(levelmin, zero);
        Input::schwarzschildify(octree, center, massmsun*msun, lboxmpc0*mpc, [=, &massmsun, &lboxmpc0, &mpc, &levelmax, &criterion](const Gravity<floating, dimension>& gravity, const uint level){return level < levelmax && std::abs(gravity.phi()) > std::abs((g*massmsun*msun)/(lboxmpc0*mpc))*std::pow(criterion, level);});
    }
    octree.shrink();
    octree.update();
    
    // Execution
    if (propagation) {
        // Initialization
        homotree.assign(ncoarse/two, zero);
        photon = Integrator::launch(center[zero], center[one], center[two], center[zero]+diameter/two, center[one], center[two]);    
        for (uint itrajectory = zero; itrajectory < ntrajectories; ++itrajectory) {
            photons[itrajectory] = Integrator::launch(microsphere, cone[rank], cone, engine, distribution);
            random[itrajectory] = distribution(engine)*two*pi;
        }
        // Loop over configurations
        for (uint ibundle = zero; ibundle < nbundlecnt; ++ibundle) {
            nbundle = nbundlemin+ibundle;
            for (uint iopening = zero; iopening < openingcnt; ++iopening) {
                opening = (iopening+one)*openingmin;
                for (uint iinterp = zero; iinterp < interpcnt; ++iinterp) {
                    interp = interpolations[iinterp];
                    for (uint istat = zero; istat < statcnt; ++istat) {
                        stat = statistics[istat];
                        // Reference
                        filename = Output::name(outputdir, outputprefix, outputsep, std::make_pair(outputint, nbundle), outputsep, std::make_pair(outputopening, opening), outputsep, interp, outputsep, std::make_pair(outputint, rank));
                        std::replace(filename.begin(), filename.end(), dot, dotc);
                        reference = Integrator::propagate<-1>(photon, nbundle, opening, real(), interp, cosmology, homotree, lboxmpch*mpc/h, nsteps*(one << (ncoarse-ncoarse/two))*two, real(), std::signbit(savemode) ? Output::name() : Output::name(filename, outputsuffix));
                        // Integration without statistics
                        if (makestat == zero) {
                            Utility::parallelize(ntrajectories, [=, &photons, &nbundle, &opening, &random, &interp, &cosmology, &octree, &lboxmpch, &mpc, &h, &nsteps, &amin, &filename, &reference, &savemode, &outputsep, &outputint](const uint i){Integrator::propagate(photons[i], nbundle, opening, random[i], interp, cosmology, octree, lboxmpch*mpc/h, nsteps, amin, std::signbit(savemode) ? Output::name() : Output::name(savemode ? Output::name(filename, outputsep, std::make_pair(outputint, i), outputsep, outputint) : Output::name(filename, outputsep, std::make_pair(outputint, i), outputsep, std::make_pair(outputint, zero)), outputsuffix), reference);});
                        // Integration with statistics
                        } else {
                            // Clear statistics arrays
                            statrefx.clear();
                            statrefy.clear();
                            statrefz.clear();
                            statx.clear();
                            staty.clear();
                            statmean.clear();
                            statstd.clear();
                            statgmean.clear();
                            statgstd.clear();
                            // Interpolation case
                            if (interp == "redshift") {
                                interpcase = zero;
                            } else if (interp == "a") {
                                interpcase = one;
                            } else if (interp == "t") {
                                interpcase = two;
                            } else if (interp == "r") {
                                interpcase = three;
                            } else {
                                interpcase = zero;
                            }
                            // Statistics case
                            if (stat == "distance") {
                                statcase = zero;
                            } else if (stat == "distance2") {
                                statcase = one;
                            } else if (stat == "homogeneous") {
                                statcase = two;
                            } else if (stat == "inhomogeneous") {
                                statcase = three;
                            } else {
                                statcase = zero;
                            }
                            // Integration
                            Utility::parallelize(ntrajectories, [=, &photons, &nbundle, &opening, &random, &interp, &cosmology, &octree, &lboxmpch, &mpc, &h, &nsteps, &amin, &filename, &reference, &savemode, &outputsep, &outputint, &mutex, &interpcase, &statcase, &statx, &staty](const uint i){
                                evolution result = Integrator::propagate(photons[i], nbundle, opening, random[i], interp, cosmology, octree, lboxmpch*mpc/h, nsteps, amin, std::signbit(savemode) ? Output::name() : Output::name(savemode ? Output::name(filename, outputsep, std::make_pair(outputint, i), outputsep, outputint) : Output::name(filename, outputsep, std::make_pair(outputint, i), outputsep, std::make_pair(outputint, zero)), outputsuffix), reference);
                                std::vector<std::vector<real> > tmp(two, std::vector<real>(result.size()));
                                for (uint j = zero; j < result.size(); ++j) {
                                    if (interpcase == zero) {
                                        tmp[zero][j] = result[j].redshift();
                                    } else if (interpcase == one) {
                                        tmp[zero][j] = result[j].t();
                                    } else if (interpcase == two) {
                                        tmp[zero][j] = result[j].a();
                                    } else if (interpcase == three) {
                                        tmp[zero][j] = std::sqrt(std::pow(result[j].x()-result[zero].x(), two)+std::pow(result[j].y()-result[zero].y(), two)+std::pow(result[j].z()-result[zero].z(), two));
                                    }
                                    if (statcase == zero) {
                                        tmp[one][j] = result[j].distance();
                                    } else if (statcase == one) {
                                        tmp[one][j] = result[j].distance()*result[j].distance();
                                    } else if (statcase == two) {
                                        tmp[one][j] = result[zero].a()/result[j].a()-one;
                                    } else if (statcase == three) {
                                        tmp[one][j] = result[j].redshift();
                                    }
                                }
                                if (result.size() > zero) {
                                    mutex.lock();
                                    statx.emplace_back(std::move(tmp[zero]));
                                    staty.emplace_back(std::move(tmp[one]));
                                    mutex.unlock();
                                }
                            });
                            // Transfer reference
                            reference.container().erase(std::remove_if(reference.container().begin(), reference.container().end(), [=, &amin](const Photon<real, dimension>& p){return std::isnormal(amin) && (p.a() < amin);}), reference.container().end());
                            statmod = std::max(one, static_cast<uint>(reference.size())/std::max(nstat, one));
                            for (uint j = zero; j < reference.size(); ++j) {
                                if (j%statmod == zero) {
                                    if (interpcase == zero) {
                                        statrefx.emplace_back(reference[j].redshift());
                                    } else if (interpcase == one) {
                                        statrefx.emplace_back(reference[j].t());
                                    } else if (interpcase == two) {
                                        statrefx.emplace_back(reference[j].a());
                                    } else if (interpcase == three) {
                                        statrefx.emplace_back(std::sqrt(std::pow(reference[j].x()-reference[zero].x(), two)+std::pow(reference[j].y()-reference[zero].y(), two)+std::pow(reference[j].z()-reference[zero].z(), two)));
                                    }
                                    if (statcase == zero) {
                                       statrefy.emplace_back(reference[j].distance());
                                    } else if (statcase == one) {
                                        statrefy.emplace_back(reference[j].distance()*reference[j].distance());
                                    } else if (statcase == two) {
                                        statrefy.emplace_back(reference[zero].a()/reference[j].a()-one);
                                    } else if (statcase == three) {
                                        statrefy.emplace_back(reference[j].redshift());
                                    }
                                    statrefz.emplace_back(reference[j].redshift());
                                }
                            }
                            // Resize
                            statlength = statrefx.size();
                            statsize = statx.size();
                            statmean.resize(statlength, real());
                            statstd.resize(statlength, real());
                            statgmean.resize(statlength, real());
                            statgstd.resize(statlength, real());
                            statmean.shrink_to_fit();
                            statstd.shrink_to_fit();
                            statgmean.shrink_to_fit();
                            statgstd.shrink_to_fit();
                            // Reinterpolate
                            Utility::parallelize(statsize, [=, &statrefx, &statx, &staty](const uint j){staty[j] = Utility::reinterpolate(statrefx, statx[j], staty[j]);});
                            // Reduction
                            MPI_Allreduce(&statsize, &statgsize, one, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
                            std::fill(statmean.begin(), statmean.end(), real());
                            for (uint j = zero; j < statsize; ++j) {
                                for (uint k = zero; k < statlength; ++k) {
                                    statmean[k] += staty[j][k];
                                }
                            }
                            MPI_Allreduce(statmean.data(), statgmean.data(), statlength, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                            for (uint k = zero; k < statlength; ++k) {
                                statmean[k] /= real(statsize);
                                statgmean[k] /= real(statgsize);
                            }
                            std::fill(statstd.begin(), statstd.end(), real());
                            for (uint j = zero; j < statsize; ++j) {
                                for (uint k = zero; k < statlength; ++k) {
                                    statstd[k] += (staty[j][k]-statgmean[k])*(staty[j][k]-statgmean[k]);
                                }
                            }
                            MPI_Allreduce(statstd.data(), statgstd.data(), statlength, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                            std::fill(statstd.begin(), statstd.end(), real());
                            for (uint j = zero; j < statsize; ++j) {
                                for (uint k = zero; k < statlength; ++k) {
                                    statstd[k] += (staty[j][k]-statmean[k])*(staty[j][k]-statmean[k]);
                                }
                            }
                            for (uint k = zero; k < statlength; ++k) {
                                statstd[k] = std::sqrt(statstd[k]/real(statsize-one));
                                statgstd[k] = std::sqrt(statgstd[k]/real(statgsize-one));
                            }
                            // Output statistics
                            filename = Output::name(outputdir, outputprefix, outputsep, outputstat, outputsep, std::make_pair(outputint, nbundle), outputsep, std::make_pair(outputopening, opening), outputsep, interp, outputsep, stat, outputsep, std::make_pair(outputint, rank));
                            std::replace(filename.begin(), filename.end(), dot, dotc);
                            stream.open(Output::name(filename, outputsuffix));
                            Output::save(stream, statrefz, statrefy, statmean, statstd, digits, statsize);
                            stream.close();
                            filename = Output::name(outputdir, outputprefix, outputsep, outputstat, outputsep, std::make_pair(outputint, nbundle), outputsep, std::make_pair(outputopening, opening), outputsep, interp, outputsep, stat);
                            std::replace(filename.begin(), filename.end(), dot, dotc);
                            if (rank == zero) {
                                stream.open(Output::name(filename, outputsuffix));
                                Output::save(stream, statrefz, statrefy, statgmean, statgstd, digits, statgsize);
                                stream.close();
                            }
                        }
                    }
                }
            }
        }
    } else if (visualization) {
        Utility::parallelize(octree.begin(), octree.end(), [=, &rvismin, &rvismax](element& e){const real x = std::get<0>(e).position<real, position, extent>(zero)-center[zero]; const real y = std::get<0>(e).position<real, position, extent>(one)-center[one]; const real z = std::get<0>(e).position<real, position, extent>(two)-center[two]; const real l = std::get<0>(e).extent<real, position, extent>()*std::sqrt(dimension)/two; const real r = std::sqrt(x*x+y*y+z*z); if ((r+l < rvismin) || (r-l > rvismax)) {std::get<0>(e).invalidate();}});
        octree.update();
        if (octree.size() > vislim) {
            octree.resize(vislim);
        }
        octree.shrink();
        filename = Output::name(outputdir, outputprefix, outputsep, outputvisu, outputsep, outputtree, outputsuffix);
        stream.open(filename);
        Output::save(stream, octree, digits);
        stream.close();
    } else if (homogeneous) {
        photon = Integrator::launch(x0, y0, z0, x, y, z);
        trajectory.resize(stepcnt);
        filename = Output::name(outputdir, outputprefix, outputsep, outputhomo, outputsep, outputrk4, outputsep, outputint, outputsuffix);
        Utility::parallelize(std::begin(trajectory), std::end(trajectory), [=, &photon](evolution& e){e.append(photon);});
        Utility::parallelize(stepcnt, [=, &trajectory, &cosmology, &octree, &lboxmpch0, &mpc, &h, &stepmin, &stepinc, &stepmul](const uint i){Integrator::integrate<-1>(trajectory[i], cosmology, octree, lboxmpch0*mpc/h, static_cast<uint>(stepmin*std::pow(stepmul, i)+i*stepinc));});
        Utility::parallelize(stepcnt, [=, &trajectory, &filename](const uint i){std::ofstream outputstream(Output::name(std::make_pair(filename, i))); Output::save(outputstream, trajectory[i], digits); outputstream.close();});
        if (savetree > zero) {
            filename = Output::name(outputdir, outputprefix, outputsep, outputhomo, outputsep, outputtree, outputsuffix);
            stream.open(filename);
            Output::save(stream, octree, digits);
            stream.close();
        }
    } else if (schwarzschild) {
        trajectory.resize(yshiftcnt);
        filename = Output::name(outputdir, outputprefix, outputsep, outputschw, outputsep, outputngp, outputsep, outputrk4, outputsep, std::make_pair(outputint, levelmin), outputsep, std::make_pair(outputint, levelmax), outputsep, outputint, outputsuffix);
        Utility::parallelize(yshiftcnt, [=, &trajectory, &yshiftmin, &yshiftmul, &yshiftinc](const uint i){real yshift = std::min(yshiftmin*std::pow(yshiftmul, i)+i*yshiftinc, diameter/two); trajectory[i].append(Integrator::launch(x0, y0+yshift, z0, x, y+yshift, z));});
        Utility::parallelize(yshiftcnt, [=, &trajectory, &cosmology, &octree, &lboxmpc0, &mpc](const uint i){Integrator::integrate<0, true>(trajectory[i], cosmology, octree, lboxmpc0*mpc, nsteps);});
        Utility::parallelize(yshiftcnt, [=, &trajectory, &filename](const uint i){std::ofstream outputstream(Output::name(std::make_pair(filename, i))); Output::save(outputstream, trajectory[i], digits); outputstream.close();});
        trajectory.clear();
        trajectory.resize(yshiftcnt);
        filename = Output::name(outputdir, outputprefix, outputsep, outputschw, outputsep, outputcic, outputsep, outputrk4, outputsep, std::make_pair(outputint, levelmin), outputsep, std::make_pair(outputint, levelmax), outputsep, outputint, outputsuffix);
        Utility::parallelize(yshiftcnt, [=, &trajectory, &yshiftmin, &yshiftmul, &yshiftinc](const uint i){real yshift = std::min(yshiftmin*std::pow(yshiftmul, i)+i*yshiftinc, diameter/two); trajectory[i].append(Integrator::launch(x0, y0+yshift, z0, x, y+yshift, z));});
        Utility::parallelize(yshiftcnt, [=, &trajectory, &cosmology, &octree, &lboxmpc0, &mpc](const uint i){Integrator::integrate<1, true>(trajectory[i], cosmology, octree, lboxmpc0*mpc, nsteps);});
        Utility::parallelize(yshiftcnt, [=, &trajectory, &filename](const uint i){std::ofstream outputstream(Output::name(std::make_pair(filename, i))); Output::save(outputstream, trajectory[i], digits); outputstream.close();});
        trajectory.clear();
        trajectory.resize(yshiftcnt);
        filename = Output::name(outputdir, outputprefix, outputsep, outputschw, outputsep, outputref, outputsep, outputrk4, outputsep, std::make_pair(outputint, levelmin), outputsep, std::make_pair(outputint, levelmax), outputsep, outputint, outputsuffix);
        Utility::parallelize(yshiftcnt, [=, &trajectory, &yshiftmin, &yshiftmul, &yshiftinc](const uint i){real yshift = std::min(yshiftmin*std::pow(yshiftmul, i)+i*yshiftinc, diameter/two); trajectory[i].append(Integrator::launch(x0, y0+yshift, z0, x, y+yshift, z));});
        Utility::parallelize(yshiftcnt, [=, &trajectory, &cosmology, &octree, &lboxmpc0, &mpc, &massmsun](const uint i){Integrator::integrate<0, true>(trajectory[i], cosmology, octree, lboxmpc0*mpc, nsteps, massmsun*msun);});
        Utility::parallelize(yshiftcnt, [=, &trajectory, &filename](const uint i){std::ofstream outputstream(Output::name(std::make_pair(filename, i))); Output::save(outputstream, trajectory[i], digits); outputstream.close();});
        if (savetree > zero) {
            filename = Output::name(outputdir, outputprefix, outputsep, outputschw, outputsep, outputtree, outputsuffix);
            stream.open(filename);
            Output::save(stream, octree, digits);
            stream.close();
        }
    } else if (tdependent) {
    } else if (test) {
        photons.resize(one);
        random.resize(one);
        photon = Integrator::launch(center[zero], center[one], center[two], center[zero]+diameter/two, center[one], center[two]);  
        photons[zero] = Integrator::launch(cone[rank].vertex()[zero], cone[rank].vertex()[one], cone[rank].vertex()[two], cone[rank].base()[zero], cone[rank].base()[one], cone[rank].base()[two]);
        random[zero] = distribution(engine)*two*pi;
        nbundle = nbundlemin;
        opening = openingmin;
        filename = Output::name(outputdir, outputprefix, outputsep, outputtest, outputsep, std::make_pair(outputint, nbundle), outputsep, std::make_pair(outputopening, opening), outputsep, interpolation, outputsep, std::make_pair(outputint, rank));
        std::replace(filename.begin(), filename.end(), dot, dotc);
        homotree.assign(ncoarse/two, zero);
        reference = Integrator::propagate<-1>(photon, nbundle, opening, real(), interpolation, cosmology, homotree, lboxmpch*mpc/h, nsteps*(one << (ncoarse-ncoarse/two))*two, real(), Output::name(filename, outputsuffix));
        homotree.clear();
        homotree.shrink();
        reference = Integrator::propagate(photons[zero], nbundle, opening, random[zero], interpolation, cosmology, octree, lboxmpch*mpc/h, nsteps, amin, std::signbit(savemode) ? Output::name() : Output::name(savemode ? Output::name(filename, outputsep, std::make_pair(outputint, zero), outputsep, outputint) : Output::name(filename, outputsep, std::make_pair(outputint, zero), outputsep, std::make_pair(outputint, zero)), outputsuffix), reference);
        if (!reference.empty()) {
            size = zero;
            count.resize(octree.size());
            Utility::parallelize(reference.begin(), reference.end(), [=, &octree, &count](const Photon<real, dimension>& p){++count[std::distance(octree.begin(), octree.locate(p.x(), p.y(), p.z()))];});
            std::for_each(count.begin(), count.end(), [=, &size](std::atomic<uint>& u){u = (u > zero) ? (++size) : (zero);});
            index.resize(size);
            Utility::parallelize(octree.size(), [=, &octree, &count, &index](const uint i){if (count[i] > zero) {index[count[i]-one] = (std::get<0>(octree[i]).level() > ncoarse) ? (std::get<0>(octree[i]).previous(zero, ncoarse)) : (std::get<0>(octree[i]));}});
            for (uint i = 0; i < size; ++i) {
                for (integer ix = -static_cast<integer>(two); ix <= static_cast<integer>(two); ++ix) {
                    for (integer iy = -static_cast<integer>(two); iy <= static_cast<integer>(two); ++iy) {
                        for (integer iz = -static_cast<integer>(two); iz <= static_cast<integer>(two); ++iz) {
                            index.emplace_back(index[i].make(ncoarse, index[i].coordinate(zero)+ix, index[i].coordinate(one)+iy, index[i].coordinate(two)+iz));
                        }
                    }
                }
            }
            std::sort(index.begin(), index.end());
            index.erase(std::unique(index.begin(), index.end()), index.end());
            size = index.size();
            Utility::parallelize(octree.begin(), octree.end(), [=, &index, &ncoarse](element& e){if (!std::binary_search(index.begin(), index.end(), (std::get<0>(e).level() > ncoarse) ? (std::get<0>(e).previous(zero, ncoarse)) : (std::get<0>(e)))) {std::get<0>(e).invalidate();}});
            octree.update();
            filename = Output::name(outputdir, outputprefix, outputsep, outputtest, outputsep, outputtree, outputsuffix);
            stream.open(filename);
            Output::save(stream, octree, digits);
            stream.close();
        }
    }
    
    // Finalization
    MPI_Finalize();*/
    return 0;
}
// -------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
