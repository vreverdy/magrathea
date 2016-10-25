/* ******************************* INTEGRATOR ******************************* */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        RAYTRACING
// TITLE :          Integrator
// DESCRIPTION :    Integration utilities for raytracing
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           integrator.h
/// \brief          Integration utilities for raytracing
/// \author         Vincent Reverdy (vince.rev@gmail.com)
/// \date           2012-2013
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef INTEGRATOR_H_INCLUDED
#define INTEGRATOR_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/



// ------------------------------ PREPROCESSOR ------------------------------ //
// Include C++
#include <iostream>
#include <iomanip>
#include <fstream>
#include <type_traits>
#include <limits>
#include <string>
#include <vector>
#include <tuple>
#include <array>
#include <cmath>
#include <random>
// Include libs
// Include project
#include "../magrathea/simplehyperoctree.h"
#include "../magrathea/simplehyperoctreeindex.h"
#include "../magrathea/constants.h"
#include "../magrathea/hypersphere.h"
#include "../magrathea/evolution.h"
#include "utility.h"
#include "gravity.h"
#include "photon.h"
#include "cone.h"
#include "output.h"
// Misc
// -------------------------------------------------------------------------- //



// --------------------------------- CLASS ---------------------------------- //
// Integration utilities for raytracing
/// \brief          Integration utilities for raytracing.
/// \details        Provides a list of integration routines for geodesics 
///                 integration.
class Integrator final
{  
    // Initialization
    /// \name           Initialization
    //@{
    public:
        template <typename Type = double, class Sphere, class Vector, typename Scalar, class Engine, class Distribution, unsigned int Dimension = Sphere::dimension(), class = typename std::enable_if<Dimension == 3>::type> static Photon<Type, Dimension> launch(const Sphere& sphere, const Cone<Vector, Scalar>& cone, Engine& engine, Distribution& distribution);
        template <typename Type = double, class Sphere, class Container, class Vector, typename Scalar, class Engine, class Distribution, unsigned int Dimension = Sphere::dimension(), class = typename std::enable_if<(Dimension == 3) && (std::is_convertible<typename std::remove_cv<typename std::remove_reference<decltype(std::declval<Container>()[0])>::type>::type, Cone<Vector, Scalar> >::value)>::type> static Photon<Type, Dimension> launch(const Sphere& sphere, const Cone<Vector, Scalar>& cone, const Container& cones, Engine& engine, Distribution& distribution);
        template <typename Type, unsigned int Dimension = 3, class = typename std::enable_if<Dimension == 3>::type> static Photon<Type, Dimension> launch(const Type xbegin, const Type ybegin, const Type zbegin, const Type xend, const Type yend, const Type zend);
        template <bool Center = false, typename Type, unsigned int Dimension, class Container = std::vector<Photon<Type, Dimension> >, class = typename std::enable_if<Dimension == 3>::type> static Container launch(const Photon<Type, Dimension>& photon, const unsigned int count, const Type angle, const Type rotation = Type());
    //@}
    
    // Computation
    /// \name           Computation
    //@{
    public:
        template <int Order = 1, class Array, class Cosmology, class Octree, class Type, class Schwarzschild = std::true_type, unsigned int Dimension = Octree::dimension(), class Data = typename std::tuple_element<1, decltype(Octree::element())>::type, class Position = decltype(Octree::position()), class Extent = decltype(Octree::extent()), class = typename std::enable_if<((std::is_arithmetic<Schwarzschild>::value) || (std::is_same<Schwarzschild, std::true_type>::value)) && (Dimension == 3) && (std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<Array>()[0])>::type>::type>::value)>::type> static Array& dphotondl(Array& output, const Array& input, const Cosmology& cosmology, const Octree& octree, const Type length, const Type dl, const Type phi, const Schwarzschild mass = Schwarzschild());
    //@}

    // Evolution
    /// \name           Evolution
    //@{
    public:
        template <int Order = 1, bool RK4 = true, bool Verbose = false, class Cosmology, class Octree, class Type, class Trajectory, class Schwarzschild = std::true_type, unsigned int Dimension = Octree::dimension(), class Element = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<Trajectory>()[0])>::type>::type, class Data = typename std::tuple_element<1, decltype(Octree::element())>::type, class Core = decltype(Element::template type<1>()), unsigned int Size = std::tuple_size<Core>::value, class Position = decltype(Octree::position()), class Extent = decltype(Octree::extent()), class = typename std::enable_if<((std::is_arithmetic<Schwarzschild>::value) || (std::is_same<Schwarzschild, std::true_type>::value)) && (Dimension == 3)>::type> static Trajectory& integrate(Trajectory& trajectory, const Cosmology& cosmology, const Octree& octree, const Type length, const unsigned int nsteps = 1, const Schwarzschild mass = Schwarzschild());
        template <int Order = 1, bool RK4 = true, bool Verbose = false, class Cosmology, class Octree, class Type, unsigned int Dimension, class Homogeneous = std::vector<Photon<Type, Dimension> >, class = typename std::enable_if<(Dimension == 3) && (Dimension == Octree::dimension())>::type> static magrathea::Evolution<Photon<Type, Dimension> > propagate(const Photon<Type, Dimension>& photon, const unsigned int count, const Type angle, const Type rotation, const std::string& interpolation, const Cosmology& cosmology, const Octree& octree, const Type length, const unsigned int nsteps = 1, const Type amin = Type(), const std::string& filenames = std::string(), const Homogeneous& homogeneous = Homogeneous());
    //@}         
    
    // Test
    /// \name           Test
    //@{
    public:
        static int example();
    //@}
};
// -------------------------------------------------------------------------- //



// ----------------------------- INITIALIZATION ----------------------------- //
// Launch a photon in a cone
/// \brief          Launch a photon in a cone.
/// \details        Launches a photon in the provided cone.
/// \tparam         Type Photon type.
/// \tparam         Dimension Number of space dimension.
/// \tparam         Vector Position vector type.
/// \tparam         Scalar Scalar data type.
/// \tparam         Engine Random engine type.
/// \tparam         Distribution Random distribution type.
/// \param[in]      sphere Sphere to pick up the center and the surface.
/// \param[in]      cone Current cone.
/// \param[in,out]  engine Random engine.
/// \param[in,out]  distribution Random distribution.
/// \return         Initial photon in the cone. 
template <typename Type, class Sphere, class Vector, typename Scalar, class Engine, class Distribution, unsigned int Dimension, class> 
Photon<Type, Dimension> Integrator::launch(const Sphere& sphere, const Cone<Vector, Scalar>& cone, Engine& engine, Distribution& distribution)
{
    // Initialization
    static const Type zero = Type();
    static const Type one = Type(1);
    static const unsigned int x = 0;
    static const unsigned int y = 1;
    static const unsigned int z = 2;
    Vector position = Vector();
    Photon<Type, Dimension> result;
    
    // Compute the photon direction
    do {
        position = sphere.template random<Dimension-1>(engine, distribution);
    } while (!cone.inside(position));
    
    // Set the photon cosmology
    result.a() = one;
    
    // Set the photon position
    result.t() = zero;
    result.x() = sphere.center(x);
    result.y() = sphere.center(y);
    result.z() = sphere.center(z);
    
    // Set the photon direction
    result.dtdl() = one;
    result.dxdl() = (position[x]-sphere.center(x))/sphere.radius();
    result.dydl() = (position[y]-sphere.center(y))/sphere.radius();
    result.dzdl() = (position[z]-sphere.center(z))/sphere.radius();

    // Finalization
    return result;
}

// Launch a photon in a serie of cones
/// \brief          Launch a photon in a serie of cones.
/// \details        Launches a photon in the provided cone with the guarantee
///                 that this cone is the closest one comparatively to a list
///                 of cones.
/// \tparam         Type Photon type.
/// \tparam         Dimension Number of space dimension.
/// \tparam         Sphere Sphere type.
/// \tparam         Container Container of cones type.
/// \tparam         Vector Position vector type.
/// \tparam         Scalar Scalar data type.
/// \tparam         Engine Random engine type.
/// \tparam         Distribution Random distribution type.
/// \param[in]      sphere Sphere to pick up the center and the surface.
/// \param[in]      cone Current cone.
/// \param[in]      cones List of cones to compare distance.
/// \param[in,out]  engine Random engine.
/// \param[in,out]  distribution Random distribution.
/// \return         Initial photon in the cone. 
template <typename Type, class Sphere, class Container, class Vector, typename Scalar, class Engine, class Distribution, unsigned int Dimension, class> 
Photon<Type, Dimension> Integrator::launch(const Sphere& sphere, const Cone<Vector, Scalar>& cone, const Container& cones, Engine& engine, Distribution& distribution)
{
    // Initialization
    static const Type zero = Type();
    static const Type one = Type(1);
    static const unsigned int x = 0;
    static const unsigned int y = 1;
    static const unsigned int z = 2;
    const unsigned int size = cones.size();
    Vector position = Vector();
    Scalar length = Scalar();
    Scalar reference = Scalar();
    Scalar distance = Scalar();  
    bool ok = false;
    Photon<Type, Dimension> result;
    
    // Compute the photon direction
    do {
        position = sphere.template random<Dimension-1>(engine, distribution);
        if (cone.inside(position)) {
            ok = true;
            for (unsigned int idim = 0; idim < Dimension; ++idim) {
                length += (cone.base(idim)-cone.vertex(idim))*(position[idim]-cone.vertex(idim));
            }
            length /= cone.template pow<2>(cone.length());
            for (unsigned int idim = 0; idim < Dimension; ++idim) {
                reference += cone.template pow<2>(position[idim]-(cone.vertex(idim)+(cone.base(idim)-cone.vertex(idim))*length));
            }
            for (unsigned int icone = 0; icone < size; ++icone) {
                if (&cones[icone] != &cone) {
                    for (unsigned int idim = 0; idim < Dimension; ++idim) {
                        length += (cones[icone].base(idim)-cone.vertex(idim))*(position[idim]-cones[icone].vertex(idim));
                    }
                    length /= cones[icone].template pow<2>(cones[icone].length());
                    for (unsigned int idim = 0; idim < Dimension; ++idim) {
                        distance += cones[icone].template pow<2>(position[idim]-(cones[icone].vertex(idim)+(cones[icone].base(idim)-cones[icone].vertex(idim))*length));
                    }
                    if (distance < reference) {
                        ok = false;
                        icone = size;
                    }
                }
            }
        }
    } while (!ok);
    
    // Set the photon cosmology
    result.a() = one;
    
    // Set the photon position
    result.t() = zero;
    result.x() = sphere.center(x);
    result.y() = sphere.center(y);
    result.z() = sphere.center(z);
    
    // Set the photon direction
    result.dtdl() = one;
    result.dxdl() = (position[x]-sphere.center(x))/sphere.radius();
    result.dydl() = (position[y]-sphere.center(y))/sphere.radius();
    result.dzdl() = (position[z]-sphere.center(z))/sphere.radius();

    // Finalization
    return result;
}

// Launch a photon going from a position to another
/// \brief          Launch a photon going from a position to another.
/// \details        Launches a photon starting from a point and going to another
///                 one.
/// \tparam         Type Photon type.
/// \tparam         Dimension Number of space dimension.
/// \param[in]      xbegin Starting x coordinate.
/// \param[in]      ybegin Starting y coordinate.
/// \param[in]      zbegin Starting z coordinate.
/// \param[in]      xend Ending x coordinate.
/// \param[in]      yend Ending y coordinate.
/// \param[in]      zend Ending z coordinate.
/// \return         Initial photon.
template <typename Type, unsigned int Dimension, class>
Photon<Type, Dimension> Integrator::launch(const Type xbegin, const Type ybegin, const Type zbegin, const Type xend, const Type yend, const Type zend)
{
    // Initialization
    static const Type zero = Type();
    static const Type one = Type(1);
    const Type norm = std::sqrt((xend-xbegin)*(xend-xbegin)+(yend-ybegin)*(yend-ybegin)+(zend-zbegin)*(zend-zbegin)); 
    Photon<Type, Dimension> result;

    // Set the photon cosmology
    result.a() = one;
    
    // Set the photon position
    result.t() = zero;
    result.x() = xbegin;
    result.y() = ybegin;
    result.z() = zbegin;
    
    // Set the photon direction
    result.dtdl() = one;
    result.dxdl() = (xend-xbegin)/norm;
    result.dydl() = (yend-ybegin)/norm;
    result.dzdl() = (zend-zbegin)/norm;
    
    // Finalization
    return result;
}

// Launch photons around one photon
/// \brief          Launch photons around one photon.
/// \details        Launches a group of photon on a cone around one photon.
/// \tparam         Center Adds the center photon in the resulting vector if
///                 true.
/// \tparam         Type Photon type.
/// \tparam         Dimension Number of space dimension.
/// \tparam         Container Container of photons type.
/// \param[in]      photon Central photon.
/// \param[in]      count Number of photons to return.
/// \param[in]      angle Half-angle at the cone vertex.
/// \param[in]      rotation Arbitrary rotation to optionally apply on the
///                 resulting circle of photons.
/// \return         Circle of photons.
template <bool Center, typename Type, unsigned int Dimension, class Container, class>
Container Integrator::launch(const Photon<Type, Dimension>& photon, const unsigned int count, const Type angle, const Type rotation)
{
    // Initialization
    static const unsigned int zero = 0;
    static const Type one = 1;
    static const Type four = 4;
    static const Type pi = four*std::atan(one);
    const Type step = (pi+pi)/Type(count+(!count));
    const Type r = std::sqrt((photon.dxdl()*photon.dxdl())+(photon.dydl()*photon.dydl())+(photon.dzdl()*photon.dzdl()));
    const Type rcos = r*std::cos(angle);
    const Type rsin = r*std::sin(angle);
    const Type costheta = std::cos(std::acos(photon.dzdl()/r));
    const Type sintheta = std::sin(std::acos(photon.dzdl()/r));
    const Type cosphi = std::cos(std::atan2(photon.dydl(), photon.dxdl()));
    const Type sinphi = std::sin(std::atan2(photon.dydl(), photon.dxdl()));
    const Type cospsi = std::cos(rotation);
    const Type sinpsi = std::sin(rotation);
    Type x = zero;
    Type y = zero;
    Type z = zero;
    Container result(count+Center, photon);
    
    // Loop over photons
    if (r > zero) {
        for (unsigned int istep = zero; istep < count; ++istep) {
            x = rsin*std::cos(istep*step);
            y = rsin*std::sin(istep*step);
            z = rcos;
            result[istep+Center].dxdl() = -cosphi*sinpsi*costheta*x-cosphi*cospsi*costheta*y+cosphi*sintheta*z+sinphi*sinpsi*y-sinphi*cospsi*x;
            result[istep+Center].dydl() = -sinphi*sinpsi*costheta*x-sinphi*cospsi*costheta*y-cosphi*sinpsi*y+cosphi*cospsi*x+sintheta*sinphi*z;
            result[istep+Center].dzdl() = sintheta*sinpsi*x+sintheta*cospsi*y+costheta*z;
        }
    }
    
    // Finalization
    return result;
}
// -------------------------------------------------------------------------- //



// ------------------------------- COMPUTATION ------------------------------ //
// Derivative of a photon
/// \brief          Derivative of a photon.
/// \details        Computes the derivative of the core components of a photon.
/// \tparam         Order Octree interpolation order :  0 for NGP, 1 for CIC
///                 or -1 for an homogeneous universe.
/// \tparam         Array Core array type.
/// \tparam         Cosmology Cosmology evolution type.
/// \tparam         Octree Octree type.
/// \tparam         Type Scalar type.
/// \tparam         Dimension Number of space dimension.
/// \tparam         Data Data type.
/// \tparam         Schwarzschild Optional schwarzschild type.
/// \tparam         Position Position of the hyperoctree center.
/// \tparam         Extent Extent of the hyperoctree.
/// \param[in,out]  output Output data.
/// \param[in]      input Input data.
/// \param[in]      cosmology Cosmology evolution.
/// \param[in]      octree Octree.
/// \param[in]      length Spatial length.
/// \param[in]      dl Lambda step.
/// \param[in]      phi Value of the previous phi.
/// \param[in]      mass Optional schwarzschild mass put in the center of the 
///                 the hyperoctree.
/// \return         Reference to the output data.
template <int Order, class Array, class Cosmology, class Octree, class Type, class Schwarzschild, unsigned int Dimension, class Data, class Position, class Extent, class> 
Array& Integrator::dphotondl(Array& output, const Array& input, const Cosmology& cosmology, const Octree& octree, const Type length, const Type dl, const Type phi, const Schwarzschild mass)
{
    // Initialization
    static const Type position = Type(Position::num)/Type(Position::den);
    static const Type extent = Type(Extent::num)/Type(Extent::den);
    static const std::array<Type, Dimension> center = std::array<Type, Dimension>({{position, position, position}});
    static const unsigned int a = 0;
    static const unsigned int t = 1;
    static const unsigned int x = 2;
    static const unsigned int y = 3;
    static const unsigned int z = 4;
    static const unsigned int dtdl = 5;
    static const unsigned int dxdl = 6;
    static const unsigned int dydl = 7;
    static const unsigned int dzdl = 8;
    static const Type zero = 0;
    static const Type one = 1;
    static const Type two = 2;
    static const Type c2 = magrathea::Constants<Type>::c2();
    static const Type g = magrathea::Constants<Type>::g();
    Data data = (Order == 0) ? (octree.ngp(input[x], input[y], input[z])) : ((Order == 1) ? (octree.cic(input[x], input[y], input[z])) : (Data()));
    const Type dphidl = ((dl > zero) || (dl < zero)) ? ((data.phi()-phi)/dl) : (phi);
    const Type dadt = Utility::interpolate(input[t], std::get<0>(cosmology), std::get<2>(cosmology));
    const Type scale = length/extent;
    const Type distance = (std::is_arithmetic<Schwarzschild>::value) ? (Utility::distance<Dimension>(center, std::array<Type, Dimension>({{input[x], input[y], input[z]}}))*scale) : (zero);

    // Exact schwarzschild mass computation
    if (std::is_arithmetic<Schwarzschild>::value) {
        data.rho() = zero;
        data.phi() = (distance > zero) ? (-(g*mass)/distance) : (std::numeric_limits<Type>::max());
        data.dphidx() = (distance > zero) ? (((g*mass)/(distance*distance))*((input[x]-center[0])*(scale/distance))) : (zero);
        data.dphidy() = (distance > zero) ? (((g*mass)/(distance*distance))*((input[y]-center[1])*(scale/distance))) : (zero);
        data.dphidz() = (distance > zero) ? (((g*mass)/(distance*distance))*((input[z]-center[2])*(scale/distance))) : (zero);
        data.a() = one;
    }
    
    // Computation
    output[a] = input[dtdl]*dadt;
    output[t] = input[dtdl];
    output[x] = input[dxdl]/scale;
    output[y] = input[dydl]/scale;
    output[z] = input[dzdl]/scale;
    output[dtdl] = -(two*dadt/input[a]*input[dtdl]*input[dtdl])-(two/c2*input[dtdl])*(data.dphidx()*input[dxdl]+data.dphidy()*input[dydl]+data.dphidz()*input[dzdl]);
    output[dxdl] = -(two*dadt/input[a]*input[dtdl]*input[dxdl])+(two/c2*dphidl*input[dxdl])-(two*data.dphidx()*input[dtdl]*input[dtdl]);
    output[dydl] = -(two*dadt/input[a]*input[dtdl]*input[dydl])+(two/c2*dphidl*input[dydl])-(two*data.dphidy()*input[dtdl]*input[dtdl]);
    output[dzdl] = -(two*dadt/input[a]*input[dtdl]*input[dzdl])+(two/c2*dphidl*input[dzdl])-(two*data.dphidz()*input[dtdl]*input[dtdl]);
    
    // Finalization
    return output;
}
// -------------------------------------------------------------------------- //



// -------------------------------- EVOLUTION ------------------------------- //
// Geodesics integration
/// \brief          Geodesics integration.
/// \details        Integrates the geodesics equation of a photon.
/// \tparam         Order Octree interpolation order :  0 for NGP, 1 for CIC
///                 or -1 for an homogeneous universe.
/// \tparam         RK4 Runge-kutta of fourth order or euler.
/// \tparam         Verbose Verbose mode for debug purposes.
/// \tparam         Cosmology Cosmology evolution type.
/// \tparam         Octree Octree type.
/// \tparam         Type Scalar type.
/// \tparam         Trajectory Trajectory type.
/// \tparam         Schwarzschild Optional schwarzschild type.
/// \tparam         Dimension Number of space dimension.
/// \tparam         Element Photon type.
/// \tparam         Data Data type.
/// \tparam         Core Core data type.
/// \tparam         Size Core size.
/// \tparam         Position Position of the hyperoctree center.
/// \tparam         Extent Extent of the hyperoctree.
/// \param[in,out]  trajectory Trajectory.
/// \param[in]      cosmology Cosmology evolution.
/// \param[in]      octree Octree.
/// \param[in]      length Spatial length in SI units.
/// \param[in]      nsteps Number of lambda steps per grid.
/// \param[in]      mass Optional schwarzschild mass put in the center of the 
///                 the hyperoctree.
/// \return         Reference to the trajectory data.
template <int Order, bool RK4, bool Verbose, class Cosmology, class Octree, class Type, class Trajectory, class Schwarzschild, unsigned int Dimension, class Element, class Data, class Core, unsigned int Size, class Position, class Extent, class> 
Trajectory& Integrator::integrate(Trajectory& trajectory, const Cosmology& cosmology, const Octree& octree, const Type length, const unsigned int nsteps, const Schwarzschild mass)
{
    // Initialization
    static const Type zero = 0;
    static const Type one = 1;
    static const Type two = 2;
    static const Type six = 6;
    static const Type c = magrathea::Constants<Type>::c();
    static const Type c2 = magrathea::Constants<Type>::c2();
    static const Type g = magrathea::Constants<Type>::g();
    static const Type position = Type(Position::num)/Type(Position::den);
    static const Type extent = Type(Extent::num)/Type(Extent::den);
    static const Type min = position-(extent/two);
    static const Type max = position+(extent/two);
    static const Data empty = Data();
    static const Data homogeneous = empty.copy().a(one);
    static const std::array<Type, Dimension> center = std::array<Type, Dimension>({{position, position, position}});
    const Type scale = length/extent;
    Type norm = Type();
    Data data = Data();
    Element photon = Element();
    Type ratio = Type();
    Type distance = Type();
    Type dl = Type();
    Type gref = Type(); 
    std::array<Core, 4> dcoredl = std::array<Core, 4>();

    // Integrate 
    if (!trajectory.empty()) {
        
        // Get initial data
        data = (Order == 0) ? (octree.ngp(trajectory.back().x(), trajectory.back().y(), trajectory.back().z())) : ((Order == 1) ? (octree.cic(trajectory.back().x(), trajectory.back().y(), trajectory.back().z())) : (homogeneous));
        if (std::is_arithmetic<Schwarzschild>::value) {
            distance = Utility::distance<Dimension>(center, std::array<Type, Dimension>({{trajectory.back().x(), trajectory.back().y(), trajectory.back().z()}}))*scale;
            data.rho() = zero;
            data.phi() = (distance > zero) ? (-(g*mass)/distance) : (std::numeric_limits<Type>::max());
            data.dphidx() = (distance > zero) ? (((g*mass)/(distance*distance))*((trajectory.back().x()-center[0])*(scale/distance))) : (zero);
            data.dphidy() = (distance > zero) ? (((g*mass)/(distance*distance))*((trajectory.back().y()-center[1])*(scale/distance))) : (zero);
            data.dphidz() = (distance > zero) ? (((g*mass)/(distance*distance))*((trajectory.back().z()-center[2])*(scale/distance))) : (zero);
            data.a() = one;
        }
        norm = std::sqrt((c2*(one+two/c2*data.phi())*trajectory.back().dtdl()*trajectory.back().dtdl())/((one-two/c2*data.phi())*(trajectory.back().dxdl()*trajectory.back().dxdl()+trajectory.back().dydl()*trajectory.back().dydl()+trajectory.back().dzdl()*trajectory.back().dzdl())));
        trajectory.back().dxdl() *= norm;
        trajectory.back().dydl() *= norm;
        trajectory.back().dzdl() *= norm;
        trajectory.back().a() = Utility::interpolate(trajectory.back().t(), std::get<0>(cosmology), std::get<1>(cosmology));
        trajectory.back().level() = std::get<0>(*octree.locate(trajectory.back().x(), trajectory.back().y(), trajectory.back().z())).level();
        trajectory.back().ah() = data.a();
        trajectory.back().rho() = data.rho();
        trajectory.back().phi() = data.phi();
        trajectory.back().dphidx() = data.dphidx();
        trajectory.back().dphidy() = data.dphidy();
        trajectory.back().dphidz() = data.dphidz();
        trajectory.back().dphidl() = zero;
        trajectory.back().laplacian() = zero;
        trajectory.back().redshift() = zero;
        trajectory.back().dsdl2() = (trajectory.back().a()*trajectory.back().a())*(-(c2*(one+two/c2*trajectory.back().phi())*trajectory.back().dtdl()*trajectory.back().dtdl())+((one-two/c2*trajectory.back().phi())*(trajectory.back().dxdl()*trajectory.back().dxdl()+trajectory.back().dydl()*trajectory.back().dydl()+trajectory.back().dzdl()*trajectory.back().dzdl())));
        trajectory.back().error() = one-((one-two/c2*trajectory.back().phi())*(trajectory.back().dxdl()*trajectory.back().dxdl()+trajectory.back().dydl()*trajectory.back().dydl()+trajectory.back().dzdl()*trajectory.back().dzdl()))/(c2*(one+two/c2*trajectory.back().phi())*trajectory.back().dtdl()*trajectory.back().dtdl());
        trajectory.back().distance() = zero;
        trajectory.back().major() = zero;
        trajectory.back().minor() = zero;
        trajectory.back().rotation() = zero;
        ratio = (Order == 0) ? (data.a()*data.a()*(scale/c)/nsteps) : ((Order == 1) ? (data.a()*data.a()*(scale/c)/nsteps) : (trajectory.back().a()*trajectory.back().a()*(scale/c)/nsteps));
        dl = std::get<0>(*octree.locate(trajectory.back().x(), trajectory.back().y(), trajectory.back().z())).template extent<Type, Position, Extent>()*ratio;
        gref = -trajectory.back().a()*c*trajectory.back().dtdl()*(one+trajectory.back().phi()/c2); 

        // Advance
        while (data != empty) {

            // Photon index
            photon.index() = trajectory.back().index()+one;
            
            // Photon core
            if (RK4) {
                dphotondl<Order>(dcoredl[0], trajectory.back().core(), cosmology, octree, length, zero, trajectory.back().dphidl(), mass);
                for (unsigned int i = 0; i < Size; ++i) {
                    photon.core(i) = trajectory.back().core(i)+dl/two*dcoredl[0][i];
                }
                dphotondl<Order>(dcoredl[1], photon.core(), cosmology, octree, length, dl/two, trajectory.back().phi(), mass);
                for (unsigned int i = 0; i < Size; ++i) {
                    photon.core(i) = trajectory.back().core(i)+dl/two*dcoredl[1][i];
                }
                dphotondl<Order>(dcoredl[2], photon.core(), cosmology, octree, length, dl/two, trajectory.back().phi(), mass);
                for (unsigned int i = 0; i < std::tuple_size<Core>::value; ++i) {
                    photon.core(i) = trajectory.back().core(i)+dl*dcoredl[2][i];
                }
                dphotondl<Order>(dcoredl[3], photon.core(), cosmology, octree, length, dl, trajectory.back().phi(), mass);
                for (unsigned int i = 0; i < Size; ++i) {
                    photon.core(i) = trajectory.back().core(i)+(dl/six)*(dcoredl[0][i]+two*dcoredl[1][i]+two*dcoredl[2][i]+dcoredl[3][i]);
                }
            } else {
                dphotondl<Order>(photon.core(), trajectory.back().core(), cosmology, octree, length, zero, trajectory.back().dphidl(), mass);
                for (unsigned int i = 0; i < Size; ++i) {
                    photon.core(i) = trajectory.back().core(i)+dl*photon.core(i);
                }
            }
            
            // Photon extra
            data = (Order == 0) ? (octree.ngp(photon.x(), photon.y(), photon.z())) : ((Order == 1) ? (octree.cic(photon.x(), photon.y(), photon.z())) : (homogeneous));
            if (std::is_arithmetic<Schwarzschild>::value) {
                distance = Utility::distance<Dimension>(center, std::array<Type, Dimension>({{photon.x(), photon.y(), photon.z()}}))*scale;
                data.rho() = zero;
                data.phi() = (distance > zero) ? (-(g*mass)/distance) : (std::numeric_limits<Type>::max());
                data.dphidx() = (distance > zero) ? (((g*mass)/(distance*distance))*((photon.x()-center[0])*(scale/distance))) : (zero);
                data.dphidy() = (distance > zero) ? (((g*mass)/(distance*distance))*((photon.y()-center[1])*(scale/distance))) : (zero);
                data.dphidz() = (distance > zero) ? (((g*mass)/(distance*distance))*((photon.z()-center[2])*(scale/distance))) : (zero);
                data.a() = one;
            }
            data = (!(photon.a() < zero) && ((photon.x() > min) && (photon.x() < max) && (photon.y() > min) && (photon.y() < max) && (photon.z() > min) && (photon.z() < max))) ? (data) : (empty);
            photon.level() = std::get<0>(*octree.locate(photon.x(), photon.y(), photon.z())).level();
            photon.ah() = data.a();
            photon.rho() = data.rho();
            photon.phi() = data.phi();
            photon.dphidx() = data.dphidx();
            photon.dphidy() = data.dphidy();
            photon.dphidz() = data.dphidz();
            photon.dphidl() = (data.phi()-trajectory.back().phi())/dl;
            photon.laplacian() = zero;
            photon.redshift() = -photon.a()*c*photon.dtdl()*(one+photon.phi()/c2)/gref-one;
            photon.dsdl2() = (photon.a()*photon.a())*((-c2*(one+two/c2*photon.phi())*photon.dtdl()*photon.dtdl())+((one-two/c2*photon.phi())*(photon.dxdl()*photon.dxdl()+photon.dydl()*photon.dydl()+photon.dzdl()*photon.dzdl())));
            photon.error() = one-((one-two/c2*photon.phi())*(photon.dxdl()*photon.dxdl()+photon.dydl()*photon.dydl()+photon.dzdl()*photon.dzdl()))/(c2*(one+two/c2*photon.phi())*photon.dtdl()*photon.dtdl());
            photon.distance() = zero;
            photon.major() = zero;
            photon.minor() = zero;
            photon.rotation() = zero;
            
            // Next step
            if (data != empty) {
                ratio = (Order == 0) ? (data.a()*data.a()*(scale/c)/nsteps) : ((Order == 1) ? (data.a()*data.a()*(scale/c)/nsteps) : (photon.a()*photon.a()*(scale/c)/nsteps));
                dl = std::get<0>(*(octree.locate(photon.x(), photon.y(), photon.z()))).template extent<Type, Position, Extent>()*ratio;
                trajectory.append(photon);
                if (Verbose) {
                    if (photon.a() > 0.99 || photon.a() < 0.04) {
                        std::cout<<std::setprecision(17)<<"photon = "<<photon<<" "<<dl<<std::endl;
                        std::cout<<"end = "<<trajectory.back().a()<<" "<<trajectory.back().dtdl()<<" "<<trajectory.back().a()*trajectory.back().dtdl()<<" "<<-photon.a()*c*photon.dtdl()*(one+photon.phi()/c2)/gref-one<<std::endl;
                        std::cout<<"----------------------------------------------------------------------------------------------------------"<<std::endl;
                    }
                }
            }
        }
        
        // Erase last element if non compatible
        if (!trajectory.empty()) {
            if (std::signbit(trajectory.back().redshift()) || (std::signbit(trajectory.back().a()))) {
                trajectory.pop();
            }
        }
    }

    // Finalization
    return trajectory;
}

// Propagation of a ray bundle
/// \brief          Propagation of a ray bundle.
/// \details        Propagates a ray bundle calling the integrator for each
///                 photon.
/// \tparam         Order Octree interpolation order :  0 for NGP, 1 for CIC
///                 or -1 for an homogeneous universe.
/// \tparam         RK4 Runge-kutta of fourth order or euler.
/// \tparam         Verbose Verbose mode for debug purposes.
/// \tparam         Cosmology Cosmology evolution type.
/// \tparam         Octree Octree type.
/// \tparam         Type Scalar type.
/// \tparam         Dimension Number of space dimension.
/// \tparam         Homogeneous Homogeneous reference.
/// \param[in]      photon Central photon initial data.
/// \param[in]      count Number of other photons to use.
/// \param[in]      angle Half-angle at the cone vertex.
/// \param[in]      rotation Arbitrary rotation to optionally apply on the
///                 resulting circle of photons.
/// \param[in]      interpolation Stop condition : redshift, a, t, r.
/// \param[in]      cosmology Cosmology evolution.
/// \param[in]      octree Octree.
/// \param[in]      length Spatial length in SI units.
/// \param[in]      nsteps Number of lambda steps per grid.
/// \param[in]      amin If different from zero, all photons should end by this
///                 value of a.
/// \param[in]      filenames File names of the output. If empty, no output. 
///                 If at least on percent sign, all trajectories are saved.
///                 Otherwise, only the central one is saved.
/// \param[in]      Homogeneous Optional homogeneous trajectory. If not provided
///                 the angular diameter distance use the inhomogeneous value of
///                 a. If provided, the homogeneous value of a for the given 
///                 radius is used.
/// \return         Central photon trajectory.
template <int Order, bool RK4, bool Verbose, class Cosmology, class Octree, class Type, unsigned int Dimension, class Homogeneous, class>
magrathea::Evolution<Photon<Type, Dimension> > Integrator::propagate(const Photon<Type, Dimension>& photon, const unsigned int count, const Type angle, const Type rotation, const std::string& interpolation, const Cosmology& cosmology, const Octree& octree, const Type length, const unsigned int nsteps, const Type amin, const std::string& filenames, const Homogeneous& homogeneous)
{
    // Initialization
    static const Type zero = 0;
    static const Type one = 1;
    static const Type two = 2;
    static const unsigned int first = 0;
    static const unsigned int center = 0;
    static const unsigned int x = 0;
    static const unsigned int y = 1;
    static const unsigned int z = 2;
    static const char percent = '%';
    static const unsigned int digits = std::numeric_limits<Type>::max_digits10;
    static const Type quarter = (Type(octree.extent().num)/Type(two*octree.extent().den))/two;
    static const Type limit = one/(two*two*two);
    std::vector<Photon<Type, Dimension> > initial = launch<true>(photon, count, angle, rotation); 
    std::vector<magrathea::Evolution<Photon<Type, Dimension> > > trajectories(initial.size());
    unsigned int ntrajectories = trajectories.size();
    unsigned int size = zero;
    std::vector<Type> last(ntrajectories);
    std::vector<std::vector<Type> > ref(ntrajectories);
    std::vector<std::array<std::vector<Type>, Dimension> > xyz(ntrajectories);
    std::array<Type, Dimension> coord = std::array<Type, Dimension>();
    std::pair<std::vector<Type>, std::vector<Type> > flrw;
    std::ofstream stream;
    
    // Integration
    for (unsigned int itrajectory = 0; itrajectory < ntrajectories; ++itrajectory) {
        trajectories[itrajectory].append(initial[itrajectory]);
        integrate<Order, RK4, Verbose>(trajectories[itrajectory], cosmology, octree, length, nsteps);
        size = trajectories[itrajectory].size();
        for (unsigned int idim = 0; idim < Dimension; ++idim) {
            xyz[itrajectory][idim].resize(size);
        }
        for (unsigned int istep = 0; istep < size; ++istep) {
            xyz[itrajectory][x][istep] = trajectories[itrajectory][istep].x();
            xyz[itrajectory][y][istep] = trajectories[itrajectory][istep].y();
            xyz[itrajectory][z][istep] = trajectories[itrajectory][istep].z();
        }
        if (!size) {
            ntrajectories = size;
        } else {
            coord[x] = trajectories[itrajectory].back().x()-trajectories[itrajectory].front().x();
            coord[y] = trajectories[itrajectory].back().y()-trajectories[itrajectory].front().y();
            coord[z] = trajectories[itrajectory].back().z()-trajectories[itrajectory].front().z();
            last[itrajectory] = std::sqrt(coord[x]*coord[x]+coord[y]*coord[y]+coord[z]*coord[z]);
            ntrajectories *= (last[itrajectory] > quarter);
            ntrajectories *= (!(std::isnormal(amin) && std::isnormal(trajectories[itrajectory].back().ah()) && (trajectories[itrajectory].back().ah() < one))) || (!(trajectories[itrajectory].back().ah() > amin));
        }
    }
    ntrajectories *= (std::abs((*std::max_element(last.begin(), last.end()))-(*std::min_element(last.begin(), last.end())))/(*std::max_element(last.begin(), last.end())) < limit);
    
    // Interpolation
    if (ntrajectories > 1) {
        if (interpolation == "redshift") {
            for (unsigned int itrajectory = 0; itrajectory < ntrajectories; ++itrajectory) {
                size = trajectories[itrajectory].size();
                ref[itrajectory].resize(size);
                for (unsigned int istep = 0; istep < size; ++istep) {
                    ref[itrajectory][istep] = trajectories[itrajectory][istep].redshift();
                }
            }
        } else if (interpolation == "a") {
            for (unsigned int itrajectory = 0; itrajectory < ntrajectories; ++itrajectory) {
                size = trajectories[itrajectory].size();
                ref[itrajectory].resize(size);
                for (unsigned int istep = 0; istep < size; ++istep) {
                    ref[itrajectory][istep] = trajectories[itrajectory][istep].a();
                }
            }
        } else if ((interpolation == "t") || (interpolation == "eta")) {
            for (unsigned int itrajectory = 0; itrajectory < ntrajectories; ++itrajectory) {
                size = trajectories[itrajectory].size();
                ref[itrajectory].resize(size);
                for (unsigned int istep = 0; istep < size; ++istep) {
                    ref[itrajectory][istep] = trajectories[itrajectory][istep].t();
                }
            }
        } else if ((interpolation == "r") || (interpolation == "radius")) {
            for (unsigned int itrajectory = 0; itrajectory < ntrajectories; ++itrajectory) {
                size = trajectories[itrajectory].size();
                ref[itrajectory].resize(size);
                for (unsigned int istep = 0; istep < size; ++istep) {
                    coord[x] = trajectories[itrajectory][istep].x()-trajectories[itrajectory][first].x();
                    coord[y] = trajectories[itrajectory][istep].y()-trajectories[itrajectory][first].y();
                    coord[z] = trajectories[itrajectory][istep].z()-trajectories[itrajectory][first].z();
                    ref[itrajectory][istep] = std::sqrt(coord[x]*coord[x]+coord[y]*coord[y]+coord[z]*coord[z]);
                }
            }
        } else {
            for (unsigned int itrajectory = 0; itrajectory < ntrajectories; ++itrajectory) {
                size = trajectories[itrajectory].size();
                ref[itrajectory].resize(size);
                for (unsigned int istep = 0; istep < size; ++istep) {
                    ref[itrajectory][istep] = trajectories[itrajectory][istep].redshift();
                }
            }
        }
    }

    // Computation of the angular diameter distance
    if (ntrajectories > 1) {
        size = trajectories[center].size();
        for (unsigned int istep = 0; istep < size; ++istep) {
            for (unsigned int itrajectory = 1; itrajectory < ntrajectories; ++itrajectory) {
                coord[x] = Utility::interpolate(ref[center][istep], ref[itrajectory], xyz[itrajectory][x])-xyz[center][x][istep];
                coord[y] = Utility::interpolate(ref[center][istep], ref[itrajectory], xyz[itrajectory][y])-xyz[center][y][istep];
                coord[z] = Utility::interpolate(ref[center][istep], ref[itrajectory], xyz[itrajectory][z])-xyz[center][z][istep];
                trajectories[center][istep].distance() += std::sqrt(coord[x]*coord[x]+coord[y]*coord[y]+coord[z]*coord[z]);
            }
            trajectories[center][istep].distance() /= (ntrajectories-1);
        }
        if (homogeneous.empty()) {
            for (unsigned int istep = 0; istep < size; ++istep) {
                trajectories[center][istep].distance() *= (length*trajectories[center][istep].a())/angle;
            }
        } else {
            for (unsigned int istep = 0; istep < size; ++istep) {
                coord[x] = xyz[center][x][istep]-xyz[center][x][first];
                coord[y] = xyz[center][y][istep]-xyz[center][y][first];
                coord[z] = xyz[center][z][istep]-xyz[center][z][first];
                ref[0][istep] = std::sqrt(coord[x]*coord[x]+coord[y]*coord[y]+coord[z]*coord[z]);
            }
            size = homogeneous.size();
            flrw.first.resize(size);
            flrw.second.resize(size);
            for (unsigned int istep = 0; istep < size; ++istep) {
                coord[x] = homogeneous[istep].x()-homogeneous[first].x();
                coord[y] = homogeneous[istep].y()-homogeneous[first].y();
                coord[z] = homogeneous[istep].z()-homogeneous[first].z();
                flrw.first[istep] = std::sqrt(coord[x]*coord[x]+coord[y]*coord[y]+coord[z]*coord[z]);
                flrw.second[istep] = homogeneous[istep].a();
            }
            size = trajectories[center].size();
            for (unsigned int istep = 0; istep < size; ++istep) {
                trajectories[center][istep].distance() *= (length*Utility::interpolate(ref[0][istep], flrw.first, flrw.second))/angle;
            }
        }
    }
    
    // Output
    if ((ntrajectories > 0) && (!filenames.empty())) {
        if (filenames.find(percent) == std::string::npos) {
            stream.open(filenames);
            Output::save(stream, trajectories[center], digits); 
            stream.close();
        } else {
            for (unsigned int itrajectory = 0; itrajectory < ntrajectories; ++itrajectory) {
                stream.open(Output::name(std::make_pair(filenames, itrajectory)));
                Output::save(stream, trajectories[itrajectory], digits); 
                stream.close();
            }
        }
    }
    
    // Finalization
    return ntrajectories > 0 ? trajectories[0] : magrathea::Evolution<Photon<Type, Dimension> >();
}
// -------------------------------------------------------------------------- //



// ---------------------------------- TEST ---------------------------------- //
// Example function
/// \brief          Example function.
/// \details        Tests and demonstrates the use of Integrator.
/// \return         0 if no error.
int Integrator::example()
{
    // Initialize
    std::cout<<"BEGIN = Integrator::example()"<<std::endl;
    std::cout<<std::boolalpha<<std::left;
    const unsigned int width = 40;
    std::random_device device;
    std::mt19937 engine(device());
    std::uniform_real_distribution<double> distribution(0, 1);
    std::array<double, 3> beg = std::array<double, 3>({{0., 0., 0.}});
    std::array<double, 3> end = std::array<double, 3>({{16., 23., 42.}});
    std::array<double, 10> array = std::array<double, 10>();
    std::array<std::vector<double>, 4> cosmology = std::array<std::vector<double>, 4>();
    magrathea::SimpleHyperOctree<double, magrathea::SimpleHyperOctreeIndex<unsigned long long int, 3>, Gravity<float, 3> > octree(0, 2);
    magrathea::HyperSphere<3> sphere = magrathea::HyperSphere<3>::unit();
    magrathea::Evolution<Photon<double, 3> > trajectory;
    Photon<double, 3> photon;
    Cone<> cone(beg, end, 0.42);
    std::vector<Cone<> > cones(3, cone);
    double one = 1;
    
    // Construction
    Integrator integrator;

    // Lifecycle and operators
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Lifecycle and operators : "                   <<std::endl;
    std::cout<<std::setw(width)<<"Integrator() : "                              ; Integrator(); std::cout<<std::endl;
    std::cout<<std::setw(width)<<"integrator = Integrator() : "                 ; integrator = Integrator(); std::cout<<std::endl;

    // Initialization
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Initialization : "                                                                  <<std::endl;
    std::cout<<std::setw(width*2)<<"integrator.launch(sphere, cones[0], engine, distribution) : "                       <<integrator.launch(sphere, cones[0], engine, distribution)<<std::endl;
    std::cout<<std::setw(width*2)<<"integrator.launch(sphere, cones[0], cones, engine, distribution) : "                <<integrator.launch(sphere, cones[0], cones, engine, distribution)<<std::endl;
    std::cout<<std::setw(width*2)<<"integrator.launch(beg[0], beg[1], beg[2], end[0], end[1], end[2]) : "               <<integrator.launch(beg[0], beg[1], beg[2], end[0], end[1], end[2])<<std::endl;
    std::cout<<std::setw(width*2)<<"integrator.launch(photon, 3, 0.42, 0.1).size() : "                                  <<integrator.launch(photon, 3, 0.42, 0.1).size()<<std::endl;
    
    // Computation
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Computation : "                                                                     <<std::endl;
    std::cout<<std::setw(width*2)<<"integrator.dphotondl(array, array, cosmology, octree, one, one, one)[0] : "         <<integrator.dphotondl(array, array, cosmology, octree, one, one, one)[0]<<std::endl;

    // Evolution
    std::cout<<std::endl;
    std::cout<<std::setw(width*3)<<"Evolution : "                                                                                                               <<std::endl;
    std::cout<<std::setw(width*3)<<"integrator.integrate(trajectory, cosmology, octree, one, one).size() : "                                                    <<integrator.integrate(trajectory, cosmology, octree, one, one).size()<<std::endl;
    std::cout<<std::setw(width*3)<<"integrator.propagate(photon, 3, 0.42, 0.1, \"a\", cosmology, octree, one, one).size()"                                      <<integrator.propagate(photon, 3, 0.42, 0.1, "a", cosmology, octree, one, one).size()<<std::endl;
    
    // Finalize
    std::cout<<std::noboolalpha<<std::right<<std::endl;
    std::cout<<"END = Integrator::example()"<<std::endl;
    return 0;
}
// -------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
#endif // INTEGRATOR_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/
