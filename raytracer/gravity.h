/* ********************************* GRAVITY ******************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        RAYTRACING
// TITLE :          Gravity
// DESCRIPTION :    Gravity cell implementation for raytracing
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           gravity.h
/// \brief          Gravity cell implementation for raytracing
/// \author         Vincent Reverdy (vince.rev@gmail.com)
/// \date           2012-2013
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef GRAVITY_H_INCLUDED
#define GRAVITY_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/



// ------------------------------ PREPROCESSOR ------------------------------ //
// Include C++
#include <iostream>
#include <iomanip>
#include <type_traits>
#include <array>
#include <utility>
// Include libs
// Include project
#include "../magrathea/abstractcontents.h"
#include "../magrathea/euleriancategory.h"
// Misc
// -------------------------------------------------------------------------- //



// ---------------------------------- CLASS --------------------------------- //
// Gravity cell implementation for raytracing
/// \brief          Gravity cell implementation for raytracing.
/// \details        A gravity cell containing the local density and the local 
///                 potential.
/// \tparam         Type Data type.
/// \tparam         Dimension Number of space dimension.
template <typename Type = double, unsigned int Dimension = 3>
class Gravity final
: public magrathea::AbstractContents<Gravity<Type, Dimension>, magrathea::EulerianCategory, Type, Type, std::array<Type, Dimension>, Type>
{
    // Setup 
    public: using magrathea::AbstractContents<Gravity<Type, Dimension>, magrathea::EulerianCategory, Type, Type, std::array<Type, Dimension>, Type>::operator=; 
    
    // Lifecycle 
    /// \name           Lifecycle 
    //@{
    public: 
        template <class... Misc> explicit inline Gravity(Misc&&... misc); 
    //@}
    
    // Data 
    /// \name           Data 
    //@{
    public: 
        template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<magrathea::AbstractContents<Gravity<Type, Dimension>, magrathea::EulerianCategory, Type, Type, std::array<Type, Dimension>, Type> >().template data<0, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type> inline Template rho(Misc&&... misc); 
        template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const magrathea::AbstractContents<Gravity<Type, Dimension>, magrathea::EulerianCategory, Type, Type, std::array<Type, Dimension>, Type> >().template data<0, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type> inline Template rho(Misc&&... misc) const; 
        template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<magrathea::AbstractContents<Gravity<Type, Dimension>, magrathea::EulerianCategory, Type, Type, std::array<Type, Dimension>, Type> >().template data<1, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type> inline Template phi(Misc&&... misc); 
        template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const magrathea::AbstractContents<Gravity<Type, Dimension>, magrathea::EulerianCategory, Type, Type, std::array<Type, Dimension>, Type> >().template data<1, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type> inline Template phi(Misc&&... misc) const; 
        template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<magrathea::AbstractContents<Gravity<Type, Dimension>, magrathea::EulerianCategory, Type, Type, std::array<Type, Dimension>, Type> >().template data<2, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type> inline Template dphidxyz(Misc&&... misc); 
        template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const magrathea::AbstractContents<Gravity<Type, Dimension>, magrathea::EulerianCategory, Type, Type, std::array<Type, Dimension>, Type> >().template data<2, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type> inline Template dphidxyz(Misc&&... misc) const; 
        template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<magrathea::AbstractContents<Gravity<Type, Dimension>, magrathea::EulerianCategory, Type, Type, std::array<Type, Dimension>, Type> >().template data<2, Values...>(0, std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type> inline Template dphidx(Misc&&... misc); 
        template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const magrathea::AbstractContents<Gravity<Type, Dimension>, magrathea::EulerianCategory, Type, Type, std::array<Type, Dimension>, Type> >().template data<2, Values...>(0, std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type> inline Template dphidx(Misc&&... misc) const; 
        template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<magrathea::AbstractContents<Gravity<Type, Dimension>, magrathea::EulerianCategory, Type, Type, std::array<Type, Dimension>, Type> >().template data<2, Values...>(1, std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type> inline Template dphidy(Misc&&... misc); 
        template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const magrathea::AbstractContents<Gravity<Type, Dimension>, magrathea::EulerianCategory, Type, Type, std::array<Type, Dimension>, Type> >().template data<2, Values...>(1, std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type> inline Template dphidy(Misc&&... misc) const; 
        template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<magrathea::AbstractContents<Gravity<Type, Dimension>, magrathea::EulerianCategory, Type, Type, std::array<Type, Dimension>, Type> >().template data<2, Values...>(2, std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type> inline Template dphidz(Misc&&... misc); 
        template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const magrathea::AbstractContents<Gravity<Type, Dimension>, magrathea::EulerianCategory, Type, Type, std::array<Type, Dimension>, Type> >().template data<2, Values...>(2, std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type> inline Template dphidz(Misc&&... misc) const; 
        template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<magrathea::AbstractContents<Gravity<Type, Dimension>, magrathea::EulerianCategory, Type, Type, std::array<Type, Dimension>, Type> >().template data<3, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type> inline Template a(Misc&&... misc); 
        template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const magrathea::AbstractContents<Gravity<Type, Dimension>, magrathea::EulerianCategory, Type, Type, std::array<Type, Dimension>, Type> >().template data<3, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type> inline Template a(Misc&&... misc) const; 
    //@}
    
    // Test
    /// \name           Test
    //@{
    public:
        static int example();
    //@}
};
// -------------------------------------------------------------------------- //



// ------------------------------- LIFECYCLE -------------------------------- //
// Explicit generic constructor 
/// \brief          Explicit generic constructor. 
/// \details        Provides a generic interface to all constructors of the 
///                 base class. 
/// \tparam         Misc (Miscellaneous types.)
/// \param[in]      misc Miscellaneous arguments.
template <typename Type, unsigned int Dimension>
template <class... Misc>
inline Gravity<Type, Dimension>::Gravity(Misc&&... misc)
: magrathea::AbstractContents<Gravity<Type, Dimension>, magrathea::EulerianCategory, Type, Type, std::array<Type, Dimension>, Type>(std::forward<Misc>(misc)...)
{
    ;
}
// -------------------------------------------------------------------------- //



// ---------------------------------- DATA ---------------------------------- //
// Access to the rho data
/// \brief          Access to the rho data. 
/// \details        Provides an access to the rho data by forwarding 
///                 parameters to the unified base accessor member. 
/// \tparam         Values List of template values. 
/// \tparam         Misc (Miscellaneous types.) 
/// \tparam         Template (Deduced template type.) 
/// \param[in]      misc Miscellaneous arguments. 
/// \return         Forwarded result. 
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Gravity<Type, Dimension>::rho(Misc&&... misc)
{
    return this->template data<0, Values...>(std::forward<Misc>(misc)...);
}

// Immutable access to the rho data
/// \brief          Immutable access to the rho data. 
/// \details        Provides an immutable access to the rho data by forwarding 
///                 parameters to the unified base accessor member. 
/// \tparam         Values List of template values. 
/// \tparam         Misc (Miscellaneous types.) 
/// \tparam         Template (Deduced template type.) 
/// \param[in]      misc Miscellaneous arguments. 
/// \return         Forwarded result. 
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Gravity<Type, Dimension>::rho(Misc&&... misc) const
{
    return this->template data<0, Values...>(std::forward<Misc>(misc)...);
}

// Access to the phi data
/// \brief          Access to the phi data. 
/// \details        Provides an access to the phi data by forwarding 
///                 parameters to the unified base accessor member. 
/// \tparam         Values List of template values. 
/// \tparam         Misc (Miscellaneous types.) 
/// \tparam         Template (Deduced template type.) 
/// \param[in]      misc Miscellaneous arguments. 
/// \return         Forwarded result. 
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Gravity<Type, Dimension>::phi(Misc&&... misc)
{
    return this->template data<1, Values...>(std::forward<Misc>(misc)...);
}

// Immutable access to the phi data
/// \brief          Immutable access to the phi data. 
/// \details        Provides an immutable access to the phi data by forwarding 
///                 parameters to the unified base accessor member. 
/// \tparam         Values List of template values. 
/// \tparam         Misc (Miscellaneous types.) 
/// \tparam         Template (Deduced template type.) 
/// \param[in]      misc Miscellaneous arguments. 
/// \return         Forwarded result. 
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Gravity<Type, Dimension>::phi(Misc&&... misc) const
{
    return this->template data<1, Values...>(std::forward<Misc>(misc)...);
}

// Access to the dphidxyz data
/// \brief          Access to the dphidxyz data. 
/// \details        Provides an access to the dphidxyz data by forwarding 
///                 parameters to the unified base accessor member. 
/// \tparam         Values List of template values. 
/// \tparam         Misc (Miscellaneous types.) 
/// \tparam         Template (Deduced template type.) 
/// \param[in]      misc Miscellaneous arguments. 
/// \return         Forwarded result. 
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Gravity<Type, Dimension>::dphidxyz(Misc&&... misc)
{
    return this->template data<2, Values...>(std::forward<Misc>(misc)...);
}

// Immutable access to the dphidxyz data
/// \brief          Immutable access to the dphidxyz data. 
/// \details        Provides an immutable access to the dphidxyz data by 
///                 forwarding parameters to the unified base accessor member. 
/// \tparam         Values List of template values. 
/// \tparam         Misc (Miscellaneous types.) 
/// \tparam         Template (Deduced template type.) 
/// \param[in]      misc Miscellaneous arguments. 
/// \return         Forwarded result. 
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Gravity<Type, Dimension>::dphidxyz(Misc&&... misc) const
{
    return this->template data<2, Values...>(std::forward<Misc>(misc)...);
}

// Access to the dphidx data
/// \brief          Access to the dphidx data. 
/// \details        Provides an access to the dphidx data by forwarding
///                 parameters to the unified base accessor member. 
/// \tparam         Values List of template values. 
/// \tparam         Misc (Miscellaneous types.) 
/// \tparam         Template (Deduced template type.) 
/// \param[in]      misc Miscellaneous arguments. 
/// \return         Forwarded result. 
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Gravity<Type, Dimension>::dphidx(Misc&&... misc)
{
    return this->template data<2, Values...>(0, std::forward<Misc>(misc)...);
}

// Immutable access to the dphidx data
/// \brief          Immutable access to the dphidx data. 
/// \details        Provides an immutable access to the dphidx data by
///                 forwarding parameters to the unified base accessor member. 
/// \tparam         Values List of template values. 
/// \tparam         Misc (Miscellaneous types.) 
/// \tparam         Template (Deduced template type.) 
/// \param[in]      misc Miscellaneous arguments. 
/// \return         Forwarded result. 
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Gravity<Type, Dimension>::dphidx(Misc&&... misc) const
{
    return this->template data<2, Values...>(0, std::forward<Misc>(misc)...);
}

// Access to the dphidy data
/// \brief          Access to the dphidy data. 
/// \details        Provides an access to the dphidy data by forwarding
///                 parameters to the unified base accessor member. 
/// \tparam         Values List of template values. 
/// \tparam         Misc (Miscellaneous types.) 
/// \tparam         Template (Deduced template type.) 
/// \param[in]      misc Miscellaneous arguments. 
/// \return         Forwarded result. 
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Gravity<Type, Dimension>::dphidy(Misc&&... misc)
{
    return this->template data<2, Values...>(1, std::forward<Misc>(misc)...);
}

// Immutable access to the dphidy data
/// \brief          Immutable access to the dphidy data. 
/// \details        Provides an immutable access to the dphidy data by
///                 forwarding parameters to the unified base accessor member. 
/// \tparam         Values List of template values. 
/// \tparam         Misc (Miscellaneous types.) 
/// \tparam         Template (Deduced template type.) 
/// \param[in]      misc Miscellaneous arguments. 
/// \return         Forwarded result. 
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Gravity<Type, Dimension>::dphidy(Misc&&... misc) const
{
    return this->template data<2, Values...>(1, std::forward<Misc>(misc)...);
}

// Access to the dphidz data
/// \brief          Access to the dphidz data. 
/// \details        Provides an access to the dphidz data by forwarding
///                 parameters to the unified base accessor member. 
/// \tparam         Values List of template values. 
/// \tparam         Misc (Miscellaneous types.) 
/// \tparam         Template (Deduced template type.) 
/// \param[in]      misc Miscellaneous arguments. 
/// \return         Forwarded result. 
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Gravity<Type, Dimension>::dphidz(Misc&&... misc)
{
    return this->template data<2, Values...>(2, std::forward<Misc>(misc)...);
}

// Immutable access to the dphidz data
/// \brief          Immutable access to the dphidz data. 
/// \details        Provides an immutable access to the dphidz data by 
///                 forwarding parameters to the unified base accessor member. 
/// \tparam         Values List of template values. 
/// \tparam         Misc (Miscellaneous types.) 
/// \tparam         Template (Deduced template type.) 
/// \param[in]      misc Miscellaneous arguments. 
/// \return         Forwarded result. 
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Gravity<Type, Dimension>::dphidz(Misc&&... misc) const
{
    return this->template data<2, Values...>(2, std::forward<Misc>(misc)...);
}

// Access to the a data
/// \brief          Access to the a data. 
/// \details        Provides an access to the a data by forwarding parameters 
///                 to the unified base accessor member. 
/// \tparam         Values List of template values. 
/// \tparam         Misc (Miscellaneous types.) 
/// \tparam         Template (Deduced template type.) 
/// \param[in]      misc Miscellaneous arguments. 
/// \return         Forwarded result. 
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Gravity<Type, Dimension>::a(Misc&&... misc)
{
    return this->template data<3, Values...>(std::forward<Misc>(misc)...);
}

// Immutable access to the a data
/// \brief          Immutable access to the a data. 
/// \details        Provides an immutable access to the a data by forwarding 
///                 parameters to the unified base accessor member. 
/// \tparam         Values List of template values. 
/// \tparam         Misc (Miscellaneous types.) 
/// \tparam         Template (Deduced template type.) 
/// \param[in]      misc Miscellaneous arguments. 
/// \return         Forwarded result. 
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Gravity<Type, Dimension>::a(Misc&&... misc) const
{
    return this->template data<3, Values...>(std::forward<Misc>(misc)...);
}
// -------------------------------------------------------------------------- //



// ---------------------------------- TEST ---------------------------------- //
// Example function
/// \brief          Example function.
/// \details        Tests and demonstrates the use of Gravity.
/// \return         0 if no error.
template <typename Type, unsigned int Dimension>
int Gravity<Type, Dimension>::example()
{
    // Initialize
    std::cout<<"BEGIN = Gravity::example()"<<std::endl;
    std::cout<<std::boolalpha<<std::left;
    const unsigned int width = 40;
    std::array<double, 3> array({{4, 8, 15}});

    // Construction
    Gravity<double, 3> gravity(16, 23, array, 42);
    
    // Lifecycle
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Lifecycle : "                                                                       <<std::endl;
    std::cout<<std::setw(width*2)<<"Gravity<>() : "                                                                     <<Gravity<>()<<std::endl;
    std::cout<<std::setw(width*2)<<"Gravity<double, 3>(16, 23, array, 42) : "                                           <<Gravity<double, 3>(16, 23, array, 42)<<std::endl;
    std::cout<<std::setw(width*2)<<"Gravity<double, 3>(gravity) : "                                                     <<Gravity<double, 3>(gravity)<<std::endl;

    // Data
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Data : "                                      <<std::endl;
    std::cout<<std::setw(width)<<"gravity.rho() = 4 : "                         <<(gravity.rho() = 4)<<std::endl;
    std::cout<<std::setw(width)<<"gravity.rho() : "                             <<gravity.rho()<<std::endl;
    std::cout<<std::setw(width)<<"gravity.phi() = 8 : "                         <<(gravity.phi() = 8)<<std::endl;
    std::cout<<std::setw(width)<<"gravity.phi() : "                             <<gravity.phi()<<std::endl;
    std::cout<<std::setw(width)<<"gravity.dphidxyz()[0] = 15 : "                <<(gravity.dphidxyz()[0] = 15)<<std::endl;
    std::cout<<std::setw(width)<<"gravity.dphidxyz()[0] : "                     <<gravity.dphidxyz()[0]<<std::endl;
    std::cout<<std::setw(width)<<"gravity.dphidx() = 16 : "                     <<(gravity.dphidx() = 16)<<std::endl;
    std::cout<<std::setw(width)<<"gravity.dphidx() : "                          <<gravity.dphidx()<<std::endl;
    std::cout<<std::setw(width)<<"gravity.dphidy() = 23 : "                     <<(gravity.dphidy() = 23)<<std::endl;
    std::cout<<std::setw(width)<<"gravity.dphidy() : "                          <<gravity.dphidy()<<std::endl;
    std::cout<<std::setw(width)<<"gravity.dphidz() = 42 : "                     <<(gravity.dphidz() = 42)<<std::endl;
    std::cout<<std::setw(width)<<"gravity.dphidz() : "                          <<gravity.dphidz()<<std::endl;
    std::cout<<std::setw(width)<<"gravity.a() = 4 : "                           <<(gravity.a() = 4)<<std::endl;
    std::cout<<std::setw(width)<<"gravity.a() : "                               <<gravity.a()<<std::endl;
    
    // Finalize
    std::cout<<std::noboolalpha<<std::right<<std::endl;
    std::cout<<"END = Gravity::example()"<<std::endl;
    return 0;
}
// -------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
#endif // GRAVITY_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/

