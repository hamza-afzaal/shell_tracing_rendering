///////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 1998-2011, Industrial Light & Magic, a division of Lucas
// Digital Ltd. LLC
// 
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// *       Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
// *       Redistributions in binary form must reproduce the above
// copyright notice, this list of conditions and the following disclaimer
// in the documentation and/or other materials provided with the
// distribution.
// *       Neither the name of Industrial Light & Magic nor the names of
// its contributors may be used to endorse or promote products derived
// from this software without specific prior written permission. 
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
///////////////////////////////////////////////////////////////////////////


#include "PyImathVec4ArrayImpl.h"
#include "PyImathExport.h"

namespace PyImath {
using namespace boost::python;
using namespace IMATH_NAMESPACE;

template PYIMATH_EXPORT class_<FixedArray<IMATH_NAMESPACE::Vec4<unsigned char> > > register_Vec4Array<unsigned char>();
template PYIMATH_EXPORT class_<FixedArray<IMATH_NAMESPACE::Vec4<short> > > register_Vec4Array<short>();
template PYIMATH_EXPORT class_<FixedArray<IMATH_NAMESPACE::Vec4<int> > > register_Vec4Array<int>();

template<> PYIMATH_EXPORT IMATH_NAMESPACE::Vec4<unsigned char> FixedArrayDefaultValue<IMATH_NAMESPACE::Vec4<unsigned char> >::value() { return IMATH_NAMESPACE::Vec4<unsigned char>(0,0,0,0); }
template<> PYIMATH_EXPORT IMATH_NAMESPACE::Vec4<short> FixedArrayDefaultValue<IMATH_NAMESPACE::Vec4<short> >::value() { return IMATH_NAMESPACE::Vec4<short>(0,0,0,0); }
template<> PYIMATH_EXPORT IMATH_NAMESPACE::Vec4<int> FixedArrayDefaultValue<IMATH_NAMESPACE::Vec4<int> >::value() { return IMATH_NAMESPACE::Vec4<int>(0,0,0,0); }
}
