// Copyright (C) 2002 Charless C. Fowlkes <fowlkes@eecs.berkeley.edu>
// Copyright (C) 2002 David R. Martin <dmartin@eecs.berkeley.edu>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
// 02111-1307, USA, or see http://www.gnu.org/copyleft/gpl.html.

#ifndef TYPES_HH
#define TYPES_HH

#include <limits.h>
#include <float.h>

namespace Util
{
  typedef signed char int8;
  typedef signed short int16;
  typedef signed int int32;
  typedef signed long long int int64;

  typedef unsigned char uint8;
  typedef unsigned short uint16;
  typedef unsigned int uint32;
  typedef unsigned long long int uint64;

  typedef unsigned int uint;

#undef INT8_MIN
#undef INT16_MIN
#undef INT32_MIN
#undef INT64_MIN

#undef INT8_MAX
#undef INT16_MAX
#undef INT32_MAX
#undef INT64_MAX

#undef UINT8_MAX
#undef UINT16_MAX
#undef UINT32_MAX
#undef UINT64_MAX

#define INT8_MIN        ((int8)  0x80)
#define INT16_MIN       ((int16) 0x8000)
#define INT32_MIN       ((int32) 0x80000000)
#define INT64_MIN       ((int64) 0x8000000000000000ll)

#define INT8_MAX        ((int8)  0x7f)
#define INT16_MAX       ((int16) 0x7fff)
#define INT32_MAX       ((int32) 0x7fffffff)
#define INT64_MAX       ((int64) 0x7fffffffffffffffll)

#define UINT8_MAX       ((uint8)  0xff)
#define UINT16_MAX      ((uint16) 0xffff)
#define UINT32_MAX      ((uint32) 0xffffffff)
#define UINT64_MAX      ((uint64) 0xffffffffffffffffull)

}

#endif                          // __types_h__
