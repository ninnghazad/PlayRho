/*
 * Original work Copyright (c) 2006-2009 Erin Catto http://www.box2d.org
 * Modified work Copyright (c) 2017 Louis Langholtz https://github.com/louis-langholtz/Box2D
 *
 * This software is provided 'as-is', without any express or implied
 * warranty.  In no event will the authors be held liable for any damages
 * arising from the use of this software.
 * Permission is granted to anyone to use this software for any purpose,
 * including commercial applications, and to alter it and redistribute it
 * freely, subject to the following restrictions:
 * 1. The origin of this software must not be misrepresented; you must not
 * claim that you wrote the original software. If you use this software
 * in a product, an acknowledgment in the product documentation would be
 * appreciated but is not required.
 * 2. Altered source versions must be plainly marked as such, and must not be
 * misrepresented as being the original software.
 * 3. This notice may not be removed or altered from any source distribution.
 */

#ifndef RayCastInput_hpp
#define RayCastInput_hpp

/// @file
/// Declaration of the RayCastInput struct.

#include <Box2D/Common/Math.hpp>

namespace box2d
{
    /// @brief Ray-cast input data.
    /// @details The ray extends from p1 to p1 + maxFraction * (p2 - p1).
    struct RayCastInput
    {
        Length2D p1; ///< Point 1.
        Length2D p2; ///< Point 2.
        RealNum maxFraction; ///< Max fraction. Unit interval value (a value between 0 and 1).
    };
    
} // namespace box2d

#endif /* RayCastInput_hpp */