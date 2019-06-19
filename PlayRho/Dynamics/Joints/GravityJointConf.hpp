/*
 * Original work Copyright (c) 2006-2007 Erin Catto http://www.box2d.org
 * Modified work Copyright (c) 2017 Louis Langholtz https://github.com/louis-langholtz/PlayRho
 *
 * This software is provided 'as-is', without any express or implied
 * warranty. In no event will the authors be held liable for any damages
 * arising from the use of this software.
 *
 * Permission is granted to anyone to use this software for any purpose,
 * including commercial applications, and to alter it and redistribute it
 * freely, subject to the following restrictions:
 *
 * 1. The origin of this software must not be misrepresented; you must not
 *	claim that you wrote the original software. If you use this software
 *	in a product, an acknowledgment in the product documentation would be
 *	appreciated but is not required.
 * 2. Altered source versions must be plainly marked as such, and must not be
 *	misrepresented as being the original software.
 * 3. This notice may not be removed or altered from any source distribution.
 */

#ifndef PLAYRHO_DYNAMICS_JOINTS_GRAVITYJOINTCONF_HPP
#define PLAYRHO_DYNAMICS_JOINTS_GRAVITYJOINTCONF_HPP

#include <PlayRho/Dynamics/Joints/JointConf.hpp>
#include <PlayRho/Common/BoundedValue.hpp>
#include <PlayRho/Common/Math.hpp>

namespace playrho {
namespace d2 {

class Body;
class GravityJoint;

/// @brief Gravity joint definition.
/// @details This requires defining an anchor point on both bodies and the non-zero
/// length of the Gravity joint. The definition uses local anchor points so that
/// the initial configuration can violate the constraint slightly. This helps when
//  saving and loading a game.
/// @warning Do not use a zero or short length.
struct GravityJointConf : public JointBuilder<GravityJointConf>
{

	/// @brief Super type.
	using super = JointBuilder<GravityJointConf>;

	PLAYRHO_CONSTEXPR inline GravityJointConf() noexcept: super{JointType::Gravity} {
		collideConnected = true;
	}

	/// @brief Copy constructor.
	GravityJointConf(const GravityJointConf& copy) = default;

	/// @brief Initializing constructor.
	/// @details Initialize the bodies, anchors, and length using the world anchors.
	GravityJointConf(NonNull<Body*> bodyA, NonNull<Body*> bodyB, Length r) noexcept;

	/// @brief Uses the given length.
	PLAYRHO_CONSTEXPR inline GravityJointConf& UseRadius(Length v) noexcept;

	/// @brief Use value for max force.
	PLAYRHO_CONSTEXPR inline GravityJointConf& UseMaxForce(NonNegative<Force> v) noexcept;

	/// @brief Uses the given frequency.
	PLAYRHO_CONSTEXPR inline GravityJointConf& UseFrequency(Positive<Frequency> v) noexcept;

	/// @brief Uses the given damping ratio.
	PLAYRHO_CONSTEXPR inline GravityJointConf& UseDampingRatio(Real v) noexcept;

	/// @brief Natural length between the anchor points.
	Length radius = 1_m;

	/// Max force.
	/// @details
	/// The maximum constraint force that can be exerted
	/// to move the candidate body. Usually you will express
	/// as some multiple of the weight (multiplier * mass * gravity).
	/// @note This may not be negative.
	NonNegative<Force> maxForce = NonNegative<Force>{0_N};

	/// @brief Mass-spring-damper frequency.
	/// @note 0 disables softness.
	Positive<Frequency> frequency = Positive<Frequency>{1_Hz};

	/// @brief Damping ratio.
	/// @note 0 = no damping, 1 = critical damping.
	Real dampingRatio = 0;
};

PLAYRHO_CONSTEXPR inline GravityJointConf& GravityJointConf::UseRadius(Length v) noexcept
{
	radius = v;
	return *this;
}

PLAYRHO_CONSTEXPR inline GravityJointConf& GravityJointConf::UseMaxForce(NonNegative<Force> v) noexcept
{
	maxForce = v;
	return *this;
}

PLAYRHO_CONSTEXPR inline GravityJointConf& GravityJointConf::UseFrequency(Positive<Frequency> v) noexcept
{
	frequency = v;
	return *this;
}

PLAYRHO_CONSTEXPR inline GravityJointConf& GravityJointConf::UseDampingRatio(Real v) noexcept
{
	dampingRatio = v;
	return *this;
}

/// @brief Gets the definition data for the given joint.
/// @relatedalso GravityJoint
GravityJointConf GetGravityJointConf(const GravityJoint& joint) noexcept;

} // namespace d2
} // namespace playrho

#endif // PLAYRHO_DYNAMICS_JOINTS_GravityJOINTCONF_HPP
