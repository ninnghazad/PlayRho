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

#ifndef PLAYRHO_DYNAMICS_JOINTS_GRAVITYJOINT_HPP
#define PLAYRHO_DYNAMICS_JOINTS_GRAVITYJOINT_HPP

#include <PlayRho/Dynamics/Joints/Joint.hpp>
#include <PlayRho/Dynamics/Joints/GravityJointConf.hpp>

namespace playrho {
namespace d2 {


/// @brief Gravity Joint.
///
/// @details A distance joint constrains two points on two bodies to remain at a
///   fixed distance from each other. You can view this as a massless, rigid rod.
///
/// @ingroup JointsGroup
///
/// @image html distanceJoint.gif
///
class GravityJoint : public Joint
{
public:

	/// @brief Is the given definition okay.
	static bool IsOkay(const GravityJointConf& data) noexcept;

	/// @brief Initializing constructor.
	/// @attention To create or use the joint within a world instance, call that world
	///   instance's create joint method instead of calling this constructor directly.
	/// @sa World::CreateJoint
	GravityJoint(const GravityJointConf& data);

	void Accept(JointVisitor& visitor) const override;
	void Accept(JointVisitor& visitor) override;
	Length2 GetAnchorA() const override;
	Length2 GetAnchorB() const override;
	Momentum2 GetLinearReaction() const override;
	AngularMomentum GetAngularReaction() const override;

	/// @brief Gets the local anchor point relative to body A's origin.
	Length2 GetLocalAnchorA() const noexcept { return {}; }

	/// @brief Gets the local anchor point relative to body B's origin.
	Length2 GetLocalAnchorB() const noexcept { return {}; }

	/// @brief Sets the natural length.
	/// @note Manipulating the length can lead to non-physical behavior when the frequency is zero.
	void SetRadius(Length radius) noexcept;

	/// @brief Gets the length.
	Length GetRadius() const noexcept;

	/// @brief Sets the natural length.
	/// @note Manipulating the length can lead to non-physical behavior when the frequency is zero.
	void SetRotate(Length rotate) noexcept;

	/// @brief Gets the length.
	Length GetRotate() const noexcept;

	/// @brief Sets the maximum force.
	void SetFactor(NonNegative<Real> force) noexcept;

	/// @brief Gets the maximum force.
	NonNegative<Real> GetFactor() const noexcept;

private:

	void InitVelocityConstraints(BodyConstraintsMap& bodies, const playrho::StepConf& step,
								 const ConstraintSolverConf&) override;
	bool SolveVelocityConstraints(BodyConstraintsMap& bodies, const playrho::StepConf& step) override;
	bool SolvePositionConstraints(BodyConstraintsMap& bodies,
								  const ConstraintSolverConf& conf) const override;

	Length m_radius{100}; ///< Radius.
	Real m_factor{10};
	bool m_rotate{false}; ///< Automatically rotate towards force of gravity.

	// Solver shared
	// Solver variables. These are only valid after InitVelocityConstraints called.
	Length2 m_rA; ///< Relative A position.
	Length2 m_rB; ///< Relative B position.

	Momentum2 m_impulse;
	Time m_lastStep;


	constexpr static auto m_inverseRadian{1.0/Radian};
};

inline void GravityJoint::SetRadius(Length radius) noexcept
{
	m_radius = radius;
}

inline Length GravityJoint::GetRadius() const noexcept
{
	return m_radius;
}

inline void GravityJoint::SetRotate(Length rotate) noexcept
{
	m_rotate = rotate;
}

inline Length GravityJoint::GetRotate() const noexcept
{
	return m_rotate;
}

inline void GravityJoint::SetFactor(NonNegative<Real> factor) noexcept
{
	m_factor = factor;
}

inline NonNegative<Real> GravityJoint::GetFactor() const noexcept
{
	return m_factor;
}

} // namespace d2
} // namespace playrho

#endif // PLAYRHO_DYNAMICS_JOINTS_MOUSEJOINT_HPP
