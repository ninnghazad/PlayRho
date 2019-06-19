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

	/// @brief Sets the maximum force.
	void SetMaxForce(NonNegative<Force> force) noexcept;

	/// @brief Gets the maximum force.
	NonNegative<Force> GetMaxForce() const noexcept;

	/// @brief Sets frequency.
	void SetFrequency(Positive<Frequency> frequency) noexcept;

	/// @brief Gets the frequency.
	Positive<Frequency> GetFrequency() const noexcept;

	/// @brief Sets the damping ratio.
	void SetDampingRatio(Real ratio) noexcept;

	/// @brief Gets damping ratio.
	Real GetDampingRatio() const noexcept;

private:

	void InitVelocityConstraints(BodyConstraintsMap& bodies, const playrho::StepConf& step,
								 const ConstraintSolverConf&) override;
	bool SolveVelocityConstraints(BodyConstraintsMap& bodies, const playrho::StepConf& step) override;
	bool SolvePositionConstraints(BodyConstraintsMap& bodies,
								  const ConstraintSolverConf& conf) const override;

	/// @brief Gets the effective mass matrix.
	Mass22 GetEffectiveMassMatrix(const BodyConstraint& body) const noexcept;

	Length m_radius; ///< Radius.
	Positive<Frequency> m_frequency = Positive<Frequency>{1_Hz}; ///< Frequency.
	Real m_dampingRatio; ///< Damping ratio.

	NonNegative<Force> m_maxForce = NonNegative<Force>{0_N}; ///< Max force.
	InvMass m_gamma = InvMass{0}; ///< Gamma.


	// Solver shared
	Momentum m_impulse = 0_Ns; ///< Impulse.
	//Momentum2 m_impulse = Momentum2{}; ///< Impulse.

	// Solver variables. These are only valid after InitVelocityConstraints called.
	InvMass m_invGamma; ///< Inverse gamma.
	LinearVelocity m_inverseDistance; ///< Inverse Distance.
	UnitVec m_u; ///< "u" directional.
	Length2 m_rA; ///< Relative A position.
	Length2 m_rB; ///< Relative B position.

	Mass m_mass; ///< Mass.
	//Mass22 m_mass; ///< 2-by-2 mass matrix in kilograms.
	LinearVelocity2 m_C; ///< Velocity constant.
};

inline void GravityJoint::SetRadius(Length radius) noexcept
{
	m_radius = radius;
}

inline Length GravityJoint::GetRadius() const noexcept
{
	return m_radius;
}

inline void GravityJoint::SetMaxForce(NonNegative<Force> force) noexcept
{
	m_maxForce = force;
}

inline NonNegative<Force> GravityJoint::GetMaxForce() const noexcept
{
	return m_maxForce;
}

inline void GravityJoint::SetFrequency(Positive<Frequency> hz) noexcept
{
	m_frequency = hz;
}

inline Positive<Frequency> GravityJoint::GetFrequency() const noexcept
{
	return m_frequency;
}

inline void GravityJoint::SetDampingRatio(Real ratio) noexcept
{
	m_dampingRatio = ratio;
}

inline Real GravityJoint::GetDampingRatio() const noexcept
{
	return m_dampingRatio;
}

/*
/// @brief Gravity Joint.
///
/// @details A gravity joint is used to make a point on a body track a
///   specified world point. This a soft constraint with a maximum
///   force. This allows the constraint to stretch and without
///   applying huge forces.
/// @note This structure is 120-bytes large (using a 4-byte Real on at least one 64-bit
///   architecture/build).
///
/// @ingroup JointsGroup
///
class GravityJoint : public Joint
{
public:

	/// @brief Is the given definition okay.
	static bool IsOkay(const GravityJointConf& def) noexcept;

	/// @brief Initializing constructor.
	/// @attention To create or use the joint within a world instance, call that world
	///   instance's create joint method instead of calling this constructor directly.
	/// @sa World::CreateJoint
	GravityJoint(const GravityJointConf& def);

	void Accept(JointVisitor& visitor) const override;
	void Accept(JointVisitor& visitor) override;
	Length2 GetAnchorA() const override;
	Length2 GetAnchorB() const override;
	Momentum2 GetLinearReaction() const override;
	AngularMomentum GetAngularReaction() const override;
	bool ShiftOrigin(const Length2 newOrigin) override;

	/// @brief Gets the local anchor B.
	Length2 GetLocalAnchorB() const noexcept;

	/// @brief Sets the gravity point.
	void SetGravity(const Length2 gravity) noexcept;

	/// @brief Gets the gravity point.
	Length2 GetTarget() const noexcept;

	/// @brief Sets the maximum force.
	void SetMaxForce(NonNegative<Force> force) noexcept;

	/// @brief Gets the maximum force.
	NonNegative<Force> GetMaxForce() const noexcept;

	/// @brief Sets the frequency.
	void SetFrequency(NonNegative<Frequency> hz) noexcept;

	/// @brief Gets the frequency.
	NonNegative<Frequency> GetFrequency() const noexcept;

	/// @brief Sets the damping ratio.
	void SetDampingRatio(NonNegative<Real> ratio) noexcept;

	/// @brief Gets the damping ratio.
	NonNegative<Real> GetDampingRatio() const noexcept;

private:
	void InitVelocityConstraints(BodyConstraintsMap& bodies, const StepConf& step,
								 const ConstraintSolverConf& conf) override;
	bool SolveVelocityConstraints(BodyConstraintsMap& bodies, const StepConf& step) override;
	bool SolvePositionConstraints(BodyConstraintsMap& bodies,
								  const ConstraintSolverConf& conf) const override;

	/// @brief Gets the effective mass matrix.
	Mass22 GetEffectiveMassMatrix(const BodyConstraint& body) const noexcept;

	Length2 m_targetA; ///< Gravity location (A).
	Length2 m_localAnchorB; ///< Local anchor B.
	NonNegative<Frequency> m_frequency = NonNegative<Frequency>{0_Hz}; ///< Frequency.
	NonNegative<Real> m_dampingRatio = NonNegative<Real>{0}; ///< Damping ratio.
	NonNegative<Force> m_maxForce = NonNegative<Force>{0_N}; ///< Max force.
	InvMass m_gamma = InvMass{0}; ///< Gamma.

	Momentum2 m_impulse = Momentum2{}; ///< Impulse.

	// Solver variables. These are only valid after InitVelocityConstraints called.
	Length2 m_rB; ///< Relative B.
	Mass22 m_mass; ///< 2-by-2 mass matrix in kilograms.
	LinearVelocity2 m_C; ///< Velocity constant.
};

inline Length2 GravityJoint::GetLocalAnchorB() const noexcept
{
	return m_localAnchorB;
}

inline Length2 GravityJoint::GetTarget() const noexcept
{
	return m_targetA;
}

inline void GravityJoint::SetMaxForce(NonNegative<Force> force) noexcept
{
	m_maxForce = force;
}

inline NonNegative<Force> GravityJoint::GetMaxForce() const noexcept
{
	return m_maxForce;
}

inline void GravityJoint::SetFrequency(NonNegative<Frequency> hz) noexcept
{
	m_frequency = hz;
}

inline NonNegative<Frequency> GravityJoint::GetFrequency() const noexcept
{
	return m_frequency;
}

inline void GravityJoint::SetDampingRatio(NonNegative<Real> ratio) noexcept
{
	m_dampingRatio = ratio;
}

inline NonNegative<Real> GravityJoint::GetDampingRatio() const noexcept
{
	return m_dampingRatio;
}
*/
} // namespace d2
} // namespace playrho

#endif // PLAYRHO_DYNAMICS_JOINTS_MOUSEJOINT_HPP
