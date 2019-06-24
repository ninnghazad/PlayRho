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

#include <PlayRho/Dynamics/Joints/GravityJoint.hpp>
#include <PlayRho/Dynamics/Joints/JointVisitor.hpp>
#include <PlayRho/Dynamics/Body.hpp>
#include <PlayRho/Dynamics/StepConf.hpp>
#include <PlayRho/Dynamics/Contacts/BodyConstraint.hpp>

namespace playrho {
namespace d2 {

// p = attached point, m = mouse point
// C = p - m
// Cdot = v
//	  = v + cross(w, r)
// J = [I r_skew]
// Identity used:
// w k % (rx i + ry j) = w * (-ry i + rx j)


bool GravityJoint::IsOkay(const GravityJointConf& def) noexcept
{
	if (!Joint::IsOkay(def))
	{
		return false;
	}
	return true;
}

GravityJoint::GravityJoint(const GravityJointConf& def):
	Joint{def},
	m_factor{def.factor},
	m_radius{def.radius}
{
}

void GravityJoint::Accept(JointVisitor& visitor) const
{
	visitor.Visit(*this);
}

void GravityJoint::Accept(JointVisitor& visitor)
{
	visitor.Visit(*this);
}

void GravityJoint::InitVelocityConstraints(
	BodyConstraintsMap& bodies,
	const StepConf& step,
	const ConstraintSolverConf&
) {
//	std::cout << __PRETTY_FUNCTION__ << std::endl;

	auto& bodyConstraintA = At(bodies, GetBodyA());
	auto& bodyConstraintB = At(bodies, GetBodyB());

	const auto invMassA = bodyConstraintA->GetInvMass();
	const auto invMassB = bodyConstraintB->GetInvMass();

	const auto posA = bodyConstraintA->GetPosition(); // Origin, not center of mass??
	auto velA = bodyConstraintA->GetVelocity(); // Velocity at origin or center of mass?

	const auto posB = bodyConstraintB->GetPosition();
	auto velB = bodyConstraintB->GetVelocity();

	const auto qA = UnitVec::Get(posA.angular);
	const auto qB = UnitVec::Get(posB.angular);

	m_rA = Rotate(Length2{}, qA);
	m_rB = Rotate(Length2{}, qB);

	// Distance between centers of mass?
	const auto deltaLocation = Length2{(posA.linear + m_rA) - (posB.linear + m_rB)};

	const auto minDistance{0.01};

	const auto uvresult = UnitVec::Get(deltaLocation[0], deltaLocation[1]);
	const auto u = std::get<UnitVec>(uvresult);
	const auto length = std::get<Length>(uvresult);
	const auto distance = length == 0_m ? minDistance:length;

	// Non-Dynamic bodies get a tiny mass of 1 - force is controlled through m_factor by user
	const auto m0 = (invMassA == 0) ? 1 : (1.0/invMassA);
	const auto m1 = (invMassB == 0) ? 1 : (1.0/invMassB);

	// Inspired by real gravity, but without G and with a shift and factor
	//m_impulse = std::max(((1.0 / (distance*distance)) - (1.0 / (m_radius*m_radius))),Real{0}) * m0 * m1 * m_factor * u;

	m_impulse = std::max(((m_factor / (distance*distance)) - (m_factor / (m_radius*m_radius))),Real{0}) * m0 * m1 * u;

	// We have to start this with 0 or SolveVelocityConstraints will not do enough iterations
	//m_lastStep = step.GetTime();
	m_lastStep = 0;

	// std::cout << "GravityJoint: " << (step.doWarmStart?"WARM":"COLD") << " " << m_impulse[0] << "x" << m_impulse[1] << " " << m_factor << " " << m0 << " " << m1
	// << " fd: " <<  std::max(((1.0 / (distance*distance)) - (1.0 / (m_radius*m_radius))),Real{0})
	// << " V: " << velA.linear[0] << "x" << velA.linear[1] << " " << velB.linear[0] << "x" << velB.linear[1]
	// << " P: " << posA.linear[0] << "x" << posA.linear[1] << " " << posB.linear[0] << "x" << posB.linear[1]
	// // << " d: " << length
	// // << " d0: " << (1.0 / (distance*distance))
	// // << " r0: " << (1.0 / (m_radius*m_radius))
	// << std::endl;


	if (step.doWarmStart)
	{
		const auto invRotInertiaA = bodyConstraintA->GetInvRotInertia();
		const auto invRotInertiaB = bodyConstraintB->GetInvRotInertia();

		const auto P = m_impulse * step.GetTime();
		const auto LA = Cross(m_rA, P) * m_inverseRadian;
		const auto LB = Cross(m_rB, P) * m_inverseRadian;
		velA -= Velocity{invMassA * P, invRotInertiaA * LA};
		velB += Velocity{invMassB * P, invRotInertiaB * LB};
	}
	else
	{
		m_impulse = Momentum2{};
	}

	bodyConstraintA->SetVelocity(velA);
	bodyConstraintB->SetVelocity(velB);
}

bool GravityJoint::SolveVelocityConstraints(BodyConstraintsMap& bodies, const StepConf& step)
{
	auto& bodyConstraintA = At(bodies, GetBodyA());
	auto& bodyConstraintB = At(bodies, GetBodyB());

	const auto invMassA = bodyConstraintA->GetInvMass();
	const auto invRotInertiaA = bodyConstraintA->GetInvRotInertia();
	const auto invMassB = bodyConstraintB->GetInvMass();
	const auto invRotInertiaB = bodyConstraintB->GetInvRotInertia();

	auto velA = bodyConstraintA->GetVelocity();
	auto velB = bodyConstraintB->GetVelocity();

	// std::cout << __PRETTY_FUNCTION__
	// << " V: " << velA.linear[0] << "x" << velA.linear[1] << " " << velB.linear[0] << "x" << velB.linear[1]
	// //<< " P: " << posA.linear[0] << "x" << posA.linear[1] << " " << posB.linear[0] << "x" << posB.linear[1]
	// << std::endl;

	const auto ovA = velA;
	const auto ovB = velB;
	const auto P = m_impulse * step.GetTime();
	const auto LA = Cross(m_rA, P) * m_inverseRadian;
	const auto LB = Cross(m_rB, P) * m_inverseRadian;
	velA -= Velocity{invMassA * P, invRotInertiaA * LA};
	velB += Velocity{invMassB * P, invRotInertiaB * LB};

	if(ovB != velB || ovA != velA) {
		bodyConstraintA->SetVelocity(velA);
		bodyConstraintB->SetVelocity(velB);
		return false;
	}

	return true;
}

bool GravityJoint::SolvePositionConstraints(BodyConstraintsMap& bodies, const ConstraintSolverConf& conf) const
{
	NOT_USED(bodies);
	NOT_USED(conf);
	return true;
}


Length2 GravityJoint::GetAnchorA() const
{
	 return GetWorldPoint(*GetBodyA(), GetLocalAnchorA());
}

Length2 GravityJoint::GetAnchorB() const
{
	 return GetWorldPoint(*GetBodyB(), GetLocalAnchorB());
}

Momentum2 GravityJoint::GetLinearReaction() const
{
	return m_impulse;
}

AngularMomentum GravityJoint::GetAngularReaction() const
{
	return AngularMomentum{0};
}

} // namespace d2
} // namespace playrho
