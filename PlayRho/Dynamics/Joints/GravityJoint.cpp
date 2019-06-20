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
	m_maxForce{def.maxForce},
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
	const auto deltaLocation = Length2{(posB.linear + m_rB) - (posA.linear + m_rA)};


	m_factor = 10;
	if(invMassA == 0) {
		// const auto mass = ComputeMassData(*GetBodyA());
		// m_massA = mass.mass;
		m_massA = 1;
	//	std::cout << "MASS override A: " << m_massA << std::endl;
	}
	if(invMassB == 0) {
		// const auto mass = ComputeMassData(*GetBodyB());
		// m_massB = mass.mass;
		m_massB = 1;
	//	std::cout << "MASS override B: " << m_massB << std::endl;
	}

	const auto uvresult = UnitVec::Get(deltaLocation[0], deltaLocation[1]);
	m_u = std::get<UnitVec>(uvresult);
	const auto length = std::get<Length>(uvresult);
	const auto distance = length == 0_m ? 0.01:length;

	m_inverseDistance = (1.0 / (distance*distance)) - (1.0 / (m_radius*m_radius));

	//m_inverseDistance = playrho::pow(std::clamp(Real{1} - distance/m_radius,Real{0},Real{1}),2);

	if (step.doWarmStart)
	{
		const auto invRotInertiaA = bodyConstraintA->GetInvRotInertia();
		const auto invRotInertiaB = bodyConstraintB->GetInvRotInertia();

		const auto m0 = (invMassA == 0) ? m_massA : (1.0/invMassA);
		const auto m1 = (invMassB == 0) ? m_massB : (1.0/invMassB);

		m_impulse = LinearAcceleration{
			//(((BigG * m_inverseDistance) * m0) * m1) * m_factor
			(((m_inverseDistance) * m0) * m1) * m_factor
			//m_inverseDistance * m_factor
		} * step.GetTime();

		const auto P = m_impulse * m_u;
		const auto LA = Cross(m_rA, -P) / Radian;
		const auto LB = Cross(m_rB, -P) / Radian;
		velA -= Velocity{invMassA * -P, invRotInertiaA * LA};
		velB += Velocity{invMassB * -P, invRotInertiaB * LB};

	}
	else
	{
		m_impulse = 0;
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

	const auto oldImpulse = m_impulse;

	const auto m0 = (invMassA == 0) ? m_massA : (1.0/invMassA);
	const auto m1 = (invMassB == 0) ? m_massB : (1.0/invMassB);

	m_impulse = LinearAcceleration{
		// (((BigG * m_inverseDistance) * m0) * m1) * m_factor
		(((m_inverseDistance) * m0) * m1) * m_factor

		//m_inverseDistance * m_factor
	} * step.GetTime();

	const auto P = m_impulse * m_u;
	const auto LA = Cross(m_rA, -P) / Radian;
	const auto LB = Cross(m_rB, -P) / Radian;
	velA -= Velocity{invMassA * -P, invRotInertiaA * LA};
	velB += Velocity{invMassB * -P, invRotInertiaB * LB};

	bodyConstraintA->SetVelocity(velA);
	bodyConstraintB->SetVelocity(velB);
	return oldImpulse == m_impulse;
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
	return m_impulse * m_u;
}

AngularMomentum GravityJoint::GetAngularReaction() const
{
	return AngularMomentum{0};
}

} // namespace d2
} // namespace playrho
