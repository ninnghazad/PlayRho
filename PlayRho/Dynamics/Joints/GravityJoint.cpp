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
	m_frequency{def.frequency},
	m_dampingRatio{def.dampingRatio},
	m_maxForce{def.maxForce},
	m_radius{def.radius}
{
	//assert(IsValid(def.gravity));
	assert(IsValid(def.dampingRatio));
}

void GravityJoint::Accept(JointVisitor& visitor) const
{
	visitor.Visit(*this);
}

void GravityJoint::Accept(JointVisitor& visitor)
{
	visitor.Visit(*this);
}

Mass22 GravityJoint::GetEffectiveMassMatrix(const BodyConstraint& body) const noexcept
{
	// K	= [(1/m1 + 1/m2) * eye(2) - skew(r1) * invI1 * skew(r1) - skew(r2) * invI2 * skew(r2)]
	//	  = [1/m1+1/m2	 0	] + invI1 * [r1.y*r1.y -r1.x*r1.y] + invI2 * [r1.y*r1.y -r1.x*r1.y]
	//		[	0	 1/m1+1/m2]		   [-r1.x*r1.y r1.x*r1.x]		   [-r1.x*r1.y r1.x*r1.x]
/*
	const auto invMass = body.GetInvMass();
	const auto invRotInertia = body.GetInvRotInertia();

	const auto exx = InvMass{invMass + (invRotInertia * Square(GetY(m_rB)) / SquareRadian) + m_gamma};
	const auto exy = InvMass{-invRotInertia * GetX(m_rB) * GetY(m_rB) / SquareRadian};
	const auto eyy = InvMass{invMass + (invRotInertia * Square(GetX(m_rB)) / SquareRadian) + m_gamma};

	InvMass22 K;
	GetX(GetX(K)) = exx;
	GetY(GetX(K)) = exy;
	GetX(GetY(K)) = exy;
	GetY(GetY(K)) = eyy;
	return Invert(K);
*/
	return {};
}

void GravityJoint::InitVelocityConstraints(
	BodyConstraintsMap& bodies,
	const StepConf& step,
	const ConstraintSolverConf&
) {
	auto& bodyConstraintA = At(bodies, GetBodyA());
	auto& bodyConstraintB = At(bodies, GetBodyB());

	const auto invMassA = bodyConstraintA->GetInvMass();
	const auto invRotInertiaA = bodyConstraintA->GetInvRotInertia(); // L^-2 M^-1 QP^2
	const auto invMassB = bodyConstraintB->GetInvMass();
	const auto invRotInertiaB = bodyConstraintB->GetInvRotInertia(); // L^-2 M^-1 QP^2

	const auto posA = bodyConstraintA->GetPosition();
	auto velA = bodyConstraintA->GetVelocity();

	const auto posB = bodyConstraintB->GetPosition();
	auto velB = bodyConstraintB->GetVelocity();

	const auto qA = UnitVec::Get(posA.angular);
	const auto qB = UnitVec::Get(posB.angular);

	m_rA = Rotate(Length2{} - bodyConstraintA->GetLocalCenter(), qA);
	m_rB = Rotate(Length2{} - bodyConstraintB->GetLocalCenter(), qB);
	const auto deltaLocation = Length2{(posB.linear + m_rB) - (posA.linear + m_rA)};

	//const auto deltaLocation = Length2{posB.linear - posA.linear};

	// Handle singularity.
	const auto uvresult = UnitVec::Get(deltaLocation[0], deltaLocation[1]);
	m_u = std::get<UnitVec>(uvresult);
	const auto length = std::min(std::get<Length>(uvresult),m_radius);
	//const auto length = std::get<Length>(uvresult);

	const auto crAu = Length{Cross(m_rA, m_u)} / Radian;
	const auto crBu = Length{Cross(m_rB, m_u)} / Radian;
	const auto invRotMassA = invRotInertiaA * Square(crAu);
	const auto invRotMassB = invRotInertiaB * Square(crBu);
	auto invMass = InvMass{invMassA + invRotMassA + invMassB + invRotMassB};

	// Compute the effective mass matrix.
	m_mass = (invMass != InvMass{0}) ? Real{1} / invMass: 0_kg;

	//const auto C = length - m_radius; // L
	const auto C = m_radius - length; // L

	// Frequency
	const auto omega = Real{2} * Pi * m_frequency;

	// Damping coefficient
	const auto d = Real{2} * m_mass * m_dampingRatio * omega; // M T^-1

	// Spring stiffness
	const auto k = m_mass * Square(omega); // M T^-2

	// magic formulas
	const auto h = step.GetTime();
	const auto gamma = Mass{h * (d + h * k)}; // T (M T^-1 + T M T^-2) = M
	m_invGamma = (gamma != 0_kg)? Real{1} / gamma: 0;
	m_bias = C * h * k * m_invGamma; // L T M T^-2 M^-1 = L T^-1

	invMass += m_invGamma;
	m_mass = (invMass != InvMass{0}) ? Real{1} / invMass: 0;

	if (step.doWarmStart)
	{
		// Scale the impulse to support a variable time step.
		m_impulse *= step.dtRatio;

		const auto P = m_impulse * m_u;

		// P is M L T^-2
		// Cross(Length2, P) is: M L^2 T^-1
		// inv rotational inertia is: L^-2 M^-1 QP^2
		// Product is: L^-2 M^-1 QP^2 M L^2 T^-1 = QP^2 T^-1
		const auto LA = Cross(m_rA, P) / Radian;
		const auto LB = Cross(m_rB, P) / Radian;
		velA -= Velocity{invMassA * P, invRotInertiaA * LA};
		velB += Velocity{invMassB * P, invRotInertiaB * LB};
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

	// Cdot = dot(u, v + cross(w, r))
	const auto vpA = velA.linear + GetRevPerpendicular(m_rA) * (velA.angular / Radian);
	const auto vpB = velB.linear + GetRevPerpendicular(m_rB) * (velB.angular / Radian);
	const auto Cdot = LinearVelocity{Dot(m_u, vpB - vpA)};

	const auto impulse = Momentum{-m_mass * (Cdot + m_bias + m_invGamma * m_impulse)};
	m_impulse += impulse;

	const auto P = impulse * m_u;
	const auto LA = Cross(m_rA, P) / Radian;
	const auto LB = Cross(m_rB, P) / Radian;
	velA -= Velocity{invMassA * P, invRotInertiaA * LA};
	velB += Velocity{invMassB * P, invRotInertiaB * LB};

	bodyConstraintA->SetVelocity(velA);
	bodyConstraintB->SetVelocity(velB);

	return impulse == 0_Ns;
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
