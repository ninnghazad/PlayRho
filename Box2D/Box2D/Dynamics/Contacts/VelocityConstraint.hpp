/*
 * Original work Copyright (c) 2006-2009 Erin Catto http://www.box2d.org
 * Modified work Copyright (c) 2016 Louis Langholtz https://github.com/louis-langholtz/Box2D
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

#ifndef VelocityConstraint_hpp
#define VelocityConstraint_hpp

#include <Box2D/Common/Math.h>

// Define <code>BOX2D_CACHE_VC_POINT_MASSES</code> to cache velocity constraint point masses
// instead of re-computing them every time these values are read. This uses an additional
// 16-bytes of memory per VelocityConstraint object.
//
// Note that with up to 4000 elements in the Tumbler test and using a library built without
// optimizations enabled, caching the masses does seem to result in faster simulations.
// It's unknown whether increasing the number of elements would eventually result in it being
// faster not to cache the mass values.
//
#define BOX2D_CACHE_VC_POINT_MASSES

namespace box2d {
	
	/// Contact velocity constraint.
	/// @note A valid contact velocity constraint must have a point count of either 1 or 2.
	/// @note This data structure is 176-bytes large (on at least one 64-bit platform) if
	///   <code>BOX2D_CACHE_VC_POINT_MASSES</code> is defined. It's 160-bytes large otherwise.
	class VelocityConstraint
	{
	public:
		using size_type = std::remove_const<decltype(MaxManifoldPoints)>::type;
		using index_type = size_t;

		/// Velocity constraint point.
		/// @note This structure is at least 36-bytes large.
		struct Point
		{
			Vec2 rA; ///< Position of body A relative to world manifold point (8-bytes).
			Vec2 rB; ///< Position of body B relative to world manifold point (8-bytes).
			float_t normalImpulse; ///< Normal impulse (4-bytes).
			float_t tangentImpulse; ///< Tangent impulse (4-bytes).
#if defined(BOX2D_CACHE_VC_POINT_MASSES)
			float_t normalMass; ///< Normal mass (4-bytes). 0 or greater.
			float_t tangentMass; ///< Tangent mass (4-bytes).
#endif
			float_t velocityBias; ///< Velocity bias (4-bytes).
		};

		/// Contact velocity constraint body data.
		/// @invariant The inverse mass is a value of zero or more.
		/// @invariant The inverse rotational inertia is a value of zero or more.
		class BodyData
		{
		public:
			BodyData() noexcept = default;
			BodyData(const BodyData& copy) noexcept = default;
			
			/// Initializing constructor.
			/// @note Behavior is undefined if the given inverse mass or given inverse rotational
			///   inertia is less than zero.
			/// @param iM Inverse mass. A value of 0 or more.
			/// @param iI Inverse rotational inertia. A value of 0 or more.
			constexpr BodyData(index_type i, float_t iM, float_t iI) noexcept:
				index{i}, invMass{iM}, invI{iI}
			{
				assert(iM >= 0);
				assert(iI >= 0);
			}
			
			index_type GetIndex() const noexcept { return index; }
			
			/// Gets the inverse mass.
			/// @return 0 or greater value.
			float_t GetInvMass() const noexcept { return invMass; }
			
			/// Gets the inverse rotational inertia.
			/// @return 0 or greater value.
			float_t GetInvRotI() const noexcept { return invI; }
			
		private:
			float_t invMass = 0; ///< Inverse mass of body. Value of 0 or greater.
			float_t invI = 0; ///< Inverse rotational interia of body. Value of 0 or greater.
			index_type index = 0; ///< Index within island of body.
		};
		
		/// Default constructor.
		/// @detail
		/// Initializes object with: a zero point count, an invalid K, an invalid normal mass,
		/// an invalid normal, invalid friction, invalid restitution, an invalid tangent speed.
		VelocityConstraint() = default;
		
		VelocityConstraint(const VelocityConstraint& copy) = default;

		VelocityConstraint& operator= (const VelocityConstraint& copy) = default;
		
		VelocityConstraint(index_type ci, float_t f, float_t r, float_t ts):
			contactIndex{ci}, friction{f}, restitution{r}, tangentSpeed{ts} {}
		
		/// Gets the count of points added to this object.
		/// @return Value between 0 and MaxManifoldPoints
		/// @sa MaxManifoldPoints.
		/// @sa AddPoint.
		size_type GetPointCount() const noexcept { return pointCount; }
		
		/// Accesses the point identified by the given index.
		/// @note Behavior is undefined if the identified point does not exist.
		/// @param index Index of the point to return. This is a value less than returned by GetPointCount().
		/// @return velocity constraint point for the given index.
		/// @sa GetPointCount.
		const Point& PointAt(size_type index) const
		{
			assert(index < pointCount);
			return points[index];
		}
		
		/// Accesses the point identified by the given index.
		/// @note Behavior is undefined if the identified point does not exist.
		/// @param index Index of the point to return. This is a value less than returned by GetPointCount().
		/// @return velocity constraint point for the given index.
		/// @sa GetPointCount.
		struct Point& PointAt(size_type index)
		{
			assert(index < pointCount);
			return points[index];
		}
		
		void ClearPoints() noexcept { pointCount = 0; }
		
		/// Adds the given point to this contact velocity constraint object.
		/// @detail Adds up to MaxManifoldPoints points. To find out how many points have already
		///   been added, call GetPointCount().
		/// @param val Velocity constraint point value to add.
		/// @note Behavior is undefined if an attempt is made to add more than MaxManifoldPoints points.
		/// @sa GetPointCount().
		void AddPoint(const Point& val);
		
		void RemovePoint() noexcept;
		
		/// Sets this object's K value.
		/// @param value A position constraint dependent value or the zero matrix (Mat22_zero).
		void SetK(const Mat22& value) noexcept;
		
		/// Gets the "K" value.
		/// @note This value is only valid if SetK had previously been called with a valid value.
		/// @warning Behavior is undefined if called before SetK was called.
		/// @return "K" value previously set.
		Mat22 GetK() const noexcept;
		
		Mat22 GetNormalMass() const noexcept;
		
		/// Gets the contact index.
		/// @note This value can only be set via the initializing constructor.
		/// @return Index of the associated contact (the index of the contact that this constraint is for).
		index_type GetContactIndex() const noexcept { return contactIndex; }
		
		/// Gets the combined friction of the associated contact.
		float_t GetFriction() const noexcept { return friction; }
		
		/// Gets the combined restitution of the associated contact.
		float_t GetRestitution() const noexcept { return restitution; }
		
		/// Gets the tangent speed of the associated contact.
		float_t GetTangentSpeed() const noexcept { return tangentSpeed; }
		
		Vec2 normal = GetInvalid<Vec2>(); ///< Normal of the world manifold.
		
		BodyData bodyA; ///< Body A contact velocity constraint data.
		BodyData bodyB; ///< Body B contact velocity constraint data.
		
	private:
		float_t friction = GetInvalid<float_t>(); ///< Friction coefficient (4-bytes). Usually in the range of [0,1].
		float_t restitution = GetInvalid<float_t>(); ///< Restitution coefficient (4-bytes).
		float_t tangentSpeed = GetInvalid<float_t>(); ///< Tangent speed (4-bytes).
		
		index_type contactIndex = GetInvalid<index_type>(); ///< Index of the contact that this constraint is for (typically 8-bytes).
		
		// K and normalMass fields are only used for the block solver.
		Mat22 K = GetInvalid<Mat22>(); ///< Block solver "K" info (only used by block solver, 16-bytes).
		Mat22 normalMass = GetInvalid<Mat22>(); ///< Block solver "normal mass" info (only used by block solver, 16-bytes).
		
		Point points[MaxManifoldPoints]; ///< Velocity constraint points array (at least 72-bytes).
		size_type pointCount = 0; ///< Point count (at least 1-byte).
	};
	
	inline void VelocityConstraint::AddPoint(const Point& val)
	{
		assert(pointCount < MaxManifoldPoints);
		points[pointCount] = val;
		++pointCount;
	}
	
	inline void VelocityConstraint::RemovePoint() noexcept
	{
		assert(pointCount > 0);
		--pointCount;
	}
	
	inline void VelocityConstraint::SetK(const Mat22& value) noexcept
	{
		assert(IsValid(value));
		K = value;
		normalMass = Invert(value);
	}
	
	/// Gets the "K" value.
	/// @return "K" value or the invalid Mat22 if no other value has been set.
	/// @sa SetK.
	inline Mat22 VelocityConstraint::GetK() const noexcept
	{
		return K;
	}
	
	/// Gets the "normal mass" value.
	/// @return "normal mass" value or the invalid Mat22 if no other value has been set.
	/// @sa SetK.
	inline Mat22 VelocityConstraint::GetNormalMass() const noexcept
	{
		return normalMass;
	}
	
	inline Vec2 GetNormal(const VelocityConstraint& vc)
	{
		return vc.normal;
	}

	inline Vec2 GetTangent(const VelocityConstraint& vc)
	{
		return GetFwdPerpendicular(vc.normal);
	}

	inline float_t GetInverseMass(const VelocityConstraint& vc)
	{
		return vc.bodyA.GetInvMass() + vc.bodyB.GetInvMass();
	}

	inline float_t ComputeNormalMassAtPoint(const VelocityConstraint& vc, VelocityConstraint::size_type index)
	{
		const auto value = GetInverseMass(vc)
			+ (vc.bodyA.GetInvRotI() * Square(Cross(vc.PointAt(index).rA, GetNormal(vc))))
			+ (vc.bodyB.GetInvRotI() * Square(Cross(vc.PointAt(index).rB, GetNormal(vc))));
		return (value != 0)? float_t{1} / value : float_t{0};
	}
	
	inline float_t ComputeTangentMassAtPoint(const VelocityConstraint& vc, VelocityConstraint::size_type index)
	{
		const auto value = GetInverseMass(vc)
			+ (vc.bodyA.GetInvRotI() * Square(Cross(vc.PointAt(index).rA, GetTangent(vc))))
			+ (vc.bodyB.GetInvRotI() * Square(Cross(vc.PointAt(index).rB, GetTangent(vc))));
		return (value != 0)? float_t{1} / value : float_t{0};
	}
	
	inline float_t GetNormalMassAtPoint(const VelocityConstraint& vc, VelocityConstraint::size_type index)
	{
#if defined(BOX2D_CACHE_VC_POINT_MASSES)
		return vc.PointAt(index).normalMass;
#else
		return ComputeNormalMassAtPoint(vc, index);
#endif
	}
	
	inline float_t GetTangentMassAtPoint(const VelocityConstraint& vc, VelocityConstraint::size_type index)
	{
#if defined(BOX2D_CACHE_VC_POINT_MASSES)
		return vc.PointAt(index).normalMass;
#else
		return ComputeTangentMassAtPoint(vc, index);
#endif
	}

	inline Vec2 GetNormalImpulses(const VelocityConstraint& vc)
	{
		return Vec2{vc.PointAt(0).normalImpulse, vc.PointAt(1).normalImpulse};
	}
	
	inline void SetNormalImpulses(VelocityConstraint& vc, const Vec2 impulses)
	{
		vc.PointAt(0).normalImpulse = impulses[0];
		vc.PointAt(1).normalImpulse = impulses[1];
	}
	
	inline Vec2 GetTangentImpulses(const VelocityConstraint& vc)
	{
		return Vec2{vc.PointAt(0).tangentImpulse, vc.PointAt(1).tangentImpulse};
	}

	inline void SetTangentImpulses(VelocityConstraint& vc, const Vec2 impulses)
	{
		vc.PointAt(0).tangentImpulse = impulses[0];
		vc.PointAt(1).tangentImpulse = impulses[1];
	}

} // namespace box2d

#endif /* VelocityConstraint_hpp */