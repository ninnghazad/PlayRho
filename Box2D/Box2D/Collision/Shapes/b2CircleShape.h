/*
* Copyright (c) 2006-2009 Erin Catto http://www.box2d.org
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

#ifndef B2_CIRCLE_SHAPE_H
#define B2_CIRCLE_SHAPE_H

#include <Box2D/Collision/Shapes/b2Shape.h>

/// A circle shape.
class b2CircleShape : public b2Shape
{
public:
	constexpr b2CircleShape(): b2Shape(e_circle, 0.0f) {}

	constexpr explicit b2CircleShape(float32 radius, const b2Vec2& position) noexcept:
		b2Shape(e_circle, radius), m_p(position) {}

	b2CircleShape(const b2CircleShape&) = default;

	/// Implement b2Shape.
	b2Shape* Clone(b2BlockAllocator* allocator) const override;

	/// @see b2Shape::GetChildCount
	int32 GetChildCount() const override;

	/// Implement b2Shape.
	bool TestPoint(const b2Transform& transform, const b2Vec2& p) const override;

	/// Implement b2Shape.
	bool RayCast(b2RayCastOutput* output, const b2RayCastInput& input,
				const b2Transform& transform, int32 childIndex) const override;

	/// @see b2Shape::ComputeAABB
	void ComputeAABB(b2AABB* aabb, const b2Transform& transform, int32 childIndex) const override;

	/// @see b2Shape::ComputeMass
	void ComputeMass(b2MassData* massData, float32 density) const override;

	b2Vec2 GetPosition() const noexcept { return m_p; }

	void SetPosition(const b2Vec2& value) noexcept { m_p = value; }

private:
	/// Position
	b2Vec2 m_p = b2Vec2_zero;
};

#endif
