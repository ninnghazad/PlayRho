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

#ifndef COLLISION_PROCESSING_H
#define COLLISION_PROCESSING_H

#include <algorithm>

namespace box2d {

// This test shows collision processing and tests
// deferred body destruction.
class CollisionProcessing : public Test
{
public:
	CollisionProcessing()
	{
		// Ground body
		{
			EdgeShape shape;
			shape.Set(Vec2(-50.0f, 0.0f), Vec2(50.0f, 0.0f));

			FixtureDef sd;
			sd.shape = &shape;;

			BodyDef bd;
			Body* ground = m_world->CreateBody(bd);
			ground->CreateFixture(sd);
		}

		float_t xLo = -5.0f, xHi = 5.0f;
		float_t yLo = 2.0f, yHi = 35.0f;

		// Small triangle
		Vec2 vertices[3];
		vertices[0] = Vec2(-1.0f, 0.0f);
		vertices[1] = Vec2(1.0f, 0.0f);
		vertices[2] = Vec2(0.0f, 2.0f);

		PolygonShape polygon;
		polygon.Set(Span<const Vec2>{vertices, 3});

		FixtureDef triangleShapeDef;
		triangleShapeDef.shape = &polygon;
		triangleShapeDef.density = 1.0f;

		BodyDef triangleBodyDef;
		triangleBodyDef.type = BodyType::Dynamic;
		triangleBodyDef.position = Vec2(RandomFloat(xLo, xHi), RandomFloat(yLo, yHi));

		Body* body1 = m_world->CreateBody(triangleBodyDef);
		body1->CreateFixture(triangleShapeDef);

		// Large triangle (recycle definitions)
		vertices[0] *= 2.0f;
		vertices[1] *= 2.0f;
		vertices[2] *= 2.0f;
		polygon.Set(Span<const Vec2>{vertices, 3});

		triangleBodyDef.position = Vec2(RandomFloat(xLo, xHi), RandomFloat(yLo, yHi));

		Body* body2 = m_world->CreateBody(triangleBodyDef);
		body2->CreateFixture(triangleShapeDef);
		
		// Small box
		polygon.SetAsBox(1.0f, 0.5f);

		FixtureDef boxShapeDef;
		boxShapeDef.shape = &polygon;
		boxShapeDef.density = 1.0f;

		BodyDef boxBodyDef;
		boxBodyDef.type = BodyType::Dynamic;
		boxBodyDef.position = Vec2(RandomFloat(xLo, xHi), RandomFloat(yLo, yHi));

		Body* body3 = m_world->CreateBody(boxBodyDef);
		body3->CreateFixture(boxShapeDef);

		// Large box (recycle definitions)
		polygon.SetAsBox(2.0f, 1.0f);
		boxBodyDef.position = Vec2(RandomFloat(xLo, xHi), RandomFloat(yLo, yHi));
		
		Body* body4 = m_world->CreateBody(boxBodyDef);
		body4->CreateFixture(boxShapeDef);

		// Small circle
		CircleShape circle;
		circle.SetRadius(float_t{1});

		FixtureDef circleShapeDef;
		circleShapeDef.shape = &circle;
		circleShapeDef.density = 1.0f;

		BodyDef circleBodyDef;
		circleBodyDef.type = BodyType::Dynamic;
		circleBodyDef.position = Vec2(RandomFloat(xLo, xHi), RandomFloat(yLo, yHi));

		Body* body5 = m_world->CreateBody(circleBodyDef);
		body5->CreateFixture(circleShapeDef);

		// Large circle
		circle.SetRadius(circle.GetRadius() * 2);
		circleBodyDef.position = Vec2(RandomFloat(xLo, xHi), RandomFloat(yLo, yHi));

		Body* body6 = m_world->CreateBody(circleBodyDef);
		body6->CreateFixture(circleShapeDef);
	}

	void Step(Settings& settings, Drawer& drawer) override
	{
		Test::Step(settings, drawer);

		// We are going to destroy some bodies according to contact
		// points. We must buffer the bodies that should be destroyed
		// because they may belong to multiple contact points.
		const int32 k_maxNuke = 6;
		Body* nuke[k_maxNuke];
		int32 nukeCount = 0;

		// Traverse the contact results. Destroy bodies that
		// are touching heavier bodies.
		for (int32 i = 0; i < m_pointCount; ++i)
		{
			ContactPoint* point = m_points + i;

			Body* body1 = point->fixtureA->GetBody();
			Body* body2 = point->fixtureB->GetBody();
			const auto mass1 = GetMass(*body1);
			const auto mass2 = GetMass(*body2);

			if (mass1 > 0.0f && mass2 > 0.0f)
			{
				if (mass2 > mass1)
				{
					nuke[nukeCount++] = body1;
				}
				else
				{
					nuke[nukeCount++] = body2;
				}

				if (nukeCount == k_maxNuke)
				{
					break;
				}
			}
		}

		// Sort the nuke array to group duplicates.
		std::sort(nuke, nuke + nukeCount);

		// Destroy the bodies, skipping duplicates.
		int32 i = 0;
		while (i < nukeCount)
		{
			Body* b = nuke[i++];
			while (i < nukeCount && nuke[i] == b)
			{
				++i;
			}

			if (b != m_bomb)
			{
				m_world->Destroy(b);
			}
		}
	}

	static Test* Create()
	{
		return new CollisionProcessing;
	}
};

} // namespace box2d

#endif