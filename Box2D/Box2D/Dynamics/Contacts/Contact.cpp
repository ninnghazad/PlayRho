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

#include <Box2D/Dynamics/Contacts/Contact.hpp>
#include <Box2D/Collision/Collision.hpp>
#include <Box2D/Collision/Distance.hpp>
#include <Box2D/Collision/DistanceProxy.hpp>
#include <Box2D/Collision/Shapes/Shape.hpp>
#include <Box2D/Dynamics/Body.hpp>
#include <Box2D/Dynamics/Fixture.hpp>
#include <Box2D/Dynamics/World.hpp>
#include <Box2D/Dynamics/StepConf.hpp>

using namespace box2d;

using ContactCreateFcn = Contact* (Fixture* fixtureA, child_count_t indexA,
                                       Fixture* fixtureB, child_count_t indexB);
using ContactDestroyFcn = void (Contact* contact);

static inline Manifold::Conf GetManifoldConf(const StepConf& conf)
{
    auto manifoldConf = Manifold::Conf{};
    manifoldConf.tolerance = conf.tolerance;
    manifoldConf.targetDepth = conf.targetDepth;
    manifoldConf.maxCirclesRatio = conf.maxCirclesRatio;
    return manifoldConf;
}

static inline DistanceConf GetDistanceConf(const StepConf& conf)
{
    DistanceConf distanceConf;
    distanceConf.maxIterations = conf.maxDistanceIters;
    return distanceConf;
}

Contact* Contact::Create(Fixture& fixtureA, child_count_t indexA,
                         Fixture& fixtureB, child_count_t indexB)
{
    return new Contact{&fixtureA, indexA, &fixtureB, indexB};
}

void Contact::Destroy(Contact* contact)
{
    const auto fixtureA = contact->GetFixtureA();
    const auto fixtureB = contact->GetFixtureB();

    if ((contact->m_manifold.GetPointCount() > 0) &&
        !fixtureA->IsSensor() && !fixtureB->IsSensor())
    {
        // Contact may have been keeping accelerable bodies of fixture A or B from moving.
        // Need to awaken those bodies now in case they are again movable.
        SetAwake(*fixtureA);
        SetAwake(*fixtureB);
    }

    delete contact;
}

Contact::Contact(Fixture* fA, child_count_t iA, Fixture* fB, child_count_t iB):
    m_fixtureA{fA}, m_fixtureB{fB}, m_indexA{iA}, m_indexB{iB},
    m_friction{MixFriction(fA->GetFriction(), fB->GetFriction())},
    m_restitution{MixRestitution(fA->GetRestitution(), fB->GetRestitution())}
{
    assert(fA && fB);
    assert(fA != fB);
    assert(fA->GetShape());
    assert(fB->GetShape());
    assert(fA->GetBody() != fB->GetBody());
}

void Contact::Update(const StepConf& conf, ContactListener* listener)
{
    const auto oldManifold = m_manifold;

    // Note: do not assume the fixture AABBs are overlapping or are valid.
    const auto oldTouching = IsTouching();
    auto newTouching = false;

    const auto fixtureA = GetFixtureA();
    const auto indexA = GetChildIndexA();
    const auto fixtureB = GetFixtureB();
    const auto indexB = GetChildIndexB();
    const auto shapeA = fixtureA->GetShape();
    const auto xfA = fixtureA->GetBody()->GetTransformation();
    const auto shapeB = fixtureB->GetShape();
    const auto xfB = fixtureB->GetBody()->GetTransformation();
    const auto childA = shapeA->GetChild(indexA);
    const auto childB = shapeB->GetChild(indexB);

    // NOTE: The touching state returned by the TestOverlap function
    //    **SHOULD** agree with that returned from the CollideShapes function.

    const auto sensor = fixtureA->IsSensor() || fixtureB->IsSensor();
    if (sensor)
    {
        newTouching = TestOverlap(childA, xfA, childB, xfB, GetDistanceConf(conf));

        // Check agreement between TestOverlap and CollideShapes.
        assert(newTouching == (CollideShapes(childA, xfA, childB, xfB,
                                             GetManifoldConf(conf)).GetPointCount() > 0));
        // Sensors don't generate manifolds.
        m_manifold = Manifold{};
    }
    else
    {
        auto newManifold = CollideShapes(childA, xfA, childB, xfB, GetManifoldConf(conf));

        const auto old_point_count = oldManifold.GetPointCount();
        const auto new_point_count = newManifold.GetPointCount();

        newTouching = new_point_count > 0;

        // Check agreement between TestOverlap and CollideShapes.
        assert(newTouching == TestOverlap(childA, xfA, childB, xfB, GetDistanceConf(conf)));

        // Match old contact ids to new contact ids and copy the
        // stored impulses to warm start the solver.
        for (auto i = decltype(new_point_count){0}; i < new_point_count; ++i)
        {
            const auto new_cf = newManifold.GetContactFeature(i);
            for (auto j = decltype(old_point_count){0}; j < old_point_count; ++j)
            {
                if (new_cf == oldManifold.GetContactFeature(j))
                {
                    newManifold.SetContactImpulses(i, oldManifold.GetContactImpulses(j));
                    break;
                }
            }
        }

        // Ideally this method is **NEVER** called unless a dependency changed such
        // that the following assertion is **ALWAYS** valid.
        //assert(newManifold != oldManifold);

        m_manifold = newManifold;

#ifdef MAKE_CONTACT_PROCESSING_ORDER_DEPENDENT
        const auto bodyA = fixtureA->GetBody();
        const auto bodyB = fixtureB->GetBody();

        assert(bodyA);
        assert(bodyB);

        /*
         * The following code creates an ordering dependency in terms of update processing
         * over a container of contacts. It also puts this method into the situation of
         * modifying bodies which adds race potential in a multi-threaded mode of operation.
         * Lastly, without this code, the step-statistics show a world getting to sleep in
         * less TOI position iterations.
         */
        if (newTouching != oldTouching)
        {
            bodyA->SetAwake();
            bodyB->SetAwake();
        }
#endif
    }

    if (!oldTouching && newTouching)
    {
        SetTouching();
        if (listener)
        {
            listener->BeginContact(*this);
        }
    }
    else if (oldTouching && !newTouching)
    {
        UnsetTouching();
        if (listener)
        {
            listener->EndContact(*this);
        }
    }

    if (!sensor && newTouching)
    {
        if (listener)
        {
            listener->PreSolve(*this, oldManifold);
        }
    }
}

// Free functions...

bool box2d::HasSensor(const Contact& contact) noexcept
{
    return contact.GetFixtureA()->IsSensor() || contact.GetFixtureB()->IsSensor();
}

bool box2d::IsImpenetrable(const Contact& contact) noexcept
{
    const auto bA = contact.GetFixtureA()->GetBody();
    const auto bB = contact.GetFixtureB()->GetBody();
    return bA->IsImpenetrable() || bB->IsImpenetrable();
}

void box2d::SetAwake(Contact& c) noexcept
{
    SetAwake(*c.GetFixtureA());
    SetAwake(*c.GetFixtureB());
}

/// Resets the friction mixture to the default value.
void box2d::ResetFriction(Contact& contact)
{
    contact.SetFriction(MixFriction(contact.GetFixtureA()->GetFriction(), contact.GetFixtureB()->GetFriction()));
}

/// Reset the restitution to the default value.
void box2d::ResetRestitution(Contact& contact) noexcept
{
    const auto restitutionA = contact.GetFixtureA()->GetRestitution();
    const auto restitutionB = contact.GetFixtureB()->GetRestitution();
    contact.SetRestitution(MixRestitution(restitutionA, restitutionB));
}

TOIOutput box2d::CalcToi(const Contact& contact, const ToiConf conf)
{
    const auto fA = contact.GetFixtureA();
    const auto fB = contact.GetFixtureB();
    const auto bA = fA->GetBody();
    const auto bB = fB->GetBody();

    const auto proxyA = fA->GetShape()->GetChild(contact.GetChildIndexA());
    const auto proxyB = fB->GetShape()->GetChild(contact.GetChildIndexB());

    // Large rotations can make the root finder of TimeOfImpact fail, so normalize sweep angles.
    const auto sweepA = GetAnglesNormalized(bA->GetSweep());
    const auto sweepB = GetAnglesNormalized(bB->GetSweep());

    // Compute the TOI for this contact (one or both bodies are active and impenetrable).
    // Computes the time of impact in interval [0, 1]
    // Large rotations can make the root finder of TimeOfImpact fail, so normalize the sweep angles.
    return GetToiViaSat(proxyA, sweepA, proxyB, sweepB, conf);
}
