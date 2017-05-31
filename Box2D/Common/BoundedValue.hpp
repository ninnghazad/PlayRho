/*
 * Copyright (c) 2017 Louis Langholtz https://github.com/louis-langholtz/Box2D
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
 *    claim that you wrote the original software. If you use this software
 *    in a product, an acknowledgment in the product documentation would be
 *    appreciated but is not required.
 * 2. Altered source versions must be plainly marked as such, and must not be
 *    misrepresented as being the original software.
 * 3. This notice may not be removed or altered from any source distribution.
 */

#ifndef BoundedValue_hpp
#define BoundedValue_hpp

#include <stdexcept>
#include <limits>

namespace box2d {
    
    enum class LoValueCheck
    {
        Any,
        AboveZero,
        ZeroOrMore,
        AboveNegInf
    };

    enum class HiValueCheck
    {
        Any,
        BelowZero,
    	ZeroOrLess,
        BelowPosInf
    };

    template <typename T, LoValueCheck lo, HiValueCheck hi>
    class BoundedValue
    {
    public:
        using value_type = T;
        using exception_type = std::invalid_argument;
        using this_type = BoundedValue<value_type, lo, hi>;

        static constexpr LoValueCheck GetLoCheck() { return lo; }
        static constexpr HiValueCheck GetHiCheck() { return hi; }

        static constexpr void DoLoCheck(value_type value)
        {
            switch (GetLoCheck())
            {
                case LoValueCheck::Any:
                    return;
                case LoValueCheck::AboveZero:
                    if (!(value > value_type(0)))
                    {
                        throw exception_type{"value not > 0"};
                    }
                    return;
                case LoValueCheck::ZeroOrMore:
                    if (!(value >= value_type(0)))
                    {
                        throw exception_type{"value not >= 0"};
                    }
                    return;
                case LoValueCheck::AboveNegInf:
                    if (std::numeric_limits<value_type>::has_infinity)
                    {
	                    if (!(value > -std::numeric_limits<value_type>::infinity()))
                        {
                            throw exception_type{"value not > -inf"};;
                        }
                    }
                    return;
            }
        }
        
        static constexpr void DoHiCheck(value_type value)
        {
            switch (GetHiCheck())
            {
                case HiValueCheck::Any:
                    return;
                case HiValueCheck::BelowZero:
                    if (!(value < value_type(0)))
                    {
                        throw exception_type{"value not < 0"};
                    }
                    return;
                case HiValueCheck::ZeroOrLess:
                    if (!(value <= value_type(0)))
                    {
                        throw exception_type{"value not <= 0"};
                    }
                    return;
                case HiValueCheck::BelowPosInf:
                    if (std::numeric_limits<value_type>::has_infinity)
                    {
                        if (!(value < +std::numeric_limits<value_type>::infinity()))
                        {
                            throw exception_type{"value not < +inf"};;
                        }
                    }
                    return;
            }
        }

        constexpr BoundedValue(value_type value): m_value{value}
        {
            DoLoCheck(value);
            DoHiCheck(value);
        }
        
        constexpr BoundedValue& operator= (const this_type& other) noexcept
        {
            m_value = other.m_value;
            return *this;
        }

        constexpr BoundedValue& operator= (const T& value)
        {
            DoLoCheck(value);
            DoHiCheck(value);
            m_value = value;
            return *this;
        }

        constexpr operator value_type () const
        {
            return m_value;
        }

    private:
        value_type m_value;
    };
    
    // Logical operations for BoundedValue<T, lo, hi> OP BoundedValue<T, lo, hi>

    template <typename T, LoValueCheck lo, HiValueCheck hi>
    bool operator== (const BoundedValue<T, lo, hi> lhs, const BoundedValue<T, lo, hi> rhs)
    {
        return T{lhs} == T{rhs};
    }
    
    template <typename T, LoValueCheck lo, HiValueCheck hi>
    bool operator!= (const BoundedValue<T, lo, hi> lhs, const BoundedValue<T, lo, hi> rhs)
    {
        return T{lhs} != T{rhs};
    }
    
    template <typename T, LoValueCheck lo, HiValueCheck hi>
    bool operator<= (const BoundedValue<T, lo, hi> lhs, const BoundedValue<T, lo, hi> rhs)
    {
        return T{lhs} <= T{rhs};
    }
    
    template <typename T, LoValueCheck lo, HiValueCheck hi>
    bool operator>= (const BoundedValue<T, lo, hi> lhs, const BoundedValue<T, lo, hi> rhs)
    {
        return T{lhs} >= T{rhs};
    }
    
    template <typename T, LoValueCheck lo, HiValueCheck hi>
    bool operator< (const BoundedValue<T, lo, hi> lhs, const BoundedValue<T, lo, hi> rhs)
    {
        return T{lhs} < T{rhs};
    }
    
    template <typename T, LoValueCheck lo, HiValueCheck hi>
    bool operator> (const BoundedValue<T, lo, hi> lhs, const BoundedValue<T, lo, hi> rhs)
    {
        return T{lhs} > T{rhs};
    }

    // Logical operations for BoundedValue<T, lo, hi> OP T

    template <typename T, LoValueCheck lo, HiValueCheck hi>
    bool operator== (const BoundedValue<T, lo, hi> lhs, const T rhs)
    {
        return T{lhs} == T{rhs};
    }
    
    template <typename T, LoValueCheck lo, HiValueCheck hi>
    bool operator!= (const BoundedValue<T, lo, hi> lhs, const T rhs)
    {
        return T{lhs} != T{rhs};
    }
    
    template <typename T, LoValueCheck lo, HiValueCheck hi>
    bool operator<= (const BoundedValue<T, lo, hi> lhs, const T rhs)
    {
        return T{lhs} <= T{rhs};
    }
    
    template <typename T, LoValueCheck lo, HiValueCheck hi>
    bool operator>= (const BoundedValue<T, lo, hi> lhs, const T rhs)
    {
        return T{lhs} >= T{rhs};
    }
    
    template <typename T, LoValueCheck lo, HiValueCheck hi>
    bool operator< (const BoundedValue<T, lo, hi> lhs, const T rhs)
    {
        return T{lhs} < T{rhs};
    }
    
    template <typename T, LoValueCheck lo, HiValueCheck hi>
    bool operator> (const BoundedValue<T, lo, hi> lhs, const T rhs)
    {
        return T{lhs} > T{rhs};
    }

    // Logical operations for T OP BoundedValue<T, lo, hi>

    template <typename T, LoValueCheck lo, HiValueCheck hi>
    bool operator== (const T lhs, const BoundedValue<T, lo, hi> rhs)
    {
        return T{lhs} == T{rhs};
    }
    
    template <typename T, LoValueCheck lo, HiValueCheck hi>
    bool operator!= (const T lhs, const BoundedValue<T, lo, hi> rhs)
    {
        return T{lhs} != T{rhs};
    }
    
    template <typename T, LoValueCheck lo, HiValueCheck hi>
    bool operator<= (const T lhs, const BoundedValue<T, lo, hi> rhs)
    {
        return T{lhs} <= T{rhs};
    }
    
    template <typename T, LoValueCheck lo, HiValueCheck hi>
    bool operator>= (const T lhs, const BoundedValue<T, lo, hi> rhs)
    {
        return T{lhs} >= T{rhs};
    }
    
    template <typename T, LoValueCheck lo, HiValueCheck hi>
    bool operator< (const T lhs, const BoundedValue<T, lo, hi> rhs)
    {
        return T{lhs} < T{rhs};
    }
    
    template <typename T, LoValueCheck lo, HiValueCheck hi>
    bool operator> (const T lhs, const BoundedValue<T, lo, hi> rhs)
    {
        return T{lhs} > T{rhs};
    }

    // Common useful aliases...

    template <typename T>
    using NonNegative = BoundedValue<T, LoValueCheck::ZeroOrMore, HiValueCheck::Any>;

    template <typename T>
    using NonPositive = BoundedValue<T, LoValueCheck::Any, HiValueCheck::ZeroOrLess>;

    template <typename T>
    using Positive = BoundedValue<T, LoValueCheck::AboveZero, HiValueCheck::Any>;

    template <typename T>
    using Negative = BoundedValue<T, LoValueCheck::Any, HiValueCheck::BelowZero>;

    template <typename T>
    using Finite = BoundedValue<T, LoValueCheck::AboveNegInf, HiValueCheck::BelowPosInf>;
}

#endif /* BoundedValue_hpp */