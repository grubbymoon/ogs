/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MESHLIB_ELEMENTS_UTILS_H_
#define MESHLIB_ELEMENTS_UTILS_H_

#include <algorithm>
#include <vector>

#include "Elements.h"
#include "MeshLib/Node.h"

namespace
{
template <typename Element>
struct LowerOrderElement;

template <> struct LowerOrderElement<MeshLib::Hex>
{
    using type = void;
};
template <> struct LowerOrderElement<MeshLib::Prism>
{
    using type = void;
};
template <> struct LowerOrderElement<MeshLib::Pyramid>
{
    using type = void;
};
template <> struct LowerOrderElement<MeshLib::Quad>
{
    using type = void;
};
template <> struct LowerOrderElement<MeshLib::Tri>
{
    using type = void;
};
template <> struct LowerOrderElement<MeshLib::Point>
{
    using type = void;
};
template <> struct LowerOrderElement<MeshLib::Line>
{
    using type = void;
};

template <> struct LowerOrderElement<MeshLib::Hex20>
{
    using type = MeshLib::Hex;
};
template <> struct LowerOrderElement<MeshLib::Prism15>
{
    using type = MeshLib::Prism;
};
template <> struct LowerOrderElement<MeshLib::Pyramid13>
{
    using type = MeshLib::Pyramid;
};
template <> struct LowerOrderElement<MeshLib::Tet10>
{
    using type = MeshLib::Tet;
};
template <> struct LowerOrderElement<MeshLib::Quad8>
{
    using type = MeshLib::Quad;
};
template <> struct LowerOrderElement<MeshLib::Tri6>
{
    using type = MeshLib::Tri;
};
template <> struct LowerOrderElement<MeshLib::Line3>
{
    using type = MeshLib::Line;
};
}  // anonymous namespace

namespace MeshLib
{
/// Returns a vector of node pointers containing the base nodes of the elements
/// input vector.
inline std::vector<Node*> getBaseNodes(std::vector<Element*> const& elements)
{
    std::vector<Node*> base_nodes;
    base_nodes.reserve(elements.size() * 2);

    for (auto* const e : elements)
    {
        std::copy(e->getNodes(), e->getNodes() + e->getNumberOfBaseNodes(),
                  std::back_inserter(base_nodes));
    }

    return base_nodes;
}

template <typename ElementType>
Element* doSomething(ElementType const& element, int const order)
{
    using ElementTraits = LowerOrderElement<ElementType>;
    // The construction of the new element relies on the same nodes
    // order in the element and the lower order element.
    // TODO create deep copy of Nodes here.
    return new typename ElementTraits::type(
        const_cast<MeshLib::Node**>(element.getNodes()), element.getID());
}

/// Clones an element of an order p as a new element of given order r, where 0 <
/// r <= p. This cloning behaviour is used to create elements with lower order
/// (if requested) on meshes of higher (or equal) order. Given, for example, a
/// mesh with all elements being second order (Tri6, Quad8, etc.) lower order
/// elements (Tri3, Quad4, etc) are constructed.
/// If it is not possible to construct an element with requested order a nullptr
/// is returned.
inline Element* cloneWithGivenOrder(Element const* element, int const order)
{
    assert(order > 0);
    if (element->getElementsOrder() == order)  // same order, just clone.
        return element->clone();
    if (element->getElementsOrder() - 1 == order)  // one order lower
    {
        if (auto* const e = dynamic_cast<MeshLib::Hex20 const*>(element))
        {
            return doSomething(*e, order);
        }
        else if (auto* const e = dynamic_cast<MeshLib::Prism15 const*>(element))
        {
            return doSomething(*e, order);
        }
        else if (auto* const e =
                     dynamic_cast<MeshLib::Pyramid13 const*>(element))
        {
            return doSomething(*e, order);
        }
        else if (auto* const e = dynamic_cast<MeshLib::Tet10 const*>(element))
        {
            return doSomething(*e, order);
        }
        else if (auto* const e = dynamic_cast<MeshLib::Quad8 const*>(element))
        {
            return doSomething(*e, order);
        }
        else if (auto* const e = dynamic_cast<MeshLib::Tri6 const*>(element))
        {
            return doSomething(*e, order);
        }
        else if (auto* const e = dynamic_cast<MeshLib::Line3 const*>(element))
        {
            return doSomething(*e, order);
        }
    }

    // otherwise the implementation is not supporting this.
    OGS_FATAL("Not implemented.");
    return nullptr;
}
}  // namespace MeshLib

#endif  // MESHLIB_ELEMENTS_UTILS_H_
