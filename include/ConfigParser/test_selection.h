#pragma once
#include <variant>
#include "../../include/GMGPolar/test_cases.h"

using DomainGeometryVariant = std::variant<CircularGeometry, ShafranovGeometry, CzarnyGeometry, CulhamGeometry>;
