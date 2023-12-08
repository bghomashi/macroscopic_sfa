#pragma once

#include <string>
#include "maths/maths.h"
#include "maths/vec2.h"

bool Store(const std::string& filename, const dvector& domain, const c2vector& data);
bool Store(const std::string& filename, const dvector& domain, const d2vector& data);
bool Store(const std::string& filename, const dvector& domain, const cvector& data);
bool Store(const std::string& filename, const dvector& domain, const dvector& data);



