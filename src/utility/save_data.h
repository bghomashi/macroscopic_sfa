#pragma once

#include <string>
#include "maths/maths.h"
#include "maths/vec2.h"

bool Store(const std::string& filename, const dvector& domain, const c2vector& data);
bool Store(const std::string& filename, const dvector& domain, const d2vector& data);
bool Store(const std::string& filename, const dvector& domain, const cvector& data);
bool Store(const std::string& filename, const dvector& domain, const dvector& data);

bool Store(const std::string& filename, const c2vector& data);
bool Store(const std::string& filename, const d2vector& data);
bool Store(const std::string& filename, const cvector& data);
bool Store(const std::string& filename, const dvector& data);

bool StoreBinary(const std::string& filename, const c2vector& data, int index = -1);
bool StoreBinary(const std::string& filename, const dvector& data);


