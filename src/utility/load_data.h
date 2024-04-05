#pragma once

#include <string>
#include "maths/maths.h"
#include "maths/vec2.h"

bool LoadBinary(const std::string& filename, c2vector& data, int index = -1);
bool LoadBinary(const std::string& filename, dvector& data);
