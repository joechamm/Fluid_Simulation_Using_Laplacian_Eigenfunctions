#pragma once

#include "Bitmap.h"

Bitmap convertEquirectangularMapToVerticalCross(const Bitmap& equirectangularMap);
Bitmap convertVerticalCrossToCubeMapFaces(const Bitmap& b);