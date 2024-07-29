#pragma once

#include "Bitmap.h"

Bitmap convertEquirectangularMapToVerticalCross(const Bitmap& equirectangularMap, int faceSize);
Bitmap convertVerticalCrossToCubeMapFaces(const Bitmap& b);