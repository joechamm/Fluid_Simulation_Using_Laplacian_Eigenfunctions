#pragma once

#include <commonFramework/include/scene/Scene.h>
#include <commonFramework/include/scene/VtxData.h>

namespace commonFramework
{
	void mergeScene(Scene& scene, MeshData& meshData, const std::string& materialName);
}