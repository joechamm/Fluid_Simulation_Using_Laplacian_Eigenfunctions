#pragma once

#include <commonFramework/include/scene/Scene.h>
#include <commonFramework/include/scene/Material.h>
#include <commonFramework/include/scene/VtxData.h>
#include <glFramework/include/GLShader.h>
#include <glFramework/include/GLTexture.h>
#include <glFramework/include/kIdxBinds.h>

namespace glFramework
{
	using commonFramework::MeshFileHeader;
	using commonFramework::MeshData;
	using commonFramework::Scene;
	using commonFramework::MaterialDescription;
	using commonFramework::DrawData;

	class GLSceneData
	{
	public:
		GLSceneData(
			const char* meshFile,
			const char* sceneFile,
			const char* materialFile);

		std::vector<GLTexture> allMaterialTextures_;

		MeshFileHeader header_;
		MeshData meshData_;

		Scene scene_;
		std::vector<MaterialDescription> materials_;
		std::vector<DrawData> shapes_;

		void loadScene(const char* sceneFile);
	};
}