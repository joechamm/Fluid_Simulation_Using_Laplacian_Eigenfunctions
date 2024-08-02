#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <commonFramework/include/scene/vec4.h>

using std::vector;
using std::string;

using glm::mat4;
using glm::vec4;

namespace commonFramework
{

	enum MaterialFlags
	{
		sMaterialFlags_CastShadow = 0x1,
		sMaterialFlags_ReceiveShadow = 0x2,
		sMaterialFlags_Transparent = 0x4,
	};

	constexpr const uint64_t INVALID_TEXTURE = 0xFFFFFFFF;

	/**
	*	@struct						MaterialDescription
	*	@brief						The material description contains both the numeric values that define the lighting properties of the material and the set of texture indices for use in Vulkan or OpenGL. Since OpenGL texture handles are an opaque 64-bit handle, we use that type rather than a 32-bit integer. For empty textures a special guard value is used with all bits set to 1.
	*	@var emissiveColor_			for light emitting materials
	*	@var albedoColor_			ambient material color
	*	@var roughness_				roughness factor can use x or both x and y components for anisotropic roughness
	*	@var transparencyFactor_	for use with alpha-blended materials
	*	@var alphaTest_				threshold for using simple-punch through transparency rendering
	*	@var flags_					placeholder for differentiating different types of material properties
	*	@var ambientOcclusionMap_	ambient occlusion texture handle
	*	@var emissiveMap_			emissive texture handle
	*	@var albedoMap_				albedo texture handle
	*	@var metallicRoughnessMap_	metallic roughness texture handle
	*	@var normalMap_				normal map texture handle
	*	@var opacityMap_			used during conversion
	*/
	struct PACKED_STRUCT MaterialDescription final
	{
		gpuvec4 emissiveColor_ = { 0.5f, 0.5f, 0.5f, 0.5f };
		gpuvec4 albedoColor_ = { 1.5f, 1.0f, 1.0f, 1.5f };

		gpuvec4 roughness_ = { 1.5f, 1.50f, 0.5f, 5.0f };
		float transparencyFactor_ = 1.0f;
		float alphaTest_ = 0.0f;
		float metallicFactor_ = 1.0f;

		uint32_t flags_ = sMaterialFlags_CastShadow | sMaterialFlags_ReceiveShadow;

		uint64_t ambientOcclusionMap_ = INVALID_TEXTURE;
		uint64_t emissiveMap_ = INVALID_TEXTURE;
		uint64_t albedoMap_ = INVALID_TEXTURE;
		uint64_t metallicRoughnessMap_ = INVALID_TEXTURE;
		uint64_t normalMap_ = INVALID_TEXTURE;
		uint64_t opacityMap_ = INVALID_TEXTURE;
	};

	static_assert(sizeof(MaterialDescription) % 16 == 0, "MaterialDescription should be padded to 16 bytes");

	/**
	 * @brief Save the converted material data to file.
	 * @param filename name of the file to save the material data.
	 * @param materials array of material data to be saved
	 * @param files names of the texture files associated the material data
	*/
	void saveMaterials(const char* filename, const vector<MaterialDescription>& materials, const vector<string>& files);

	/**
	 * @brief Load material data from a file.
	 * @param filename name of the file where the material data is stored.
	 * @param materials array of material data to load the file contents into.
	 * @param files array to load the texture filenames into.
	*/
	void loadMaterials(const char* filename, vector<MaterialDescription>& materials, vector<string>& files);

	/**
	 * @brief Merge material lists from multiple scenes following the logic of merging in mergeScenes
	 * @param oldMaterials vector of pointers to all the materials from the input scenes
	 * @param oldTextures vector of pointers to all the texture files from the input scenes
	 * @param allMaterials output vector of all the material lists from all the input scenes
	 * @param newTextures out vector of all the texture filenames from all the input scenes
	*/
	void mergeMaterialLists(
		const vector<vector<MaterialDescription>* >& oldMaterials,
		const vector<vector<string>* >& oldTextures,
		vector<MaterialDescription>& allMaterials,
		vector<string>& newTextures
	);

}