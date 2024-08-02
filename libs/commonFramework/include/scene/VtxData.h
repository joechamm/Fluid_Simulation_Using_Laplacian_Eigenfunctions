#pragma once

#include <stdint.h>
#include <glm/glm.hpp>

#include <commonFramework/include/Utils.h>
#include <commonFramework/include/UtilsMath.h>

namespace commonFramework
{
	/**
 * @brief we need two constants to define the limits on how many LODs and vertex streams we can have in a single mesh
*/

	constexpr const uint32_t kMaxLODs = 8;
	constexpr const uint32_t kMaxStreams = 8;

	/**
	 * @brief Define an individual mesh description.
	 * We deliberately avoid using pointers that hide memory allocations and prohibit the simple saving and loading of data.
	 * We store offsets to individual data streams and LOD index buffers. They are equivalent to pointers but are more flexible and, most importantly, GPU-friendlier.
	 * All the offsets in the Mesh structure are given relative to the beginning of the data block.
	 * The LOD count, where the original mesh counts as one of the LODs, must be strictly less than kMaxLODs. This is because we do not store LOD index buffer sizes but calculate them from offsets.
	 * To calculate these sizes, we store one additional empty LOD level at the end. The number of vertex data streams is stored directly with no modifications.
	*/

	struct Mesh final
	{
		/* Number of LODs in this mesh. Strictly less than MAX_LODS, last LOD offset is used as a marker only */
		uint32_t lodCount = 1;

		/* Number of vertex data streams */
		uint32_t streamCount = 0;

		/* The total count of all previous vertices in this mesh file */
		uint32_t indexOffset = 0;

		uint32_t vertexOffset = 0;

		/* Vertex count (for all LODs) */
		uint32_t vertexCount = 0;

		/* Each mesh can potentially be displayed at different LODs. The file contains all the indices for all the LODs, and offsets to the beginning of each LOD are stored in the lodOffset array. This array contains one extra item at the end, which serves as a marker to calculate the size of the last LOD */
		uint32_t lodOffset[kMaxLODs] = { 0 };

		/* Instead of storing the sizes of each LOD, we define a little helper function to calculate their sizes */
		inline uint32_t getLODIndicesCount(uint32_t lod) const
		{
			return lodOffset[lod + 1] - lodOffset[lod];
		}

		/* streamOffset field stores offsets to all of the individual vertex data streams */
		/* IMPORTANT NOTE Besides the element size, we might want to store the element type, such as byte, short integer, or float. This information is important for performance reasons in real-world applications. To simplify the code in this book, we will not do it here */

		/* All the data "pointers" for all the streams */
		uint32_t streamOffset[kMaxStreams] = { 0 };
		/* Information about stream element (size pretty much defines everything else, the "semantics" is defined by the shader) */
		uint32_t streamElementSize[kMaxStreams] = { 0 };

		/*
		* We could have included the streamStride[] array here to allow interleaved storage of attributes. For this book we assume tightly-packed (non-interleaved) vertex attribute streams.
		*
		* Additional information, like mesh name, can be added here
		*/
	};

	struct MeshFileHeader
	{
		/* To ensure data integrity and to check the validity of the header, a magic hexadecimal value of 0x12345678 is stored in the first 4 bytes of the header */
		uint32_t magicValue;

		/* The number of different meshes in this file is stored in meshCount */
		uint32_t meshCount;

		/* For convenience, store an offset to the the beginning of the mesh data */
		uint32_t dataBlockStartOffset;

		/* store the sizes of index and vertex data in bytes to check the integrity of a mesh file */
		uint32_t indexDataSize;
		uint32_t vertexDataSize;
		/* Any additional meta data here*/
	};

	struct DrawData
	{
		uint32_t meshIndex;
		uint32_t materialIndex;
		uint32_t LOD;
		uint32_t indexOffset;
		uint32_t vertexOffset;
		uint32_t transformIndex;
	};

	struct MeshData
	{
		std::vector<uint32_t>		indexData_;
		std::vector<float>			vertexData_;
		std::vector<Mesh>			meshes_;
		std::vector<BoundingBox>	boxes_;
	};

	static_assert(sizeof(DrawData) == sizeof(uint32_t) * 6);
	static_assert(sizeof(BoundingBox) == sizeof(float) * 6);

	MeshFileHeader loadMeshData(const char* meshFilename, MeshData& out);
	void saveMeshData(const char* filename, const MeshData& m);

	void recalculateBoundingBoxes(MeshData& m);

	// Combine a list of meshes to a single mesh container
	MeshFileHeader mergeMeshData(MeshData& m, const std::vector<MeshData*> md);
} // namespace commonFramework