#pragma once

#include <vector>
#include <unordered_map>
#include <string>

#include <glm/glm.hpp>
#include <glm/ext.hpp>

using std::vector;
using std::unordered_map;
using std::string;

using glm::mat4;

namespace commonFramework
{
	constexpr const int MAX_NODE_LEVEL = 16;

	struct Hierarchy
	{
		// parent for this node (or -1 for root)
		int parent_;
		// first child for a node (or - 1)
		int firstChild_;
		// next sibling for a node (or - 1)
		int nextSibling_;
		// last added node (or - 1)
		int lastSibling_;
		// cached node level
		int level_;
	};

	struct SceneNode
	{
		int mesh_;
		int material_;
		int parent_;
		int firstChild_;
		int rightSibling_;
	};

	struct Scene
	{
		/* Local transformations for each node and global transforms and an array of 'dirty/changed' local transforms */
		vector<mat4> localTransforms_;
		vector<mat4> globalTransforms_;

		// list of nodes whose global transform must be recalculated
		vector<int> changedAtThisFrame_[MAX_NODE_LEVEL];

		// Hierarchy components
		vector<Hierarchy> hierarchy_;

		// Meshes for nodes (Node -> Mesh)
		unordered_map<uint32_t, uint32_t> meshes_;

		// Materials for nodes (Node -> Material) 
		unordered_map<uint32_t, uint32_t> materialForNode_;

		/* Useful for debugging */
		// Node names: which name is assigned to the node
		unordered_map<uint32_t, uint32_t> nameForNode_;

		// Collection of scene node names
		vector<string> names_;

		// Collection of debug material names
		vector<string> materialNames_;
	};

	int addNode(Scene& scene, int parent, int level);

	// markAsChanged() routine starts with a given node and recursively descends to each and every child node, adding it to the changedAtLevel_ arrays. 
	void markAsChanged(Scene& scene, int node);

	int findNodeByName(const Scene& scene, const std::string& name);

	inline string getNodeName(const Scene& scene, int node)
	{
		int strID = scene.nameForNode_.contains(node) ? scene.nameForNode_.at(node) : -1;
		return (strID > -1) ? scene.names_[strID] : std::string();
	}

	inline void setNodeName(Scene& scene, int node, const string& name)
	{
		uint32_t stringID = (uint32_t)scene.names_.size();
		scene.names_.push_back(name);
		scene.nameForNode_[node] = stringID;
	}

	int getNodeLevel(const Scene& scene, int n);

	void recalculateGlobalTransforms(Scene& scene);

	void loadScene(const char* filename, Scene& scene);
	void saveScene(const char* filename, const Scene& scene);

	void dumpTransformations(const char* filename, const Scene& scene);
	void printChagedNodes(const Scene& scene);

	void dumpSceneToDot(const char* filename, const Scene& scene, int* visited = nullptr);

	void mergeScenes(Scene& scene, const vector<Scene*>& scenes, const vector<mat4>& rootTransforms, const vector<uint32_t>& meshCounts, bool mergeMeshes = true, bool mergeMaterials = true);

	// Delete a collection of nodes from a scenegraph
	void deleteSceneNodes(Scene& scene, const vector<uint32_t>& nodesToDelete);
}