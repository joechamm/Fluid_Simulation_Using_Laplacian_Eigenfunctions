#pragma once

#include <helpers/RootDir.h>
#include <string>
#include <vector>

using std::string;

namespace commonFramework {

	inline string appendToRoot(const char* str)
	{
		return string(ROOT_DIR).append(str);
	}

	inline std::vector<string> getShaderFilenamesWithRoot(const char* vtxFilename, const char* fragFilename)
	{
		string vtxFilenameWRoot = appendToRoot(vtxFilename);
		string fragFilenameWRoot = appendToRoot(fragFilename);

		std::vector<string> shaderFilenames = {
			vtxFilenameWRoot,
			fragFilenameWRoot
		};

		return shaderFilenames;
	}

	inline std::vector<string> getShaderFilenamesWithRoot(const char* vtxFilename, const char*, const char* fragFilename, const char* geomFilename)
	{
		string vtxFilenameWRoot = appendToRoot(vtxFilename);
		string fragFilenameWRoot = appendToRoot(fragFilename);
		string geomFilenameWRoot = appendToRoot(geomFilename);

		std::vector<string> shaderFilenames = {
			vtxFilenameWRoot,
			fragFilenameWRoot,
			geomFilenameWRoot
		};

		return shaderFilenames;
	}
}