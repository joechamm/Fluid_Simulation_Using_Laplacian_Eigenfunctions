#pragma once

#include <glad/glad.h>
#include <glFramework/include/kIdxBinds.h>

namespace glFramework 
{
	class GLOcclusionQuery
	{
	public:
		GLOcclusionQuery(const char* label = nullptr);
		~GLOcclusionQuery();

		/** Begin the occlusion query. Until the query is ended, samples that pass the rendering pipeline are counted. */
		void beginQuery() const;

		/** Ends occlusion query and caches the result - number of samples that passed the rendering pipeline. */
		void endQuery();

		/** Gets number of samples that have passed the rendering pipeline.  */
		GLint getNumSamplesPassed() const;

		/** Helper method that returns is any samples have passed the rendering pipeline.  */
		bool anySamplesPassed() const;

	private:
		GLuint queryID_{ 0 }; // OpenGL query object ID
		GLint samplesPassed_{ 0 };
	};
}


