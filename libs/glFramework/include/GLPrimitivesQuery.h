#pragma once

#include <glad/glad.h>
#include <glFramework/include/kIdxBinds.h>

namespace glFramework
{
	class GLPrimitivesQuery
	{
	public:
		GLPrimitivesQuery(const char* label = nullptr);
		~GLPrimitivesQuery();

		/** Begin the primitives generated query. Until the query is ended, any primitives that are generated are counted. */
		void beginQuery() const;

		/** Ends primitives generated query query and caches the result - number of primitives generated. */
		void endQuery();

		/** Returns the number of primitives generated.  */
		GLint getNumPrimitivesGenerated() const;

		/** Helper method that returns true if any primitives have been generated.  */
		bool anyPrimitivesGenerated() const;

	private:
		GLuint queryID_{ 0 };
		GLint primitivesGenerated_{ 0 };
	};
}