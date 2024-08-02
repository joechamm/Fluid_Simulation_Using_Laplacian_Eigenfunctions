#include <glFramework/include/GLPrimitivesQuery.h>
#include <glFramework/include/kIdxBinds.h>
#include <cstring>

namespace glFramework
{
	GLPrimitivesQuery::GLPrimitivesQuery(const char* label)
	{
		glCreateQueries(GL_PRIMITIVES_GENERATED, 1, &queryID_);
		if (nullptr != label)
		{
			glObjectLabel(GL_QUERY, queryID_, std::strlen(label) + 1, label);
		}
	}

	GLPrimitivesQuery::~GLPrimitivesQuery()
	{
		glDeleteQueries(1, &queryID_);
	}

	void GLPrimitivesQuery::beginQuery() const
	{
		glBeginQuery(GL_PRIMITIVES_GENERATED, queryID_);
	}

	void GLPrimitivesQuery::endQuery()
	{
		glEndQuery(GL_PRIMITIVES_GENERATED);
		glGetQueryObjectiv(queryID_, GL_QUERY_RESULT, &primitivesGenerated_);
	}

	GLint GLPrimitivesQuery::getNumPrimitivesGenerated() const
	{
		return primitivesGenerated_;
	}

	bool GLPrimitivesQuery::anyPrimitivesGenerated() const
	{
		return primitivesGenerated_ > 0;
	}
}