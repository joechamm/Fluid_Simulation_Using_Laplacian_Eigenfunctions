#include <glFramework/include/GLOcclusionQuery.h>
#include <glFramework/include/kIdxBinds.h>
#include <cstring>

namespace glFramework
{
	GLOcclusionQuery::GLOcclusionQuery(const char* label)
	{
		glCreateQueries(GL_SAMPLES_PASSED, 1, &queryID_);
		if (nullptr != label)
		{
			glObjectLabel(GL_QUERY, queryID_, std::strlen(label) + 1, label);
		}
	}

	GLOcclusionQuery::~GLOcclusionQuery()
	{
		glDeleteQueries(1, &queryID_);
	}

	void GLOcclusionQuery::beginQuery() const
	{
		glBeginQuery(GL_SAMPLES_PASSED, queryID_);
	}

	void GLOcclusionQuery::endQuery()
	{
		glEndQuery(GL_SAMPLES_PASSED);
		glGetQueryObjectiv(queryID_, GL_QUERY_RESULT, &samplesPassed_);
	}

	GLint GLOcclusionQuery::getNumSamplesPassed() const
	{
		return samplesPassed_;
	}

	bool GLOcclusionQuery::anySamplesPassed() const
	{
		return samplesPassed_ > 0;
	}
}