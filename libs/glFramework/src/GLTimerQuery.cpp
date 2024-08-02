#include <glFramework/include/GLTimerQuery.h>
#include <glFramework/include/kIdxBinds.h>
#include <cstring>

namespace glFramework
{
	GLTimerQuery::GLTimerQuery(const char* beginLabel, const char* endLabel)
	{
		glCreateQueries(GL_TIMESTAMP, 2, queryID_);
		if (nullptr != beginLabel)
		{
			glObjectLabel(GL_QUERY, queryID_[0], std::strlen(beginLabel) + 1, beginLabel);
		}
		if (nullptr != endLabel)
		{
			glObjectLabel(GL_QUERY, queryID_[1], std::strlen(endLabel) + 1, endLabel);
		}
	}

	GLTimerQuery::~GLTimerQuery()
	{
		glDeleteQueries(2, queryID_);
	}

	void GLTimerQuery::beginQuery() const
	{
		glQueryCounter(queryID_[0], GL_TIMESTAMP);
	}

	void GLTimerQuery::endQuery() const
	{
		glQueryCounter(queryID_[1], GL_TIMESTAMP);
	}

	GLuint64 GLTimerQuery::getElapsedTime()
	{
		GLint stopTimerAvailable = 0;
		while (!stopTimerAvailable)
		{
			glGetQueryObjectiv(queryID_[1], GL_QUERY_RESULT_AVAILABLE, &stopTimerAvailable);
		}

		glGetQueryObjectui64v(queryID_[0], GL_QUERY_RESULT, &startTime_);
		glGetQueryObjectui64v(queryID_[1], GL_QUERY_RESULT, &endTime_);

		return (endTime_ - startTime_);
	}
}