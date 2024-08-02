#pragma once

#include <glad/glad.h>
#include <glFramework/include/kIdxBinds.h>

namespace glFramework
{
	class GLTimerQuery
	{
	public:
		GLTimerQuery(const char* beginLabel = nullptr, const char* endLabel = nullptr);
		~GLTimerQuery();

		void beginQuery() const;
		void endQuery() const;

		GLuint64 getElapsedTime();

	private:
		GLuint queryID_[2] = { 0, 0 };
		GLuint64 startTime_{ 0 };
		GLuint64 endTime_{ 0 };
	};
}