/*
 * Copyright (C) 2013-2015 Sergey Kosarevsky (sk@linderdaum.com)
 * Copyright (C) 2013-2015 Viktor Latypov (vl@linderdaum.com)
 * Based on Linderdaum Engine http://www.linderdaum.com
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must display the names 'Sergey Kosarevsky' and
 *    'Viktor Latypov' in the credits of the application, if such credits exist.
 *    The authors of this work must be notified via email (sk@linderdaum.com) in
 *    this case of redistribution.
 *
 * 3. Neither the name of copyright holders nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ``AS
 * IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL COPYRIGHT HOLDERS OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef __TRACKBALL_H__
#define __TRACKBALL_H__

#include <numeric>

#include <glm/glm.hpp>
#include <glm/ext.hpp>

using glm::mat4;
using glm::vec4;
using glm::vec3;
using glm::vec2;

/// Virtual trackball for user interaction with rotations
class VirtualTrackball
{
	vec3 pointCur_ = vec3(0.0f);
	vec3 pointPrev_ = vec3(0.0f);
	mat4 rotation_ = mat4(1.0f);
	mat4 rotationDelta_ = mat4(1.0f);
	bool isDraggingActive_ = false;

public:
	VirtualTrackball() = default;

	/* Get rotation matrix for new mouse point */
	mat4 dragTo(vec2 screenPoint, float speed, bool keyPressed)
	{
		if (keyPressed && !isDraggingActive_)
		{
			startDragging(screenPoint);

			isDraggingActive_ = keyPressed;

			return mat4(1.0f);
		}

		isDraggingActive_ = keyPressed;

		if (!keyPressed) return mat4(1.0f);

		pointCur_ = projectOnSphere(screenPoint);

		const vec3 direction = pointCur_ - pointPrev_;
		const float shift = glm::length(direction);

		mat4 rotMatrix = mat4(1.0f);

		if (shift > std::numeric_limits<float>::epsilon())
		{
			const vec3 axis = glm::cross(pointPrev_, pointCur_);
			rotMatrix = glm::rotate(mat4(1.0f), shift * speed, axis);
		}

		rotationDelta_ = rotMatrix;

		return rotMatrix;
	}

	const mat4& getRotationDelta() const
	{
		return rotationDelta_;
	}

	mat4 getRotationMatrix() const
	{
		return rotation_ * rotationDelta_;
	}

private:
	void startDragging(vec2 screenPoint)
	{
		rotation_ = rotation_ * rotationDelta_;
		rotationDelta_ = mat4(1.0f);
		pointCur_ = projectOnSphere(screenPoint);
		pointPrev_ = pointCur_;
	}

	vec3 projectOnSphere(vec2 screenPoint)
	{
		// convert to -1.0...1.0 range
		vec3 proj(
			+(2.0f * screenPoint.x - 1.0f),
			-(2.0f * screenPoint.y - 1.0f),
			0.0f
		);

		const float len = std::min(glm::length(proj), 1.0f);

		proj.z = sqrtf(1.001f - len * len);

		return glm::normalize(proj);
	}
};
#endif // !__TRACKBALL_H__