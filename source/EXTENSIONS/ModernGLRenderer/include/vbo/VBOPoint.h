#ifndef VBOPOINT_H
#define VBOPOINT_H

#include <vbo/VBObject.h>

/**
 * @brief Stores and manages vertex buffer objects (VBO) for point sprite data.
 **/
class VBOPoint
	: public VBObject
{
	public:
		virtual void render();
};

#endif // VBOPOINT_H
