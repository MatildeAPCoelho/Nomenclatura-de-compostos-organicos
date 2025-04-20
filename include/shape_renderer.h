#ifndef SHAPE_RENDERER_H
#define SHAPE_RENDERER_H

#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include "shader.h"

class ShapeRenderer {
    public:
        // constructor/destructor
        ShapeRenderer(Shader *shader);
        ~ShapeRenderer();

        void drawRectangle(glm::vec2 position, glm::vec2 size, float rotate = 0.0f, glm::vec4 color = glm::vec4(1.0f));

    private:
        // render state
        Shader *shader;
        unsigned int quadVAO;
        void initRenderData();
};
#endif