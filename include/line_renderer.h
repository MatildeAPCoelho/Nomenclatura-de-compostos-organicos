#ifndef LINE_RENDERER_H
#define LINE_RENDERER_H

#include "shader.h"

class LineRenderer {
public:
    LineRenderer(Shader* shader);
    ~LineRenderer();

    void drawLine(glm::vec2 pos1, glm::vec2 pos2, glm::vec4 color, float width);

private:
    Shader *shader;
    unsigned int lineVAO;
    void initRenderData();
};

#endif