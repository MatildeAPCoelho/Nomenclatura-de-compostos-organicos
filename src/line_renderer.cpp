#include "line_renderer.h"
#include "shader.h"

LineRenderer::LineRenderer(Shader *shader) {
    this->shader = shader;
    this->initRenderData();
}

LineRenderer::~LineRenderer() {
    glDeleteVertexArrays(1, &lineVAO);
}

void LineRenderer::drawLine(glm::vec2 pos1, glm::vec2 pos2, glm::vec4 color, float width) {
    this->shader->Use();
    
    this->shader->SetVector2f("pos1", pos1);
    this->shader->SetVector2f("pos2", pos2);
    this->shader->SetVector4f("color", color);

    glBindVertexArray(lineVAO);
    glEnable(GL_LINE_SMOOTH);
    glLineWidth(width);
    glDrawArrays(GL_LINES, 0, 2);
    glBindVertexArray(0);
    glDisable(GL_LINE_SMOOTH);
}

void LineRenderer::initRenderData() {
    // configure VAO/VBO
    unsigned int VBO;
    float vertices[] = {
        0.0f, 0.0f, 
        1.0f, 1.0f
    };

    glGenVertexArrays(1, &lineVAO);
    glGenBuffers(1, &VBO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    glBindVertexArray(lineVAO);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}