#include "shape_renderer.h"

ShapeRenderer::ShapeRenderer(Shader *shader) {
    this->shader = shader;
    this->initRenderData();
}

ShapeRenderer::~ShapeRenderer() {
    glDeleteVertexArrays(1, &quadVAO);
}

void ShapeRenderer::drawRectangle(glm::vec2 position, glm::vec2 size, float rotate, glm::vec4 color) {
    shader->Use();

    // prepare transformations
    glm::mat4 model = glm::mat4(1.0f);
    model = glm::translate(model, glm::vec3(position, 0.0f));
    model = glm::scale(model, glm::vec3(size, 1.0f));
    model = glm::translate(model, glm::vec3(0.5f, 0.5f, 0.0f)); // move origin back
    model = glm::rotate(model, glm::radians(rotate), glm::vec3(0.0f, 0.0f, 1.0f));
    model = glm::translate(model, glm::vec3(-0.5f, -0.5f, 0.0f)); // move origin of rotation to center of quad

    shader->SetMatrix4("model", model);

    // render texture quad
    shader->SetVector4f("color", color);
    
    glBindVertexArray(this->quadVAO);
    glDrawArrays(GL_TRIANGLES, 0, 6);
    glBindVertexArray(0);
}

void ShapeRenderer::initRenderData() {
    // configure VAO/VBO
    unsigned int VBO;
    float vertices[] = {
        0.0f, 1.0f,
        1.0f, 0.0f,
        0.0f, 0.0f,

        0.0f, 1.0f,
        1.0f, 1.0f,
        1.0f, 0.0f
    };

    glGenVertexArrays(1, &quadVAO);
    glGenBuffers(1, &VBO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    glBindVertexArray(quadVAO);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}