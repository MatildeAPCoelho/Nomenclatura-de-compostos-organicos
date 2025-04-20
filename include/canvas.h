#ifndef CANVAS_H
#define CANVAS_H

#include "line_renderer.h"
#include "molecule.h"
#include "molecule_renderer.h"

#include <GLFW/glfw3.h>
#include <vector>

class Canvas {
public:
    int cursor = GLFW_CROSSHAIR_CURSOR;

    std::vector<Molecule*> molecules;

    glm::vec2 windowSize;
    bool keys[1024] = {0};
    bool heldKeys[1024] = {0};
    bool mouseKeys[8] = {0};
    glm::vec2 mousePos;

    Canvas(glm::vec2 windowSize);
    ~Canvas();

    void init();
    void processInput();
    void render();

private:
    MoleculeRenderer* renderer;
    LineRenderer* lineRenderer;
    ShapeRenderer* shapeRenderer;

    Atom* hoveredAtom = nullptr;
};

#endif