#ifndef MOLECULE_RENDERER_H
#define MOLECULE_RENDERER_H

#include "line_renderer.h"
#include "molecule.h"
#include "shape_renderer.h"

class MoleculeRenderer {
public:
    MoleculeRenderer(LineRenderer* lineRenderer, ShapeRenderer* shapeRenderer);

    void drawMolecule(Molecule* molecule); // options
    void drawAtom(Atom* atom);

private:
    LineRenderer* lineRenderer;
    ShapeRenderer* shapeRenderer;
};

#endif