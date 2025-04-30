#ifndef MOLECULE_RENDERER_H
#define MOLECULE_RENDERER_H

#include "line_renderer.h"
#include "molecule.h"
#include "shape_renderer.h"

class MoleculeRenderer {
public:
    MoleculeRenderer(LineRenderer* lineRenderer, ShapeRenderer* shapeRenderer);

    void drawMolecule(Molecule* molecule, Atom* hoveredAtom = nullptr, std::pair<Atom*, Atom*> hoveredBond = std::pair<Atom*, Atom*>{nullptr, nullptr});
    void drawAtom(Atom* atom);

private:
    LineRenderer* lineRenderer;
    ShapeRenderer* shapeRenderer;
};

#endif