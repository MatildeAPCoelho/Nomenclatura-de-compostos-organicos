#include "molecule_renderer.h"
#include "atom.h"
#include "shape_renderer.h"

MoleculeRenderer::MoleculeRenderer(LineRenderer* lineRenderer, ShapeRenderer* shapeRenderer) {
    this->lineRenderer = lineRenderer;
    this->shapeRenderer = shapeRenderer;
}

glm::vec3 getColor(unsigned int species) {
    switch (species) {
        case OXYGEN: // red
            return glm::vec3(1.0f, 0.0f, 0.0f);

        case NITROGEN: // blue
            return glm::vec3(0.0f, 0.0f, 1.0f);
            break;

        case SULFUR: // yellow
            return glm::vec3(1.0f, 1.0f, 0.0f);

        case FLUORINE: // green
        case CHLORINE:
        case BROMINE:
        case IODINE:
            return glm::vec3(0.0f, 1.0f, 0.0f);
        default: // white
          return glm::vec3(1.0f);
    }
}

void MoleculeRenderer::drawMolecule(Molecule* molecule) {
    bool alreadyDrawn = false;
    for (int i = 0; i < molecule->atoms.size(); i++) {
        for (int j = 0; j < molecule->atoms[i]->bonds.size(); j++) {
            alreadyDrawn = false;
            for (int k = 0; k < i; k++) {
                if (molecule->atoms[i]->bonds[j] == molecule->atoms[j]) {
                    alreadyDrawn = true;
                    break;
                }
            }
            if (alreadyDrawn)
                continue;

            glm::vec3 color;
            if (molecule->atoms[i]->species == CARBON)
                color = getColor(molecule->atoms[i]->bonds[j]->species);
            else
                color = getColor(molecule->atoms[i]->species);

            bool inChain = false;
            for (int l = 0; l < molecule->principalChain.chain.size(); l++) {
                if (molecule->principalChain.chain[l] == molecule->atoms[i] || molecule->principalChain.chain[l] == molecule->atoms[i]->bonds[j]) {
                    if (inChain) {
                        color = glm::vec3(1.0f, 0.0f, 1.0f);
                    }
                    else {
                        inChain = true;
                    }
                }
            }

            lineRenderer->drawLine(molecule->atoms[i]->position, molecule->atoms[i]->bonds[j]->position, glm::vec4(color, 1.0f));
        }
    }
}

void MoleculeRenderer::drawAtom(Atom* atom) {
    glm::vec2 size = glm::vec2(10.0f);
    glm::vec2 position = atom->position - size/2.0f;
    glm::vec3 color = getColor(atom->species);

    shapeRenderer->drawRectangle(position, size, 0.0f, glm::vec4(color, 1.0f));
}