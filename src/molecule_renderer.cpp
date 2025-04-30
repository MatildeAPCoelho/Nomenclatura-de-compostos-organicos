#include "molecule_renderer.h"
#include "atom.h"
#include "shape_renderer.h"
#include <glm/geometric.hpp>

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

void MoleculeRenderer::drawMolecule(Molecule* molecule, Atom* hoveredAtom, std::pair<Atom*, Atom*> hoveredBond) {
    glm::vec3 color;
    bool inChain;
    unsigned int bondOrder;
    float width;
    glm::vec2 bondDir;
    glm::vec2 offset;

    for (int i = 0; i < molecule->atoms.size(); i++) {
        for (int j = 0; j < molecule->atoms[i]->bonds.size(); j++) {
            // already drawn bond
            for (int k = 0; k < i; k++) {
                if (molecule->atoms[i]->bonds[j] == molecule->atoms[k]) {
                    goto loopjEnd;
                }
            }
            for (int k = 0; k < j; k++) {
                if (molecule->atoms[i]->bonds[j] == molecule->atoms[i]->bonds[k]) {
                    goto loopjEnd;
                }
            }

            // get bond order
            bondOrder = 1;
            for (int k = j+1; k < molecule->atoms[i]->bonds.size(); k++) {
                if (molecule->atoms[i]->bonds[j] == molecule->atoms[i]->bonds[k]) {
                    bondOrder++;
                }
            }

            if (molecule->atoms[i]->species == CARBON)
                color = getColor(molecule->atoms[i]->bonds[j]->species);
            else
                color = getColor(molecule->atoms[i]->species);

            inChain = false;
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

            if (hoveredBond.first == molecule->atoms[i] && hoveredBond.second == molecule->atoms[i]->bonds[j] || 
                hoveredBond.first == molecule->atoms[i]->bonds[j] && hoveredBond.second == molecule->atoms[i]) {
                width = 5;
            }
            else {
                width = 3;
            }

            
            lineRenderer->drawLine(molecule->atoms[i]->position, molecule->atoms[i]->bonds[j]->position, glm::vec4(color, 1.0f), width);
            if (bondOrder > 1) {
                bondDir = molecule->atoms[i]->position - molecule->atoms[i]->bonds[j]->position;
                if (bondDir.y == 0) {
                    offset = {0.0f, 1.0f};
                }
                else {
                    offset = {1, -bondDir.x/bondDir.y};
                    offset = glm::normalize(offset);
                }
                for (int k = 1; 2*k < bondOrder; k ++) {
                    lineRenderer->drawLine(molecule->atoms[i]->position + offset*10.0f*(float)k, molecule->atoms[i]->bonds[j]->position + offset*10.0f*(float)k, glm::vec4(color, 1.0f), width);
                    lineRenderer->drawLine(molecule->atoms[i]->position - offset*10.0f*(float)k, molecule->atoms[i]->bonds[j]->position - offset*10.0f*(float)k, glm::vec4(color, 1.0f), width);
                }
                if (bondOrder % 2 == 0) {
                    lineRenderer->drawLine(molecule->atoms[i]->position + offset*10.0f*(float)bondOrder/2.0f, molecule->atoms[i]->bonds[j]->position + offset*10.0f*(float)bondOrder/2.0f, glm::vec4(color, 1.0f), width);
                }
            }

            loopjEnd:;
        }
    }

    if (hoveredAtom != nullptr) {
        drawAtom(hoveredAtom);
    }
}

void MoleculeRenderer::drawAtom(Atom* atom) {
    glm::vec2 size = glm::vec2(10.0f);
    glm::vec2 position = atom->position - size/2.0f;
    glm::vec3 color = getColor(atom->species);

    shapeRenderer->drawRectangle(position, size, 0.0f, glm::vec4(color, 1.0f));
}