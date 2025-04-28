#include "canvas.h"
#include "line_renderer.h"
#include "molecule.h"
#include "resource_manager.h"
#include "shape_renderer.h"

#include <GLFW/glfw3.h>
#include <iostream>

// utility functions
Molecule* getAtomMolecule(std::vector<Molecule*> molecules, Atom* atom);
std::vector<Atom*> getConnectedAtoms(Atom* atom, std::vector<Atom*> viewedAtoms);

Canvas::Canvas(glm::vec2 windowSize) : windowSize(windowSize) {
}

Canvas::~Canvas() {
    delete renderer;
    for (int i = 0; i < molecules.size(); i++) {
        for (int j = 0; j < molecules[i]->atoms.size(); j++) {
            delete molecules[i]->atoms[j];
        }
        delete molecules[i];
    }
}

void Canvas::init() {
    ResourceManager::LoadShader("line", "shaders/line.vs", "shaders/color.frag", nullptr);
    ResourceManager::LoadShader("shape", "shaders/shape.vs", "shaders/color.frag", nullptr);
    lineRenderer = new LineRenderer(ResourceManager::GetShader("line"));
    shapeRenderer = new ShapeRenderer(ResourceManager::GetShader("shape"));
    renderer = new MoleculeRenderer(lineRenderer, shapeRenderer);
}

void Canvas::processInput() {
    // drawing molecules
    static bool movingAtomsMode = false;
    static bool drawingLine = false;
    static unsigned int drawingSpecies = CARBON;
    static Atom* lineOrigin;
    static Atom* lineEnd;

    if (mouseKeys[0]) {
        if (!drawingLine && (!movingAtomsMode || hoveredAtom != nullptr)) {
            drawingLine = true;

            if (hoveredAtom == nullptr) {
                lineOrigin = new Atom(drawingSpecies, nullptr, mousePos);
                molecules.push_back(new Molecule(std::vector<Atom*>{lineOrigin}));
            }
            else {
                lineOrigin = hoveredAtom;
            }

            if (movingAtomsMode)
                lineEnd = hoveredAtom;
            else
                lineEnd = new Atom(drawingSpecies, lineOrigin, mousePos);
            
            getAtomMolecule(molecules, lineOrigin)->atoms.push_back(lineEnd);
        }
        if (drawingLine) {
            lineEnd->position = mousePos;
        }
    }
    else {
        if (drawingLine) {
            if (hoveredAtom != lineEnd && !movingAtomsMode) {
                Molecule* molecule = getAtomMolecule(molecules, lineEnd);
                
                // delete new atom
                for (int i = 0; i < molecule->atoms.size(); i++) {
                    if (molecule->atoms[i] == lineEnd) {
                        molecule->atoms.erase(molecule->atoms.begin() + i);
                    }
                }
                for (int i = 0; i < lineOrigin->bonds.size(); i++) {
                    if (lineOrigin->bonds[i] == lineEnd) {
                        lineOrigin->bonds.erase(lineOrigin->bonds.begin() + i);
                    }
                }

                // create bond
                if (hoveredAtom != lineOrigin) {
                    lineOrigin->bonds.push_back(hoveredAtom);
                    hoveredAtom->bonds.push_back(lineOrigin);

                    // join molecules
                    Molecule* joinedMolecule = getAtomMolecule(molecules, hoveredAtom);
                    if (molecule != joinedMolecule) {
                        molecule->atoms.insert(molecule->atoms.end(), joinedMolecule->atoms.begin(), joinedMolecule->atoms.end());
                        for (int i = 0; i < molecules.size(); i++) {
                            if (molecules[i] == joinedMolecule) {
                                molecules.erase(molecules.begin() + i);
                                delete joinedMolecule;
                                break;
                            }
                        }
                    }
                }

                delete lineEnd;
            }

            drawingLine = false;
        }
    }
    // change drawing atom species
    if (keys[GLFW_KEY_C]) {
        drawingSpecies = CARBON;
        if (hoveredAtom != nullptr)
            hoveredAtom->species = CARBON;
    }
    else if (keys[GLFW_KEY_N]) {
        drawingSpecies = NITROGEN;
        if (hoveredAtom != nullptr)
            hoveredAtom->species = NITROGEN;
    }
    else if (keys[GLFW_KEY_O]) {
        drawingSpecies = OXYGEN;
        if (hoveredAtom != nullptr)
            hoveredAtom->species = OXYGEN;
    }
    else if (keys[GLFW_KEY_F]) {
        drawingSpecies = FLUORINE;
        if (hoveredAtom != nullptr)
            hoveredAtom->species = FLUORINE;
    }
    else if (keys[GLFW_KEY_B]) {
        drawingSpecies = BROMINE;
        if (hoveredAtom != nullptr)
            hoveredAtom->species = BROMINE;
    }
    else if (keys[GLFW_KEY_S]) {
        drawingSpecies = SULFUR;
        if (hoveredAtom != nullptr)
            hoveredAtom->species = SULFUR;
    }
    // delete atom
    else if (keys[GLFW_KEY_DELETE]) {
        if (hoveredAtom != nullptr) {
            Molecule* molecule = getAtomMolecule(molecules, hoveredAtom);

            // delete atom from molecule
            for (int i = 0; i < molecule->atoms.size(); i++) {
                if (molecule->atoms[i] == hoveredAtom) {
                    molecule->atoms.erase(molecule->atoms.begin() + i);
                    break;
                }
            }

            // delete atom from bonds
            for (int i = 0; i < hoveredAtom->bonds.size(); i++) {
                for (int j = 0; j < hoveredAtom->bonds[i]->bonds.size(); j++) {
                    if (hoveredAtom->bonds[i]->bonds[j] == hoveredAtom) {
                        hoveredAtom->bonds[i]->bonds.erase(hoveredAtom->bonds[i]->bonds.begin() + j);
                        break;
                    }
                }
            }

            // separate molecules
            for (int i = 0; i < hoveredAtom->bonds.size(); i++) {
                std::vector<Atom*> connectedAtoms = getConnectedAtoms(hoveredAtom->bonds[i], std::vector<Atom *>{hoveredAtom->bonds[i]});
                
                for (int j = i+1; j < hoveredAtom->bonds.size(); j++) {
                    for (int k = 0; k < connectedAtoms.size(); k++) {
                        // connected
                        if (hoveredAtom->bonds[j] == connectedAtoms[k]) {
                            hoveredAtom->bonds.erase(hoveredAtom->bonds.begin() + j);
                            j--;
                            break;
                        }
                    }
                }
                // create new molecule
                molecules.push_back(new Molecule(connectedAtoms));
            }
            
            // delete original molecule
            for (int i = 0; i < molecules.size(); i++) {
                if (molecules[i] == molecule) {
                    molecules.erase(molecules.begin() + i);
                }
            }

            delete hoveredAtom;
            hoveredAtom = nullptr;
        }
    }


    // commands
    // get IUPAC name (QOBI version)
    if (keys[GLFW_KEY_ENTER]) {
        if (!heldKeys[GLFW_KEY_ENTER]) {
            for (int i = 0; i < molecules.size(); i++) {
                std::cout << molecules[i]->getName() << std::endl;;
            }
        }
        heldKeys[GLFW_KEY_ENTER] = true;
    }
    else if (keys[GLFW_KEY_SPACE]) {
        if (!heldKeys[GLFW_KEY_SPACE]) {
            movingAtomsMode = !movingAtomsMode;
            if (movingAtomsMode)
                cursor = GLFW_HAND_CURSOR;
            else
                cursor = GLFW_CROSSHAIR_CURSOR;
        }
        heldKeys[GLFW_KEY_SPACE] = true;
    }
    else {
        heldKeys[GLFW_KEY_ENTER] = false;
        heldKeys[GLFW_KEY_SPACE] = false;
    }


    hoveredAtom = nullptr;
    for (int i = 0; i < molecules.size(); i++) {
        for (int j = 0; j < molecules[i]->atoms.size(); j++) {
            Atom* atom = molecules[i]->atoms[j];
            if (std::abs(mousePos.x - atom->position.x) <= 10 && std::abs(mousePos.y - atom->position.y) <= 10) {
                hoveredAtom = atom;
                break;
            }
        }
        if (hoveredAtom != nullptr && (!drawingLine || hoveredAtom != lineEnd))
            break;
    }
}

void Canvas::render() {
    ResourceManager::GetShader("line")->SetVector2f("windowSize", windowSize, true);
    ResourceManager::GetShader("shape")->SetVector2f("windowSize", windowSize, true);

    for (int i = 0; i < molecules.size(); i++) {
        renderer->drawMolecule(molecules[i]);
    }
    if (hoveredAtom != nullptr) {
        renderer->drawAtom(hoveredAtom);
    }
}

// utility functions
Molecule* getAtomMolecule(std::vector<Molecule*> molecules, Atom* atom) {
    for (int i = 0; i < molecules.size(); i++) {
        for (int j = 0; j < molecules[i]->atoms.size(); j++) {
            if (molecules[i]->atoms[j] == atom) {
                return molecules[i];
            }
        }
    }
    return nullptr;
}

std::vector<Atom*> getConnectedAtoms(Atom* atom, std::vector<Atom*> viewedAtoms) {
    for (int i = 0; i < atom->bonds.size(); i++) {
        for (int j = 0; j < viewedAtoms.size(); j++) {
            if (atom->bonds[i] == viewedAtoms[j]) {
                goto loopiEnd;
            }
        }
        viewedAtoms.push_back(atom->bonds[i]);
        viewedAtoms = getConnectedAtoms(atom->bonds[i], viewedAtoms);
        loopiEnd:;
    }
    return viewedAtoms;
}