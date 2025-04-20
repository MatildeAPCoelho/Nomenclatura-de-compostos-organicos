#include "canvas.h"
#include "line_renderer.h"
#include "molecule.h"
#include "resource_manager.h"
#include "shape_renderer.h"

#include <GLFW/glfw3.h>
#include <iostream>

// utility functions
Molecule* getAtomMolecule(std::vector<Molecule*> molecules, Atom* atom);

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
    

    std::vector<Atom*> atoms;
    // exercício 8.a)
    atoms.push_back(new Atom(CARBON, nullptr, glm::vec2(0.1, 0.5)*windowSize));
    atoms.push_back(new Atom(CARBON, atoms[0], glm::vec2(0.2, 0.6)*windowSize));
    atoms.push_back(new Atom(OXYGEN, atoms[1], glm::vec2(0.3, 0.5)*windowSize));
    atoms.push_back(new Atom(CARBON, atoms[2], glm::vec2(0.4, 0.6)*windowSize)); // 3
    atoms.push_back(new Atom(CARBON, atoms[3], glm::vec2(0.4, 0.7)*windowSize)); // 4
    atoms.push_back(new Atom(CARBON, atoms[4], glm::vec2(0.3, 0.8)*windowSize));
    atoms.push_back(new Atom(OXYGEN, atoms[5], glm::vec2(0.2, 0.7)*windowSize));
    atoms[5]->bonds.push_back(atoms[6]);
    atoms[6]->bonds.push_back(atoms[5]);
    atoms.push_back(new Atom(CARBON, atoms[4], glm::vec2(0.5, 0.8)*windowSize));
    atoms.push_back(new Atom(CARBON, atoms[7], glm::vec2(0.5, 0.9)*windowSize));
    atoms.push_back(new Atom(CARBON, atoms[3], glm::vec2(0.5, 0.5)*windowSize)); // 9
    atoms.push_back(new Atom(CARBON, atoms[9], glm::vec2(0.6, 0.6)*windowSize)); // 10
    atoms.push_back(new Atom(BROMINE, atoms[10], glm::vec2(0.6, 0.7)*windowSize));
    atoms.push_back(new Atom(CARBON, atoms[10], glm::vec2(0.7, 0.5)*windowSize)); // 12
    atoms.push_back(new Atom(NITROGEN, atoms[12], glm::vec2(0.7, 0.4)*windowSize));
    atoms.push_back(new Atom(CARBON, atoms[12], glm::vec2(0.8, 0.6)*windowSize)); // 14
    atoms.push_back(new Atom(NITROGEN, atoms[14], glm::vec2(0.9, 0.5)*windowSize));
    atoms.push_back(new Atom(OXYGEN, atoms[14], glm::vec2(0.8, 0.7)*windowSize));
    atoms[14]->bonds.push_back(atoms[16]);
    atoms[16]->bonds.push_back(atoms[14]);
    atoms.push_back(new Atom(CARBON, atoms[9], glm::vec2(0.5, 0.4)*windowSize)); // 17
    atoms.push_back(new Atom(CARBON, atoms[17], glm::vec2(0.6, 0.3)*windowSize));
    atoms[17]->bonds.push_back(atoms[18]);
    atoms[18]->bonds.push_back(atoms[17]);
    atoms.push_back(new Atom(CARBON, atoms[18], glm::vec2(0.6, 0.2)*windowSize)); // 19
    atoms.push_back(new Atom(SULFUR, atoms[19], glm::vec2(0.7, 0.1)*windowSize));
    atoms.push_back(new Atom(CARBON, atoms[19], glm::vec2(0.5, 0.1)*windowSize));
    atoms[19]->bonds.push_back(atoms[21]);
    atoms[21]->bonds.push_back(atoms[19]);
    atoms.push_back(new Atom(CARBON, atoms[21], glm::vec2(0.4, 0.2)*windowSize));
    atoms.push_back(new Atom(CARBON, atoms[22], glm::vec2(0.4, 0.3)*windowSize));
    atoms[22]->bonds.push_back(atoms[23]);
    atoms[23]->bonds.push_back(atoms[22]);
    atoms.push_back(new Atom(OXYGEN, atoms[23], glm::vec2(0.3, 0.4)*windowSize));
    atoms[17]->bonds.push_back(atoms[23]);
    atoms[23]->bonds.push_back(atoms[17]);
    
    

    //molecules.push_back(new Molecule(atoms));
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
                }

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

                delete lineEnd;
            }

            drawingLine = false;
        }
    }
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


    // commands
    // get IUPAC name (QOBI version)
    if (keys[GLFW_KEY_ENTER]) {
        if (!heldKeys[GLFW_KEY_ENTER]) {
            for (int i = 0; i < molecules.size(); i++) {
                molecules[i]->getName();
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