#ifndef MOLECULE_H
#define MOLECULE_H

#include "atom.h"

#include <vector>
#include <glm/glm.hpp>
#include <string>

enum substituent {
    CHAIN,
    CARBOXYL,    // OOH
    CARBOALKOXY, // OOR
    HALOFORMYL,  // OXR
    CARBOXAMIDE, // ONRR
    NITRILE,     // =-N
    ALDEHYDE,    // =O
    KETONE,      // =O R
    HYDROXYL,    // OH
    SULFHYDRYL,  // SH
    AMINE        // NRR
};

class Molecule {
public:
    std::vector<Atom*> atoms;
    Chain principalChain;

    Molecule(std::vector<Atom*> atoms);
    ~Molecule();

    std::string getName(bool getPIN = true);
};

#endif