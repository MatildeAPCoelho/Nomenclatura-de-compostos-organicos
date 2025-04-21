#ifndef ATOM_H
#define ATOM_H

#include <vector>
#include <map>
#include <glm/glm.hpp>

struct Chain;

enum species {
    CARBON,
    NITROGEN,
    OXYGEN,
    FLUORINE,
    SULFUR,
    CHLORINE,
    BROMINE,
    IODINE
};

class Atom {
    public:
    std::vector<Atom*> bonds;
    unsigned int species;
    
    glm::vec2 position = glm::vec2(0.0f);

    Atom (unsigned int species, Atom* bond, glm::vec2 position = glm::vec2(0.0f));

    Chain getPrincipalChain (std::vector<Atom*> precedingAtoms, std::map<Atom*, std::vector<unsigned int>> groups, int principalGroup, std::vector<Atom*> ignoredAtoms = {});
};

struct Chain {
    std::vector<unsigned int> principalGroups;
    std::vector<unsigned int> multipleBonds;
    std::vector<unsigned int> doubleBonds;
    std::vector<std::pair<unsigned int, unsigned int>> substituents;
    std::vector<Atom*> chain;
    std::vector<Chain> branches;
    bool isCyclic = false;
    Atom* cyclicGroup = nullptr;
};

std::vector<Chain> compareBranches(std::vector<Chain> chains);
std::vector<Chain> compareChains(std::vector<Chain> chains);
Chain joinChains(Chain one, Chain two);
Chain getOrderedChain(Chain chain);

#endif