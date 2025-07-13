#ifndef ATOM_H
#define ATOM_H

#include <vector>
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

    Chain getPrincipalChain (std::vector<Atom*> precedingAtoms, std::vector<Atom*> principalGroupAtoms, int principalGroup, std::vector<Atom*> ignoredAtoms = {});
};

struct Chain {
    std::vector<unsigned int> principalGroups;
    std::vector<unsigned int> multipleBonds;
    std::vector<unsigned int> doubleBonds;
    std::vector<unsigned int> substituents;
    std::vector<Atom*> chain;
    std::vector<Chain> branches;
    bool isCyclic = false;
    Atom* cyclicGroup = nullptr;
    bool carbox = false;
    unsigned int halogen = 0;
};

Chain addBranch(Chain chain, Chain branch);
std::vector<Chain> compareBranches(std::vector<Chain> chains);
template <typename T> std::vector<Chain> compareBranchesPropertySize(std::vector<Chain> chains, std::vector<T> (*property)(Chain));
Chain invertChain(Chain chain);
Chain joinChains(Chain one, Chain two);

#endif