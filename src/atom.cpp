#include "atom.h"
#include <vector>

Atom::Atom (unsigned int species, Atom* bond, glm::vec2 position) {
    this->species = species;
    this->position = position;
    if (bond != nullptr) {
        this->bonds.push_back(bond);
        bond->bonds.push_back(this);
    }
}

template <typename T>
unsigned int getPropertySize(Chain chain, std::vector<T> (*property)(Chain)) {
    int size = property(chain).size();
    if (chain.branches.size() != 0)
        size += getPropertySize(chain.branches[0], property);
    return size;
}

Chain addBranch(Chain chain, Chain branch) {
    chain.principalGroups.insert(chain.principalGroups.begin(), branch.principalGroups.begin(), branch.principalGroups.end());
    chain.multipleBonds.insert(chain.multipleBonds.begin(), branch.multipleBonds.begin(), branch.multipleBonds.end());
    chain.doubleBonds.insert(chain.doubleBonds.begin(), branch.doubleBonds.begin(), branch.doubleBonds.end());
    chain.chain.insert(chain.chain.begin(), branch.chain.begin(), branch.chain.end());
    chain.substituents.insert(chain.substituents.begin(), branch.substituents.begin(), branch.substituents.end());
    chain.branches = branch.branches;
    return chain;
}

template <typename T>
std::vector<Chain> compareBranchesPropertySize(std::vector<Chain> chains, std::vector<T> (*property)(Chain)) {
    for (int i = 0; i < chains.size(); i++) {
        while (chains[i].branches.size() > 0) {
            for (int j = 1; j < chains[i].branches.size(); j++) {
                chains.push_back(addBranch(chains[i], chains[i].branches[j]));
            }
            chains[i] = (addBranch(chains[i], chains[i].branches[0]));
        }
        if (i == 0) {
            continue;
        }

        if (getPropertySize(chains[0], property) > getPropertySize(chains[i], property)) {
            chains.erase(chains.begin() + i);
            i--;
        }
        else if (getPropertySize(chains[0], property) < getPropertySize(chains[i], property)) {
            while (i > 0) {
                chains.erase(chains.begin());
                i--;
            }
        }
    }
    return chains;
}

std::vector<Chain> compareBranches(std::vector<Chain> chains) {
    // greatest number of principal groups
    chains = compareBranchesPropertySize<unsigned int>(chains, [](Chain c) {return c.principalGroups;});
    
    // greatest number of multiple bonds
    chains = compareBranchesPropertySize<unsigned int>(chains, [](Chain c) {return c.multipleBonds;});
    
    // greatest number of atoms
    chains = compareBranchesPropertySize<Atom*>(chains, [](Chain c) {return c.chain;});
    
    // greatest number of double bonds
    chains = compareBranchesPropertySize<unsigned int>(chains, [](Chain c) {return c.doubleBonds;});
    
    return chains;
}

Chain invertChain(Chain chain) {
    Chain invertedChain;
    
    if (chain.isCyclic) {
        int zeros = 0;
        for (int i = (int)chain.principalGroups.size()-1; i >= 0; i--)
            if (chain.principalGroups[i] == 0)
                zeros += 1;
            else
                invertedChain.principalGroups.push_back(chain.chain.size() - chain.principalGroups[i]);
        for (int i = 0; i < zeros; i++) {
            invertedChain.principalGroups.push_back(0);
        }

        zeros = 0;
        for (int i = (int)chain.substituents.size()-1; i >= 0; i--)
            if (chain.substituents[i] == 0)
                zeros += 1;
            else
                invertedChain.substituents.push_back(chain.chain.size() - chain.substituents[i]);
        for (int i = 0; i < zeros; i++) {
            invertedChain.substituents.push_back(0);
        }

        for (int i = (int)chain.multipleBonds.size()-1; i >= 0; i--)
            invertedChain.multipleBonds.push_back(chain.chain.size() - chain.multipleBonds[i] - 1);
        for (int i = (int)chain.doubleBonds.size()-1; i >= 0; i--)
            invertedChain.doubleBonds.push_back(chain.chain.size() - chain.doubleBonds[i] - 1);
        for (int i = (int)chain.chain.size()-2; i >= 0; i--)
            invertedChain.chain.push_back(chain.chain[i]);
        invertedChain.chain.push_back(chain.chain[(int)chain.chain.size()-1]);
    }
    else {
        for (int i = (int)chain.principalGroups.size()-1; i >= 0; i--)
            invertedChain.principalGroups.push_back(chain.chain.size() - chain.principalGroups[i] - 1);
        for (int i = (int)chain.substituents.size()-1; i >= 0; i--)
            invertedChain.substituents.push_back(chain.chain.size() - chain.substituents[i] - 1);
        for (int i = (int)chain.multipleBonds.size()-1; i >= 0; i--)
            invertedChain.multipleBonds.push_back(chain.chain.size() - chain.multipleBonds[i] - 2);
        for (int i = (int)chain.doubleBonds.size()-1; i >= 0; i--)
            invertedChain.doubleBonds.push_back(chain.chain.size() - chain.doubleBonds[i] - 2);
        for (int i = (int)chain.chain.size()-1; i >= 0; i--)
            invertedChain.chain.push_back(chain.chain[i]);
    }

    invertedChain.isCyclic = chain.isCyclic;
    return invertedChain;
}

Chain joinChains(Chain one, Chain two) {
    Chain joinedChain;
    // principalGroups
    for (int i = 0; i < two.principalGroups.size(); i++) {
        joinedChain.principalGroups.push_back(one.chain.size() + two.principalGroups[i] - 1);
    }
    joinedChain.principalGroups.insert(joinedChain.principalGroups.end(), one.principalGroups.begin(), one.principalGroups.end());
    
    // substituents
    for (int i = 0; i < two.substituents.size(); i++) {
        joinedChain.substituents.push_back(one.chain.size() + two.substituents[i] - 1);
    }
    joinedChain.substituents.insert(joinedChain.substituents.end(), one.substituents.begin(), one.substituents.end());
    
    // multipleBonds
    for (int i = 0; i < two.multipleBonds.size(); i++) {
        joinedChain.multipleBonds.push_back(one.chain.size() + two.multipleBonds[i] - 1);
    }
    joinedChain.multipleBonds.insert(joinedChain.multipleBonds.end(), one.multipleBonds.begin(), one.multipleBonds.end());
    
    // doubleBonds
    for (int i = 0; i < two.doubleBonds.size(); i++) {
        joinedChain.doubleBonds.push_back(one.chain.size() + two.doubleBonds[i] - 1);
    }
    joinedChain.doubleBonds.insert(joinedChain.doubleBonds.end(), one.doubleBonds.begin(), one.doubleBonds.end());
    
    // chain
    joinedChain.chain.insert(joinedChain.chain.end(), two.chain.begin(), two.chain.end());
    joinedChain.chain.insert(joinedChain.chain.end(), one.chain.begin(), one.chain.end());
    
    // branches
    joinedChain.branches = two.branches;
    
    return joinedChain;
}

Chain Atom::getPrincipalChain (std::vector<Atom*> precedingAtoms, std::vector<Atom*> principalGroupAtoms, int principalGroup, std::vector<Atom*> ignoredAtoms) {
    std::vector<Chain> possibleChains;

    // get branches
    Chain chain;
    int nonCarbonSubstituents = 0;
    for (int i = 0; i < bonds.size(); i++) {
        if (bonds[i]->species != CARBON) {
            if (i == 0)
                nonCarbonSubstituents += 1;

            for (int j = 0; j < i; j++) {
                if (bonds[i] == bonds[j]) {
                    break;
                }
                if (j == i - 1) {
                    nonCarbonSubstituents += 1;
                }
            }
            continue;
        }
        if (precedingAtoms.size() > 0 && bonds[i] == precedingAtoms[precedingAtoms.size()-1]) { // double bond behind
            continue;
        }

        bool skip = false;
        for (int j = 0; j < ignoredAtoms.size(); j++) {
            if (bonds[i] == ignoredAtoms[j]) {
                skip = true;
                nonCarbonSubstituents += 1;
                break;
            }
        }
        if (skip)
            continue;

        // multiple bond
        skip = false;
        for (int j = 0; j < i; j++) {
            if (bonds[i] == bonds[j]) {
                skip = true;

                for (int k = 0; k < possibleChains.size(); k++) {
                    if (possibleChains[k].chain.size() > 0 && possibleChains[k].chain.back() == bonds[i]) {
                        // triple bond
                        if (possibleChains[k].doubleBonds.size() != 0 && possibleChains[k].doubleBonds.back() == precedingAtoms.size()) {
                            possibleChains[k].doubleBonds.erase(possibleChains[k].doubleBonds.end() - 1);
                        }
                        // double bond
                        else {
                            possibleChains[k].multipleBonds.push_back(precedingAtoms.size());
                            possibleChains[k].doubleBonds.push_back(precedingAtoms.size());
                        }
                    }
                }
                break;
            }
        }
        if (skip)
            continue;

        // add substituents to cyclic chain
        if (possibleChains.size() == 1 && possibleChains[0].isCyclic) {
            if (precedingAtoms.size() > 0 || bonds[i] != possibleChains[0].chain[0]) {
                possibleChains[0].substituents.push_back(precedingAtoms.size());
            }
            continue;
        }
        
        // the principal chain is cyclic
        if (precedingAtoms.size() > 0 && precedingAtoms[0] == bonds[i]) {
            chain = Chain();
            chain.chain.push_back(this);
            chain.isCyclic = true;

            // look for multiple bonds
            for (int j = i+1; j < bonds.size(); j++) {
                if (bonds[i] == bonds[j]) {
                    // triple bond
                    if (chain.doubleBonds.size() != 0 && chain.doubleBonds.back() == precedingAtoms.size()) {
                        chain.doubleBonds.erase(chain.doubleBonds.end() - 1);
                    }
                    // double bond
                    else {
                        chain.multipleBonds.push_back(precedingAtoms.size());
                        chain.doubleBonds.push_back(precedingAtoms.size());
                    }
                }
                else if (bonds[j] != precedingAtoms[precedingAtoms.size()-1]) {
                    chain.substituents.push_back(precedingAtoms.size());
                }
            }

            // principal groups
            for (int j = 0; j < principalGroupAtoms.size(); j++) {
                if (principalGroupAtoms[j] == this) {
                    chain.principalGroups.push_back(precedingAtoms.size());
                }
            }

            // substituents
            for (int j = 0; j < possibleChains.size(); j++) {
                chain.substituents.push_back(precedingAtoms.size());
            }

            return chain;
        }

        // the principal chain has a cyclic substituent
        for (int j = 1; j + 1 < precedingAtoms.size(); j++) {
            if (precedingAtoms[j] == bonds[i]) {
                chain = Chain();
                chain.cyclicGroup = precedingAtoms[j-1];
                return chain;
            }
        }

        // branch
        std::vector<Atom*> currentChain = precedingAtoms;
        currentChain.push_back(this);
        chain = bonds[i]->getPrincipalChain(currentChain, principalGroupAtoms, principalGroup, ignoredAtoms);

        if (chain.isCyclic) {
            for (int j = 0; j < possibleChains.size(); j++) {
                chain.substituents.push_back(precedingAtoms.size());
            }
            possibleChains.clear();
            possibleChains.push_back(chain);
            continue;
        }
        if (chain.cyclicGroup != nullptr) {
            if (this == chain.cyclicGroup) {
                chain.cyclicGroup = nullptr;
                chain.substituents.push_back(precedingAtoms.size());
            }
            else {
                return chain;
            }
        }

        if (chain.chain.size() != 0) {
            possibleChains.push_back(chain);
        }
        else {
            nonCarbonSubstituents += 1;
        }
    }

    if (possibleChains.size() == 0) {
        chain = Chain();
        possibleChains.push_back(chain);
    }

    Chain principalChain;

    std::vector<Chain> branches = compareBranches(possibleChains);
    std::vector<unsigned int> substituents;
    
    if (precedingAtoms.size() == 0) {
        for (int i = 2; i < (int)possibleChains.size(); i++) {
            substituents.push_back(precedingAtoms.size());
        }
    }
    else {
        for (int i = 1; i < possibleChains.size(); i++) {
            substituents.push_back(precedingAtoms.size());
        }
    }

    if (precedingAtoms.size() == 0 && possibleChains.size() > 1) {
        bool skipBranches = false;
        for (int i = 0; i < branches.size(); i++) {
            if (branches[i].chain[branches[i].chain.size()-1] != branches[0].chain[branches[0].chain.size()-1]) {
                skipBranches = true;
                break;
            }
        }

        if (skipBranches) {
            principalChain.branches = branches;
        }
        else {
            for (int i = 0; i < possibleChains.size(); i++) {
                if (possibleChains[i].chain[possibleChains[i].chain.size()-1] == branches[0].chain[branches[0].chain.size()-1]) {
                    possibleChains.erase(possibleChains.begin() + i);
                    i--;
                }
            }
            principalChain.branches = compareBranches(possibleChains);
            principalChain.branches.insert(principalChain.branches.begin(), branches.begin(), branches.end());
        }
    }
    else {
        if (branches.size() == 1) {
            principalChain = branches[0];
        }
        else if (possibleChains.size() == 1) {
            principalChain = possibleChains[0];
        }
        else {
            principalChain.branches = branches;
        }
    }
    
    // add current atom to chain
    principalChain.chain.push_back(this);
    principalChain.substituents.insert(principalChain.substituents.end(), substituents.begin(), substituents.end());

    // add principal groups at current atom
    for (int i = 0; i < principalGroupAtoms.size(); i++) {
        if (principalGroupAtoms[i] == this) {
            principalChain.principalGroups.push_back(precedingAtoms.size());
        }
    }

    for (int i = 0; i < nonCarbonSubstituents; i++) {
        principalChain.substituents.push_back(precedingAtoms.size());
    }

    return principalChain;
}


