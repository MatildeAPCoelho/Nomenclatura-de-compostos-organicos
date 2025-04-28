#include "atom.h"
#include <vector>
#include <iostream>
Atom::Atom (unsigned int species, Atom* bond, glm::vec2 position) {
    this->species = species;
    this->position = position;
    if (bond != nullptr) {
        this->bonds.push_back(bond);
        bond->bonds.push_back(this);
    }
}

template <typename T>
unsigned int getPropertySize(Chain chain, std::vector<T> (*comparator)(Chain)) {
    int size = comparator(chain).size();
    if (chain.branches.size() != 0)
        size += getPropertySize(chain.branches[0], comparator);
    return size;
}

template <typename T>
std::vector<Chain> compareBranchesPropertySize(std::vector<Chain> chains, std::vector<T> (*comparator)(Chain)) {
    for (int i = 1; i < chains.size(); i++) {
        if (getPropertySize(chains[0], comparator) > getPropertySize(chains[i], comparator)) {
            chains.erase(chains.begin() + i);
            i--;
        }
        else if (getPropertySize(chains[0], comparator) < getPropertySize(chains[i], comparator)) {
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

Chain addBranch(Chain chain, Chain branch) {
    chain.principalGroups.insert(chain.principalGroups.end(), branch.principalGroups.begin(), branch.principalGroups.end());
    chain.multipleBonds.insert(chain.multipleBonds.end(), branch.multipleBonds.begin(), branch.multipleBonds.end());
    chain.doubleBonds.insert(chain.doubleBonds.end(), branch.doubleBonds.begin(), branch.doubleBonds.end());
    chain.chain.insert(chain.chain.begin(), branch.chain.begin(), branch.chain.end());
    chain.branches = branch.branches;
    return chain;
}

Chain invertChain(Chain chain) {
    Chain invertedChain;
    for (int i = 0; i < chain.principalGroups.size(); i++)
        invertedChain.principalGroups.push_back(chain.chain.size() - chain.principalGroups[i] - 1);
    for (int i = 0; i < chain.multipleBonds.size(); i++)
        invertedChain.multipleBonds.push_back(chain.chain.size() - chain.multipleBonds[i] - 2);
    for (int i = 0; i < chain.doubleBonds.size(); i++)
        invertedChain.doubleBonds.push_back(chain.chain.size() - chain.doubleBonds[i] - 2);
    for (int i = 0; i < chain.chain.size(); i++)
        invertedChain.chain.push_back(chain.chain[chain.chain.size()-i-1]);

    return invertedChain;
}

Chain joinChains(Chain one, Chain two) {
    Chain joinedChain;
    // principalGroups
    for (int i = 0; i < one.principalGroups.size(); i++) {
        joinedChain.principalGroups.push_back(one.chain.size() - one.principalGroups[i] - 1);
    }
    for (int i = 0; i < two.principalGroups.size(); i++) {
        joinedChain.principalGroups.push_back(one.chain.size() + two.principalGroups[i] - 1);
    }
    // multipleBonds
    for (int i = 0; i < one.multipleBonds.size(); i++) {
        joinedChain.multipleBonds.push_back(one.chain.size() - one.multipleBonds[i] - 2);
    }
    for (int i = 0; i < two.multipleBonds.size(); i++) {
        joinedChain.multipleBonds.push_back(one.chain.size() + two.multipleBonds[i] - 1);
    }
    // doubleBonds
    for (int i = 0; i < one.doubleBonds.size(); i++) {
        joinedChain.doubleBonds.push_back(one.chain.size() - one.doubleBonds[i] - 2);
    }
    for (int i = 0; i < two.doubleBonds.size(); i++) {
        joinedChain.doubleBonds.push_back(one.chain.size() + two.doubleBonds[i] - 1);
    }
    // chain
    for (int i = 0; i < one.chain.size(); i++) {
        joinedChain.chain.push_back(one.chain[i]);
    }
    for (int i = 0; i < two.chain.size(); i++) {
        joinedChain.chain.push_back(two.chain[two.chain.size() - i - 1]);
    }
    
    return joinedChain;
}

template <typename T>
std::vector<Chain> compareBranchesPropertyPosition(std::vector<Chain> chains, std::vector<T> (*property)(Chain)) {
    for (int i = 1; i < chains.size(); i++) {
        int j = 0;
        while (true) {
            if (j == property(chains[0]).size()) {
                if (chains[0].branches.size() == 0)
                    break;
                
                chains[0] = addBranch(chains[0], chains[0].branches[0]);
                for (int k = 1; k < chains[0].branches.size(); k++)
                    chains.push_back(addBranch(chains[0], chains[0].branches[k]));
                continue;
            }
            if (j == property(chains[i]).size()) {
                if (chains[i].branches.size() == 0)
                    break;

                chains[i] = addBranch(chains[i], chains[i].branches[0]);
                for (int k = 1; k < chains[i].branches.size(); k++)
                    chains.push_back(addBranch(chains[i], chains[i].branches[k]));
                continue;
            }

            if (property(chains[0])[j] < property(chains[i])[j]) {
                chains.erase(chains.begin() + i);
                i--;
            }
            else if (property(chains[0])[j] > property(chains[i])[j]) {
                while (i > 0) {
                    chains.erase(chains.begin());
                    i--;
                }
            }

            j++;
        }
    }
    return chains;
}

template <typename T>
std::vector<Chain> compareBranchesPropertyPositionBackward(std::vector<Chain> chains, std::vector<T> (*property)(Chain)) {
    for (int i = 1; i < chains.size(); i++) {
        while (chains[0].branches.size() > 0) {
            chains[0] = addBranch(chains[0], chains[0].branches[0]);
            for (int k = 1; k < chains[0].branches.size(); k++)
                chains.push_back(addBranch(chains[0], chains[0].branches[k]));
        }
        while (chains[i].branches.size() > 0) {
            chains[i] = addBranch(chains[i], chains[i].branches[0]);
            for (int k = 1; k < chains[i].branches.size(); k++)
                chains.push_back(addBranch(chains[i], chains[i].branches[k]));
        }

        for (int j = 0; j < property(chains[0]).size(); j++) {
            if (chains[0].chain.size() - property(chains[0])[property(chains[0]).size()-j-1] < chains[i].chain.size() - property(chains[i])[property(chains[i]).size()-j-1]) {
                chains.erase(chains.begin() + i);
                i--;
                break;
            }
            else if (chains[0].chain.size() - property(chains[0])[property(chains[0]).size()-j-1] > chains[i].chain.size() - property(chains[i])[property(chains[i]).size()-j-1]) {
                while (i > 0) {
                    chains.erase(chains.begin());
                    i--;
                }
                break;
            }
        }
    }
    return chains;
}

std::vector<Chain> compareChains(std::vector<Chain> chains) {
    chains = compareBranches(chains);

    // lower numbers for suffixes
    chains = compareBranchesPropertyPosition<unsigned int>(chains, [](Chain c) {return c.principalGroups;});

    // lower numbers for multiple bonds
    chains = compareBranchesPropertyPosition<unsigned int>(chains, [](Chain c) {return c.multipleBonds;});

    // lower numbers for double bonds
    chains = compareBranchesPropertyPosition<unsigned int>(chains, [](Chain c) {return c.doubleBonds;});

    // Maior nº de prefixos (substituintes)
    // Nºs mais baixos para os prefixos da cadeia fundamental
    // Nºs mais baixos para os 1ºs prefixos
    return chains;
}

Chain getOrderedChain(Chain chain) {
    std::vector<Chain> bestChains;

    if (chain.chain.size() > 1) {
        std::vector<Chain> forwardChains;
        std::vector<Chain> backwardChains;
        if (chain.branches.size() == 0) {
            forwardChains = std::vector<Chain>{chain};
            backwardChains = std::vector<Chain>{chain};
        }
        else {
            for (int i = 0; i < chain.branches.size(); i++) {
                forwardChains.push_back(addBranch(chain, chain.branches[i]));
            }
            // lower numbers for principal groups
            forwardChains = compareBranchesPropertyPosition<unsigned int>(forwardChains, [](Chain c) {return c.principalGroups;});
            // lower numbers for multiple bonds
            forwardChains = compareBranchesPropertyPosition<unsigned int>(forwardChains, [](Chain c) {return c.multipleBonds;});
            // lower numbers for double bonds
            forwardChains = compareBranchesPropertyPosition<unsigned int>(forwardChains, [](Chain c) {return c.doubleBonds;});

            while (forwardChains[0].branches.size() > 0){
                forwardChains[0] = addBranch(forwardChains[0], forwardChains[0].branches[0]);
            }

            // Maior nº de prefixos (substituintes)
            // Nºs mais baixos para os prefixos da cadeia fundamental
            // Nºs mais baixos para os 1ºs prefixos
            for (int i = 0; i < chain.branches.size(); i++) {
                backwardChains.push_back(addBranch(chain, chain.branches[i]));
            }
            backwardChains = compareBranchesPropertyPositionBackward<unsigned int>(backwardChains, [](Chain c) {return c.principalGroups;});
            backwardChains = compareBranchesPropertyPositionBackward<unsigned int>(backwardChains, [](Chain c) {return c.multipleBonds;});
            backwardChains = compareBranchesPropertyPositionBackward<unsigned int>(backwardChains, [](Chain c) {return c.doubleBonds;});
            
            while (backwardChains[0].branches.size() > 0){
                backwardChains[0] = addBranch(backwardChains[0], backwardChains[0].branches[0]);
            }
        }

        backwardChains[0] = invertChain(backwardChains[0]);
        bestChains = std::vector<Chain>{backwardChains[0], forwardChains[0]};
    }
    
    else {
        if (chain.branches.size() == 0)
            return chain;
        
        std::vector<Chain> forwardChains = compareBranches(chain.branches);
        std::vector<Chain> backwardChains = chain.branches;
        if (forwardChains.size() != chain.branches.size()) {
            for (int i = 0; i < backwardChains.size(); i++) {
                for (int j = 0; j < forwardChains.size(); j++) {
                    if (backwardChains[i].chain == forwardChains[j].chain) {
                        backwardChains.erase(backwardChains.begin() + i);
                        i--;
                    }
                }
            }
        }
        
        std::vector<Chain> firstChains = backwardChains;
        std::vector<Chain> secondChains = forwardChains;
        
        firstChains = compareBranchesPropertyPositionBackward<unsigned int>(firstChains, [](Chain c) {return c.principalGroups;});
        firstChains = compareBranchesPropertyPositionBackward<unsigned int>(firstChains, [](Chain c) {return c.multipleBonds;});
        firstChains = compareBranchesPropertyPositionBackward<unsigned int>(firstChains, [](Chain c) {return c.doubleBonds;});
        
        secondChains = compareBranchesPropertyPosition<unsigned int>(secondChains, [](Chain c) {return c.principalGroups;});
        secondChains = compareBranchesPropertyPosition<unsigned int>(secondChains, [](Chain c) {return c.multipleBonds;});
        secondChains = compareBranchesPropertyPosition<unsigned int>(secondChains, [](Chain c) {return c.doubleBonds;});

        for (int i = 0; i < firstChains.size(); i++) {
            for (int j = 0; j < secondChains.size(); j++) {
                if (firstChains[i].chain != secondChains[j].chain) {
                    firstChains[i] = addBranch(chain, firstChains[i]);
                    bestChains.push_back(joinChains(firstChains[i], secondChains[j]));
                }
            }
        }

        firstChains = backwardChains;
        firstChains = compareBranchesPropertyPosition<unsigned int>(firstChains, [](Chain c) {return c.principalGroups;});
        firstChains = compareBranchesPropertyPosition<unsigned int>(firstChains, [](Chain c) {return c.multipleBonds;});
        firstChains = compareBranchesPropertyPosition<unsigned int>(firstChains, [](Chain c) {return c.doubleBonds;});

        secondChains = forwardChains;
        secondChains = compareBranchesPropertyPositionBackward<unsigned int>(secondChains, [](Chain c) {return c.principalGroups;});
        secondChains = compareBranchesPropertyPositionBackward<unsigned int>(secondChains, [](Chain c) {return c.multipleBonds;});
        secondChains = compareBranchesPropertyPositionBackward<unsigned int>(secondChains, [](Chain c) {return c.doubleBonds;});

        for (int i = 0; i < firstChains.size(); i++) {
            for (int j = 0; j < secondChains.size(); j++) {
                if (firstChains[i].chain != secondChains[j].chain) {
                    firstChains[i] = addBranch(chain, firstChains[i]);
                    bestChains.push_back(joinChains(firstChains[i], secondChains[j]));
                }
            }
        }
    }

    bestChains = compareBranchesPropertyPosition<unsigned int>(bestChains, [](Chain c) {return c.principalGroups;});
    bestChains = compareBranchesPropertyPosition<unsigned int>(bestChains, [](Chain c) {return c.multipleBonds;});
    bestChains = compareBranchesPropertyPosition<unsigned int>(bestChains, [](Chain c) {return c.doubleBonds;});

    return bestChains[0];
}

Chain Atom::getPrincipalChain (std::vector<Atom*> precedingAtoms, std::map<Atom*, std::vector<unsigned int>> groups, int principalGroup, std::vector<Atom*> ignoredAtoms) {
    std::vector<Chain> possibleChains;

    // get branches
    Chain chain;
    for (int i = 0; i < bonds.size(); i++) {
        if (bonds[i]->species != CARBON) {
            continue;
        }
        if (precedingAtoms.size() > 0 && bonds[i] == precedingAtoms[precedingAtoms.size()-1]) { // double bond behind
            continue;
        }

        bool skip = false;
        for (int j = 0; j < ignoredAtoms.size(); j++) {
            if (bonds[i] == ignoredAtoms[j]) {
                skip = true;
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

        
        // the principal chain is cyclic
        if (precedingAtoms.size() > 0 && precedingAtoms[0] == bonds[i]) {
            chain = Chain();
            chain.chain.push_back(this);
            chain.isCyclic = true;
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
        chain = bonds[i]->getPrincipalChain(currentChain, groups, principalGroup, ignoredAtoms);

        if (chain.isCyclic) {
            possibleChains.clear();
            possibleChains.push_back(chain);
            break;
        }
        if (chain.cyclicGroup != nullptr) {
            if (this == chain.cyclicGroup) {
                chain.cyclicGroup = nullptr;
            }
            else {
                return chain;
            }
        }

        possibleChains.push_back(chain);
    }

    if (possibleChains.size() == 0) {
        chain = Chain();
        possibleChains.push_back(chain);
    }

    std::vector<Chain> branches = compareBranches(possibleChains);

    if (precedingAtoms.size() == 0 && branches.size() == 1) {
        for (int i = 0; i < possibleChains.size(); i++) {
            if (possibleChains[i].chain == branches[0].chain) {
                possibleChains.erase(possibleChains.begin() + i);
                break;
            }
        }
        
        std::vector<Chain> secondaryBranches = compareBranches(possibleChains);
        for (int i = 0; i < secondaryBranches.size(); i++)
            branches.push_back(secondaryBranches[i]);
    }

    // add current atom to chain
    Chain principalChain;
    if (branches.size() == 1) {
        principalChain = branches[0];
    }
    else {
        principalChain.branches = branches;
    }
    principalChain.chain.push_back(this);

    // add principal groups at current atom
    if (groups.find(this) != groups.end()) {
        for (int i = 0; i < groups[this].size(); i++) {
            if (groups[this][i] == principalGroup) {
                principalChain.principalGroups.push_back(precedingAtoms.size());
            }
        }
    }

    return principalChain;
}
