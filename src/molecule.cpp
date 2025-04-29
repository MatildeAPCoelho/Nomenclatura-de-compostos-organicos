#include "molecule.h"
#include "atom.h"
#include <iostream>
#include <string>
#include <vector>

#define LOG(msg) std::cout << __LINE__ << " | " << msg << std::endl

Molecule::Molecule(std::vector<Atom*> atoms) : atoms(atoms) {
}

Molecule::~Molecule() {
}

unsigned int getCarbonyl(std::vector<std::pair<Atom*, unsigned int>> bonds) {
    if (bonds.size() < 2) {
        return ALDEHYDE;
    }

    for (int k = 0; k < bonds.size(); k++) {
        switch (bonds[k].first->species) {
            case OXYGEN:
                if (bonds[k].first->bonds.size() == 2 && bonds[k].first->bonds[0] != bonds[k].first->bonds[1]) {
                    return CARBOALKOXY;
                }
                return CARBOXYL;
            
            case NITROGEN:
                return CARBOXAMIDE;
            
            case FLUORINE:
            case CHLORINE:
            case BROMINE:
            case IODINE:
                return HALOFORMYL;
            
            default:
                break;
        }
    }
    return KETONE;
}

std::map<Atom*, std::vector<unsigned int>> Molecule::getFunctionalGroups() {
    std::map<Atom*, std::vector<unsigned int>> groups;

    for (int i = 0; i < atoms.size(); i++) {
        if (atoms[i]->species != CARBON)
            continue;
        
        std::vector<std::pair<Atom*, unsigned int>> bonds;

        bool added = false;
        for (int j = 0; j < atoms[i]->bonds.size(); j++) {
            added = false;
            for (int k = 0; k < bonds.size(); k++) {
                if (atoms[i]->bonds[j] == bonds[k].first) {
                    bonds[k].second += 1;
                    added = true;
                    break;
                }
            }
            if (!added) {
                bonds.push_back(std::pair{atoms[i]->bonds[j], 1});
            }
        }

        // =O
        for (int j = 0; j < bonds.size(); j++) {
            if (bonds[j].second == 2 && bonds[j].first->species == OXYGEN) {
                bonds.erase(bonds.cbegin() + j);
                groups[atoms[i]] = std::vector<unsigned int>{getCarbonyl(bonds)};
                break;
            }
            if (bonds[j].second == 3 && bonds[j].first->species == NITROGEN) {
                groups[atoms[i]] = std::vector<unsigned int>{NITRILE};
                break;
            }
        }
        if (groups.find(atoms[i]) != groups.end()) {
            continue;
        }

        // -O / -N / -S

        for (int j = 0; j < bonds.size(); j++) {
            if (bonds[j].second == 1) {
                switch (bonds[j].first->species) {
                    case NITROGEN:
                        groups[atoms[i]].push_back(AMINE);
                        break;
                    
                    case OXYGEN:
                        if (bonds[j].first->bonds.size() == 1) {
                            groups[atoms[i]].push_back(HYDROXYL);
                        }
                        break;
                    
                    case SULFUR:
                        if (bonds[j].first->bonds.size() == 1) {
                            groups[atoms[i]].push_back(SULFHYDRYL);
                        }
                        break;
                    
                    default:
                        break;
                }
            }
        }
    }
    return groups;
}

std::map<int, std::string> greekPrefix = {
    {1, ""},
    {2, "di"},
    {3, "tri"},
    {4, "tetra"},
    {5, "penta"},
    {6, "hexa"},
    {7, "hepta"},
    {8, "octa"},
    {9, "nona"},
    {10, "deca"}
};

std::string nameCarbonChain(Chain chain, bool isSubstituent = false) {
    std::string name;

    std::map<int, std::string> length = {
        {1, "met"},
        {2, "et"},
        {3, "prop"},
        {4, "but"},
        {5, "pent"},
        {6, "hex"},
        {7, "hept"},
        {8, "oct"},
        {9, "non"},
        {10, "dec"},
        {11, "undec"},
        {12, "dodec"},
        {13, "tridec"},
        {20, "icos"}
        // rules...
    };

    if (length.find(chain.chain.size()) != length.end()) {
        name = length[chain.chain.size()];
    }
    else {
        name = std::to_string(chain.chain.size());
    }

    if (chain.isCyclic) {
        if (name[0] == 'h') {
            name = "ciclo-" + name;
        }
        else {
            name = "ciclo" + name;
        }
    }


    if (chain.multipleBonds.size() == 0) {
        name += "an";
    }
    else {
        // double bonds
        if (chain.doubleBonds.size() > 0) {
            if (chain.doubleBonds.size() > 1) {
                name += "a";
            }
            name += "-" + std::to_string(chain.doubleBonds[0]+1);

            for (int i = 1; i < chain.doubleBonds.size(); i++) {
                name += "," + std::to_string(chain.doubleBonds[i]+1);
            }

            if (greekPrefix.find(chain.doubleBonds.size()) != greekPrefix.end()) {
                name += "-" + greekPrefix[chain.doubleBonds.size()] + "en";
            }
            else {
                name += "-" + std::to_string(chain.doubleBonds.size()) + "en";
            }
        }

        // triple bonds
        std::vector<int> tripleBonds;
        for (int i = 0; i < chain.multipleBonds.size(); i++) {
            if (chain.multipleBonds.size() == chain.doubleBonds.size()) {
                break;
            }
            for (int j = 0; j < chain.doubleBonds.size(); j++) {
                if (chain.multipleBonds[i] == chain.doubleBonds[j]) {
                    goto loopiEnd;
                }
            }
            tripleBonds.push_back(chain.multipleBonds[i]);
            loopiEnd:;
        }
        if (tripleBonds.size() > 0) {
            if (tripleBonds.size() > 1) {
                if (chain.doubleBonds.size() > 0)
                    name += "o";
                else
                    name += "a";
            }

            name += "-" + std::to_string(tripleBonds[0]+1);
            for (int i = 1; i < tripleBonds.size(); i++) {
                name += "," + std::to_string(tripleBonds[i]+1);
            }

            if (greekPrefix.find(tripleBonds.size()) != greekPrefix.end()) {
                name += "-" + greekPrefix[tripleBonds.size()] + "in";
            }
            else {
                name += "-" + std::to_string(tripleBonds.size()) + "in";
            }
        }
    }

    return name;
}

std::string getChainPrefix(Chain chain) {
    std::map<std::string, std::vector<unsigned int>> prefixes;
    
    // getFunctionalGroups

    // look for substituents
    for (int i = 0; i < chain.chain.size(); i++) {
        for (int j = 0; j < chain.chain[i]->bonds.size(); j++) {
            bool skip = false;
            for (int k = 0; k < chain.chain.size(); k++) {
                if (chain.chain[i]->bonds[j] == chain.chain[k]) {
                    skip = true;
                    break;
                }
            }
            if (skip)
                continue;

            switch (chain.chain[i]->bonds[j]->species) {
                // secondary chain
                case CARBON:
                    // getPrincipalChain
                    // getPrefix
                    // -il-
                    break;

                case OXYGEN: // ester??
                    // eter
                    if (chain.chain[i]->bonds[j]->bonds.size() == 2) {
                        // getPrincipalChain
                        // getPrefix
                        // -oxi-
                    }
                    break;

                case SULFUR:
                    // sulfanil
                    if (chain.chain[i]->bonds[j]->bonds.size() == 2) {
                        // getPrincipalChain
                        // getPrefix
                        // -sulfanil-
                    }
                    break;
                
                case NITROGEN: // amide??
                    // amine
                    for (int i = 0; i < chain.chain[i]->bonds[j]->bonds.size(); i++) {
                        // getPrincipalChain
                        // getPrefix
                    }
                    // -amino-
                    break;

                // halogens
                case FLUORINE:
                    // fluoro-
                    break;
                case CHLORINE:
                    // cloro-
                    break;
                case BROMINE:
                    // bromo-
                    break;
                case IODINE:
                    // iodo-
                    break;
                default:
                    break;
            }
        }
    }
    
    // order prefixes
    // functional group ->
    // not functional group -> getPrincipalChain (principalGroup = 0; ignore atom from other chain); getPrefix + carbon + "il"
    std::string suffix = nameCarbonChain(chain);

    return suffix;
}

std::string Molecule::getName() {
    // find functional groups
    std::map<Atom*, std::vector<unsigned int>> groups = getFunctionalGroups();

    // print functional groups
    std::cout << "-----" << std::endl;
    for (std::map<Atom*, std::vector<unsigned int> >::const_iterator it = groups.begin(); it != groups.end(); ++it) {
        std::cout << it->first->bonds.size() << " : ";
        for (int j = 0; j < it->second.size(); j++) {
            std::cout << it->second[j] << " ";
        }
        std::cout << std::endl;
    }
    //

    // determine principal group
    unsigned int principalGroup = 0;
    std::vector<Atom*> principalGroupAtoms;
    for (std::map<Atom*, std::vector<unsigned int> >::const_iterator it = groups.begin(); it != groups.end(); ++it) {
        for (int i = 0; i < it->second.size(); i++) {
            if (principalGroup == it->second[i])
                principalGroupAtoms.push_back(it->first);

            else if (principalGroup > it->second[i] || principalGroup == 0) {
                principalGroup = it->second[i];
                principalGroupAtoms.clear();
                principalGroupAtoms.push_back(it->first);
            }
        }
    }

    // get possible chains with principal group
    std::vector<Chain> possibleChains;
    for (int i = 0; i < principalGroupAtoms.size(); i++) {
        Chain chain = principalGroupAtoms[i]->getPrincipalChain(std::vector<Atom*>(), groups, principalGroup);
        chain = getOrderedChain(chain);
        possibleChains.push_back(chain);
    }
    // carbox-
    bool carbox = false;
    for (int i = 0; i < possibleChains.size(); i++) {
        if (possibleChains[i].chain.size() > 1) {
            break;
        }
        if (i == possibleChains.size() - 1) {
            carbox = true;
        }
    }
    if (carbox) {
        possibleChains.clear();
    }
    if (principalGroup == CARBOXYL || principalGroup == CARBOALKOXY || principalGroup == HALOFORMYL || principalGroup == CARBOXAMIDE || principalGroup == NITRILE) {
        std::vector<Atom*> ignoredAtoms;
        for (int i = 0; i < principalGroupAtoms.size(); i++) {
            for (int j = 0; j < principalGroupAtoms[i]->bonds.size(); j++) {
                if (principalGroupAtoms[i]->bonds[j]->species == CARBON) {
                    groups.erase(principalGroupAtoms[i]);
                    groups[principalGroupAtoms[i]->bonds[j]].push_back(principalGroup);
                    
                    ignoredAtoms.push_back(principalGroupAtoms[i]);

                    principalGroupAtoms[i] = principalGroupAtoms[i]->bonds[j];
                    break;
                }
            }
        }
        
        for (int i = 0; i < principalGroupAtoms.size(); i++) {
            Chain chain = principalGroupAtoms[i]->getPrincipalChain(std::vector<Atom*>(), groups, principalGroup, ignoredAtoms);
            chain = getOrderedChain(chain);
            possibleChains.push_back(chain);
        }
    }

    // possible chains without principle group
    if (possibleChains.size() == 0) {
        for (int i = 0; i < atoms.size(); i++) {
            if (atoms[i]->species == CARBON) {
                Chain chain = atoms[i]->getPrincipalChain(std::vector<Atom*>(), groups, principalGroup);
                chain = getOrderedChain(chain);
                possibleChains.push_back(chain);
            }
        }
    }
    
    // determine principal chain
    principalChain = compareChains(possibleChains)[0];
    
    // name substituents
    std::string name = getChainPrefix(principalChain);

    // get suffix
    
    return name;
}