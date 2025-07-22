#include "molecule.h"
#include "atom.h"

#include <cctype>
#include <string>
#include <vector>
#include <map>

Molecule::Molecule(std::vector<Atom*> atoms) : atoms(atoms) {
}

Molecule::~Molecule() {
}

// function initialization
std::string getRadicalName(Chain chain, std::vector<Atom*> principalGroupAtoms, unsigned int principalGroup, std::vector<Atom*> ignoredAtoms = std::vector<Atom*>(), std::vector<Atom*> carboxyAtoms = std::vector<Atom*>());
std::vector<Chain> compareChains(std::vector<Chain> chains, std::vector<Atom*> principalGroupAtoms, unsigned int principalGroup, std::vector<Atom*> ignoredAtoms = std::vector<Atom*>());
Chain getOrderedChain(Chain chain, std::vector<Atom*> principalGroupAtoms, unsigned int principalGroup, std::vector<Atom*> ignoredAtoms = std::vector<Atom*>());

std::pair<unsigned int, unsigned int> getCarbonyl(std::vector<std::pair<Atom*, unsigned int>> bonds) {
    int carbonNumber = 0;
    std::pair<unsigned int, unsigned int> group = {0, 0};

    for (int k = 0; k < bonds.size(); k++) {
        switch (bonds[k].first->species) {
            case CARBON:
                carbonNumber += 1;
                break;
            
            case OXYGEN:
                if (bonds[k].first->bonds.size() == 2 && bonds[k].first->bonds[0] != bonds[k].first->bonds[1]) {
                    if (group.first == 0 || group.first > CARBOALKOXY) {
                        group = {CARBOALKOXY, 0};
                    }
                    break;
                }
                return {CARBOXYL, 0};
            
            case NITROGEN:
                if (group.first == 0 || group.first > CARBOXAMIDE) {
                    group = {CARBOXAMIDE, 0};
                }
                break;
            
            case FLUORINE:
            case CHLORINE:
            case BROMINE:
            case IODINE:
                if (group.first == 0 || group.first > HALOFORMYL) {
                    return {HALOFORMYL, bonds[k].first->species};
                }
                break;
            
            default:
                break;
        }
    }

    if (group.first == 0) {
        if (carbonNumber == 2)
            return {KETONE, 0};
        else
            return {ALDEHYDE, 0};
    }

    return group;
}

std::vector<std::pair<Atom*, unsigned int>> getBondTypes(Atom* atom, std::vector<Atom*> ignoredAtoms = {}) {
    std::vector<std::pair<Atom*, unsigned int>> bonds;

    for (int i = 0; i < atom->bonds.size(); i++) {
        for (int j = 0; j < ignoredAtoms.size(); j++) {
            if (atom->bonds[i] == ignoredAtoms[j]) {
                goto loopEnd;
            }
        }
        for (int j = 0; j < bonds.size(); j++) {
            if (atom->bonds[i] == bonds[j].first) {
                bonds[j].second += 1;
                goto loopEnd;
            }
        }
        bonds.push_back(std::pair{atom->bonds[i], 1});
        loopEnd:;
    }

    return bonds;
}

struct principalGroupInfo {
    std::vector<Atom*> principalGroupAtoms;
    unsigned int principalGroup;
    unsigned int halogen;  
};

principalGroupInfo getPrincipalGroup(std::vector<Atom*> atoms) {
    std::vector<Atom*> principalGroupAtoms = {};
    unsigned int principalGroup = 0;
    unsigned int halogen = 0;

    for (int i = 0; i < atoms.size(); i++) {
        if (atoms[i]->species != CARBON)
            continue;
        
        std::vector<std::pair<Atom*, unsigned int>> bonds = getBondTypes(atoms[i]);

        // =O
        unsigned int group = 0;
        for (int j = 0; j < bonds.size(); j++) {
            if (bonds[j].second == 2 && bonds[j].first->species == OXYGEN) {
                bonds.erase(bonds.begin() + j);
                std::pair<unsigned int, unsigned int> carbonyl = getCarbonyl(bonds);
                group = getCarbonyl(bonds).first;
                
                if (group == HALOFORMYL) {
                    if (principalGroup == 0 || principalGroup > HALOFORMYL || halogen > carbonyl.second && principalGroup  == HALOFORMYL) {
                        halogen = carbonyl.second;
                        principalGroupAtoms.clear();
                    }
                    else if (halogen < carbonyl.second) {
                        group = 0;
                        break;
                    }
                }
                
                break;
            }
            if (bonds[j].second == 3 && bonds[j].first->species == NITROGEN) {
                group = NITRILE;
                break;
            }
        }
        if (group != 0) {
            if (group == principalGroup) {
                principalGroupAtoms.push_back(atoms[i]);
            }
            else if (group < principalGroup || principalGroup == 0) {
                principalGroupAtoms = {atoms[i]};
                principalGroup = group;
            }
            continue;
        }

        // -O / -N / -S

        for (int j = 0; j < bonds.size(); j++) {
            group = 0;
            if (bonds[j].second == 1) {
                switch (bonds[j].first->species) {
                    case NITROGEN:
                        group = AMINE;
                        break;
                    
                    case OXYGEN:
                        if (bonds[j].first->bonds.size() == 1) {
                            group = HYDROXYL;
                        }
                        break;
                    
                    case SULFUR:
                        if (bonds[j].first->bonds.size() == 1) {
                            group = SULFHYDRYL;
                        }
                        break;
                    
                    default:
                        break;
                }
            }

            if (group != 0) {
                if (group == principalGroup) {
                    principalGroupAtoms.push_back(atoms[i]);
                }
                else if (group < principalGroup || principalGroup == 0) {
                    principalGroupAtoms = {atoms[i]};
                    principalGroup = group;
                }
            }
        }
    }

    principalGroupInfo info;
    info.principalGroupAtoms = principalGroupAtoms;
    info.principalGroup = principalGroup;
    info.halogen = halogen;
    return info;
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
            if (chain.chain.size() == 2 || chain.isCyclic && chain.multipleBonds.size() == 1 && chain.substituents.size() == 0) {
                name += "en";
            }
            else {
                if (chain.doubleBonds.size() > 1) {
                    name += "a";
                }
                
                name += "-" + std::to_string(chain.doubleBonds[chain.doubleBonds.size()-1]+1);

                for (int i = (int)chain.doubleBonds.size()-2; i >= 0; i--) {
                    name += "," + std::to_string(chain.doubleBonds[i]+1);
                }

                if (greekPrefix.find(chain.doubleBonds.size()) != greekPrefix.end()) {
                    name += "-" + greekPrefix[chain.doubleBonds.size()] + "en";
                }
                else {
                    name += "-" + std::to_string(chain.doubleBonds.size()) + "en";
                }
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
            if (tripleBonds.size() > 1 || chain.doubleBonds.size() > 0) {
                if (chain.doubleBonds.size() > 0)
                    name += "o";
                else
                    name += "a";
            }

            if (chain.chain.size() == 2) {
                name += "in";
            }
            else {
                name += "-" + std::to_string(tripleBonds[tripleBonds.size()-1]+1);
                for (int i = (int)tripleBonds.size()-2; i >= 0; i--) {
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
    }

    return name;
}

std::vector<Atom*> crossBridge(Atom* initialAtom, Atom* bridgeAtom) {
    std::vector<Atom*> atoms;

    for (int i = 0; i < bridgeAtom->bonds.size(); i++) {
        if (bridgeAtom->bonds[i] != initialAtom) {
            atoms.push_back(bridgeAtom->bonds[i]);
        }
    }
    return atoms;
}

int compareStrings(std::string one, std::string two) {
    int ii = 0;
    for (int i = 0; i < one.size(); i++, ii++) {
        while (!std::isalpha(one[i]) || one[i] == 'N') {
            i++;
        }
        while (!std::isalpha(two[ii]) || two[ii] == 'N') {
            ii++;
        }
        
        if (ii == two.size() || one[i] > two[ii]) {
            return 1;
        }
        if (one[i] < two[ii] || i == one.size()-1 && ii < two.size()-1) {
            return -1;
        }
    }
    return 0;
}

std::vector<std::string> orderStrings(std::vector<std::string> strings) {
    std::vector<std::string> orderedStrings;
    orderedStrings.push_back(strings[0]);

    for (int i = 1; i < strings.size(); i++) {
        for (int j = 0; j < orderedStrings.size(); j++) {
            if (compareStrings(strings[i], orderedStrings[j]) < 0) {
                orderedStrings.insert(orderedStrings.begin() + j, strings[i]);
                break;
            }
            if (j == orderedStrings.size()-1) {
                orderedStrings.push_back(strings[i]);
                break;
            }
        }
    }
    return orderedStrings;
}

std::vector<std::string> orderPrefixes(std::map<std::string, std::vector<std::string>> prefixes) {
    std::vector<std::string> orderedPrefixes;
    for (auto const& it : prefixes) {
        for (int i = 0; i < orderedPrefixes.size(); i++) {
            if (compareStrings(it.first, orderedPrefixes[i]) < 0) {
                orderedPrefixes.insert(orderedPrefixes.begin() + i, it.first);
                break;
            }
            if (i == orderedPrefixes.size()-1) {
                orderedPrefixes.push_back(it.first);
                break;
            }
        }
        if (orderedPrefixes.size() == 0) {
            orderedPrefixes.push_back(it.first);
        }
    }
    return orderedPrefixes;
}

std::string getSideChainPrefix(Atom* initialAtom, std::vector<Atom*> ignoredAtoms) {
    std::vector<Atom*> secondaryChainGroups = {initialAtom};
    std::vector<Atom*> carboxyAtoms;

    Chain sideChain = initialAtom->getPrincipalChain(std::vector<Atom*>(), secondaryChainGroups, CHAIN, ignoredAtoms);
    sideChain = getOrderedChain(sideChain, secondaryChainGroups, CHAIN, ignoredAtoms);
    principalGroupInfo info = getPrincipalGroup(sideChain.chain);
    
    if (info.principalGroup == CARBOXYL && sideChain.chain.size() == 1) {
        return "carboxi";
    }

    while (info.principalGroup == CARBOXYL) {
        ignoredAtoms.insert(ignoredAtoms.end(), info.principalGroupAtoms.begin(), info.principalGroupAtoms.end());
        for (int k = 0; k < info.principalGroupAtoms.size(); k++) {
            for (int l = 0; l < info.principalGroupAtoms[k]->bonds.size(); l++) {
                if (info.principalGroupAtoms[k]->bonds[l]->species == CARBON) {
                    carboxyAtoms.push_back(info.principalGroupAtoms[k]->bonds[l]);
                    break;
                }
            }
        }

        sideChain = initialAtom->getPrincipalChain(std::vector<Atom*>(), secondaryChainGroups, CHAIN, ignoredAtoms);
        sideChain = getOrderedChain(sideChain, secondaryChainGroups, CHAIN, ignoredAtoms);
        info = getPrincipalGroup(sideChain.chain);
    }

    std::string prefix = getRadicalName(sideChain, secondaryChainGroups, 0, ignoredAtoms, carboxyAtoms);

    return prefix;
}

std::map<std::string, std::vector<std::string>> getChainPrefixes(Chain chain, std::vector<Atom*> principalGroupAtoms, unsigned int principalGroup, std::vector<Atom*> ignoredAtoms = std::vector<Atom*>()) {
    std::map<std::string, std::vector<std::string>> prefixes;
    Chain sideChain;
    std::string prefix;
    
    // name substituents
    for (int i = (int)chain.chain.size()-1; i >= 0; i--) {
        std::vector<std::pair<Atom*, unsigned int>> bonds;
        Atom* atom;
        std::vector<Atom*> atoms;
        std::vector<Atom*> secondaryChainGroups;

        bonds = getBondTypes(chain.chain[i], ignoredAtoms);

        if (principalGroup < KETONE && principalGroup != CHAIN && !chain.carbox) {
            for (int j = 0; j < principalGroupAtoms.size(); j++) {
                if (principalGroupAtoms[j] == chain.chain[i]) {
                    bonds.clear();
                    break;
                }
            }
        }

        std::vector<Atom*> sideChainIgnoredAtoms;
        std::vector<Atom*> carboxyAtoms;
        principalGroupInfo info;
        for (int j = 0; j < bonds.size(); j++) {
            for (int k = 0; k < chain.chain.size(); k++) {
                if (bonds[j].first == chain.chain[k]) {
                    goto loopjEnd;
                }
            }
            for (int k = 0; k < ignoredAtoms.size(); k++) {
                if (bonds[j].first == ignoredAtoms[k]) {
                    goto loopjEnd;
                }
            }
            if (principalGroup < KETONE && principalGroup != CHAIN && chain.carbox) {
                if (getPrincipalGroup({bonds[j].first}).principalGroup == principalGroup) {
                    continue;
                }
            }

            switch (bonds[j].first->species) {
                case CARBON:
                    // secondary chain
                    prefix = getSideChainPrefix(bonds[j].first, {chain.chain[i]});
                    prefixes[prefix].push_back(std::to_string(chain.chain.size()-i));
                    break;

                case OXYGEN:
                    if (bonds[j].second == 1) {
                        // eter
                        if (bonds[j].first->bonds.size() == 2) {
                            atom = crossBridge(chain.chain[i], bonds[j].first)[0];

                            std::string sideChainPrefix = getSideChainPrefix(atom, {bonds[j].first});

                            std::vector<std::string> carbonStrings = {"metil", "etil", "propil", "butil", "fenil"};
                            for (int k = 0; k < carbonStrings.size(); k++) {
                                if (sideChainPrefix.size() >= carbonStrings[k].size() && sideChainPrefix.substr(sideChainPrefix.size()-carbonStrings[k].size(), carbonStrings[k].size()) == carbonStrings[k]) {
                                    sideChainPrefix = sideChainPrefix.substr(0, sideChainPrefix.size()-2);
                                    prefix = sideChainPrefix + "oxi";
                                    goto skipLoopk;
                                }
                                if (sideChainPrefix.size() > carbonStrings[k].size() && sideChainPrefix.substr(sideChainPrefix.size()-carbonStrings[k].size()-1, carbonStrings[k].size()) == carbonStrings[k]) {
                                    sideChainPrefix = sideChainPrefix.substr(0, sideChainPrefix.size()-3);
                                    prefix = sideChainPrefix + "oxi)";
                                    goto skipLoopk;
                                }
                            }
                            prefix = "(" + sideChainPrefix + "oxi)";
                            skipLoopk:;

                            prefixes[prefix].push_back(std::to_string(chain.chain.size()-i));
                        }
                        // hydroxyl
                        else {
                            if (principalGroup != HYDROXYL) {
                                prefixes["hidroxi"].push_back(std::to_string(chain.chain.size()-i));
                            }
                        }
                    }
                    // ketone
                    else if (bonds[j].second == 2 && principalGroup != KETONE) {
                        prefixes["oxo"].push_back(std::to_string(chain.chain.size()-i));
                    }
                    break;
                
                case NITROGEN:
                    // amine
                    if (bonds[j].second == 1 && principalGroup != AMINE) {
                        // primary
                        if (bonds[j].first->bonds.size() == 1) {
                            prefixes["amino"].push_back(std::to_string(chain.chain.size()-i));
                            break;
                        }

                        // secondary
                        atoms = crossBridge(chain.chain[i], bonds[j].first);

                        std::string sideChainPrefix0 = getSideChainPrefix(atoms[0], {bonds[j].first});
                        
                        if (atoms.size() == 1) {
                            prefix = "(" + sideChainPrefix0 + "amino)";
                            prefixes[prefix].push_back(std::to_string(chain.chain.size()-i));
                            break;
                        }

                        // tertiary
                        std::string sideChainPrefix1 = getSideChainPrefix(atoms[1], {bonds[j].first});
                        
                        int kk = 0;
                        for (int k = 0; k < sideChainPrefix1.size(); k++, kk++) {
                            while (!std::isalpha(sideChainPrefix0[k])) {
                                k++;
                            }
                            while (!std::isalpha(sideChainPrefix1[kk])) {
                                kk++;
                            }
                            
                            // equal prefixes
                            if (k == sideChainPrefix0.size()-1 && kk == sideChainPrefix1.size()-1) {
                                prefix = "(" + greekPrefix[2] + sideChainPrefix0 + "amino)";
                                break;
                            }
                            if (k == sideChainPrefix0.size()-1 || sideChainPrefix0[k] < sideChainPrefix1[kk]) {
                                if (sideChainPrefix1[0] != '(') {
                                    sideChainPrefix1 = "(" + sideChainPrefix1 + ")";
                                }
                                prefix = "(" + sideChainPrefix0 + sideChainPrefix1 + "amino)";
                                break;
                            }
                            if (kk == sideChainPrefix1.size()-1 || sideChainPrefix0[k] > sideChainPrefix1[kk]) {
                                if (sideChainPrefix0[0] != '(') {
                                    sideChainPrefix0 = "(" + sideChainPrefix0 + ")";
                                }
                                prefix = "(" + sideChainPrefix1 + sideChainPrefix0 + "amino)";
                                break;
                            }
                        }
                        prefixes[prefix].push_back(std::to_string(chain.chain.size()-i));
                    }
                    // nitrile
                    else if (bonds[j].second == 3) {
                        prefixes["ciano"].push_back(std::to_string(chain.chain.size()-i)); // not really...
                    }
                    break;

                case SULFUR:
                    // sulfanil
                    if (bonds[j].first->bonds.size() == 2) {
                        if (bonds[j].first->bonds.size() == 2) {
                            atom = crossBridge(chain.chain[i], bonds[j].first)[0];

                            std::string sideChainPrefix = getSideChainPrefix(atom, {bonds[j].first});

                            prefix = "(" + sideChainPrefix + "sulfanil)";

                            prefixes[prefix].push_back(std::to_string(chain.chain.size()-i));
                        }
                    }
                    // thiol
                    else {
                        if (principalGroup != SULFHYDRYL) {
                            prefixes["sulfanil"].push_back(std::to_string(chain.chain.size()-i));
                        }
                    }
                    break;

                // halogens
                case FLUORINE:
                    prefixes["fluoro"].push_back(std::to_string(chain.chain.size()-i));
                    break;
                case CHLORINE:
                    prefixes["cloro"].push_back(std::to_string(chain.chain.size()-i));
                    break;
                case BROMINE:
                    prefixes["bromo"].push_back(std::to_string(chain.chain.size()-i));
                    break;
                case IODINE:
                    prefixes["iodo"].push_back(std::to_string(chain.chain.size()-i));
                    break;
                default:
                    break;
            }

            loopjEnd:;
        }
    }
    return prefixes;
}

std::string getRadicalName(Chain chain, std::vector<Atom*> principalGroupAtoms, unsigned int principalGroup, std::vector<Atom*> ignoredAtoms, std::vector<Atom*> carboxyAtoms) {
    std::map<std::string, std::vector<std::string>> prefixes;
    unsigned int yl = 0;

    Chain sideChain;
    std::string prefix;

    for (int i = (int)chain.chain.size()-1; i >= 0; i--) {
        std::vector<std::pair<Atom*, unsigned int>> bonds;
        Atom* atom;
        std::vector<Atom*> atoms;

        // -yl
        if (principalGroup == CHAIN) {
            for (int j = 0; j < principalGroupAtoms.size(); j++) {
                if (principalGroupAtoms[j] == chain.chain[i]) {
                    yl = chain.chain.size()-i;
                    break;
                }
            }
        }
    }

    prefixes = getChainPrefixes(chain, principalGroupAtoms, principalGroup, ignoredAtoms);

    // carboxy
    for (int i = (int)chain.chain.size()-1; i >= 0; i--) {
        for (int j = 0; j < carboxyAtoms.size(); j++) {
            if (chain.chain[i] == carboxyAtoms[j]) {
                prefixes["carboxi"].push_back(std::to_string(chain.chain.size() - i));        
            }
        }
    }


    // order prefixes
    std::vector<std::string> orderedPrefixes = orderPrefixes(prefixes);

    for (int i = 0; i < orderedPrefixes.size(); i++) {
        std::string prefixName = "";
        if (chain.chain.size() > 1) {
            prefixName += prefixes[orderedPrefixes[i]][0];
            for (int j = 1; j < prefixes[orderedPrefixes[i]].size(); j++) {
                prefixName += "," + prefixes[orderedPrefixes[i]][j];
            }
            prefixName += "-";
        }
        if (orderedPrefixes[i][0] == 'h' && prefixes[orderedPrefixes[i]].size() > 1) {
            orderedPrefixes[i] = prefixName + greekPrefix[prefixes[orderedPrefixes[i]].size()] + "-" + orderedPrefixes[i];
        }
        else {
            orderedPrefixes[i] = prefixName + greekPrefix[prefixes[orderedPrefixes[i]].size()] + orderedPrefixes[i];
        }
    }

    std::string name;
    if (orderedPrefixes.size() > 0) {
        name += orderedPrefixes[0];
        for (int i = 1; i < orderedPrefixes.size(); i++) {
            if (chain.chain.size() > 1) {
                name += "-";
            }
            name += orderedPrefixes[i];
        }
    }

    std::string carbonChain = nameCarbonChain(chain);
    if (carbonChain[0] == 'h') {
        carbonChain = '-' + carbonChain;
    }
    name += carbonChain;

    if (name.substr(name.size()-2, 2) == "an" && yl == 1) {
        name = name.substr(0, name.size()-2);
    }
    
    if (yl != 1 || chain.multipleBonds.size() > 0 && (chain.chain.size() != 2 || chain.substituents.size() > 1)) {
        name += "-" + std::to_string(yl) + "-";
    }

    name += "il";

    // fenil
    if (name.size() >= std::string("ciclo-hexa-1,3,5-trien-1-il").size() && 
    (name.substr(name.size() - std::string("ciclo-hexa-1,3,5-trien-1-il").size(), std::string("ciclo-hexa-1,3,5-trien-1-il").size()) == "ciclo-hexa-1,3,5-trien-1-il" || 
    name.substr(name.size() - std::string("ciclo-hexa-1,3,5-trien-1-il").size(), std::string("ciclo-hexa-1,3,5-trien-1-il").size()) == "ciclo-hexa-2,4,6-trien-1-il")) {
        name = name.substr(0, name.size() - std::string("ciclo-hexa-1,3,5-trien-1-il").size()) + "fenil";
    }
    
    // carbonil / formil / carbamoil
    if (name.size() >= std::string("oxometil").size() && name.substr(name.size() - std::string("oxometil").size(), std::string("oxometil").size()) == "oxometil"){
        if (name == "oxometil") {
            name = "formil";
            return name;
        }
        name = name.substr(0, name.size() - std::string("oxometil").size());

        if (name == "amino") {
            name = "carbamoil";
            return name;
        }

        if (name.size() >= std::string("amino)").size() && name.substr(name.size() - std::string("amino)").size(), std::string("amino").size()) == "amino") {
            name = name.substr(1, name.size() - std::string("amino").size() - 2) + "carbamoil";
        }
        else {
            name += "carbonil";
        }
    }
    
    if (prefixes.size() > 0 || chain.multipleBonds.size() > 0 && chain.chain.size() != 2 && name != "fenil" || yl != 1) {
        name = "(" + name + ")";
    }

    if (name == "(cianometil)") {
        name = "ciano";
    }

    return name;
}

std::map<std::string, std::vector<std::string>> getOxyPrefix(std::map<std::string, std::vector<std::string>> prefixes, Atom* originAtom, std::string position) {
    std::vector<std::pair<Atom*, unsigned int>> bonds = getBondTypes(originAtom, std::vector<Atom*>{});
    std::vector<Atom*> secondaryChainGroups;

    for (int i = 0; i < bonds.size(); i++) {
        if (bonds[i].first->species == OXYGEN && bonds[i].second == 1) {
            std::vector<Atom*> atoms = crossBridge(originAtom, bonds[i].first);
            
            for (int j = 0; j < atoms.size(); j++) {
                secondaryChainGroups.push_back(atoms[j]);
                Chain sideChain = atoms[j]->getPrincipalChain(std::vector<Atom*>(), secondaryChainGroups, CHAIN, std::vector<Atom*>{bonds[i].first});
                sideChain = getOrderedChain(sideChain, secondaryChainGroups, CHAIN, std::vector<Atom*>{bonds[i].first});
                std::string sideChainPrefix = getRadicalName(sideChain, secondaryChainGroups, CHAIN, std::vector<Atom*>{bonds[i].first});

                prefixes[sideChainPrefix].insert(prefixes[sideChainPrefix].begin(), 1, position);
            }
        }
    }

    return prefixes;
}

std::vector<std::vector<std::string>> getAminePrefix(Atom* originAtom) {
    std::vector<std::vector<std::string>> prefixes;
    std::vector<std::pair<Atom*, unsigned int>> bonds = getBondTypes(originAtom, std::vector<Atom*>{});
    std::vector<Atom*> secondaryChainGroups;

    for (int i = 0; i < bonds.size(); i++) {
        if (bonds[i].first->species == NITROGEN && bonds[i].second == 1) {
            std::vector<Atom*> atoms = crossBridge(originAtom, bonds[i].first);
            
            if(atoms.size() > 0)
                prefixes.push_back({});

            for (int j = 0; j < atoms.size(); j++) {
                secondaryChainGroups.push_back(atoms[j]);
                Chain sideChain = atoms[j]->getPrincipalChain(std::vector<Atom*>(), secondaryChainGroups, CHAIN, std::vector<Atom*>{bonds[i].first});
                sideChain = getOrderedChain(sideChain, secondaryChainGroups, CHAIN, std::vector<Atom*>{bonds[i].first});
                std::string sideChainPrefix = getRadicalName(sideChain, secondaryChainGroups, CHAIN, std::vector<Atom*>{bonds[i].first});

                prefixes[prefixes.size()-1].push_back(sideChainPrefix);
            }
        }
    }

    return prefixes;
}

std::map<std::string, std::vector<std::string>> addAminePrefixes(std::map<std::string, std::vector<std::string>> prefixes, std::vector<std::pair<std::string, std::vector<std::vector<std::string>>>> aminePrefixes) {
    int maxSize = 0;

    // order amine groups
    for (int i = 0; i < aminePrefixes.size(); i++) {
        for (int j = 0; j < aminePrefixes[i].second.size(); j++) {
            aminePrefixes[i].second[j] = orderStrings(aminePrefixes[i].second[j]);
        }
        
        if (aminePrefixes[i].second.size() > maxSize) {
            maxSize = aminePrefixes[i].second.size();
        }
        
        for (int j = 0; j < aminePrefixes[i].second.size(); j++) {
            for (int k = 0; k < j; k++) {
                if (aminePrefixes[i].second[j].size() < aminePrefixes[i].second[k].size() || aminePrefixes[i].second[j] == aminePrefixes[i].second[k]) {
                    aminePrefixes[i].second.insert(aminePrefixes[i].second.begin() + k, aminePrefixes[i].second.begin() + j, aminePrefixes[i].second.begin() + j + 1);
                    aminePrefixes[i].second.erase(aminePrefixes[i].second.begin() + j + 1);
                    break;
                }
                for (int l = 0; l < aminePrefixes[i].second[j].size(); l++) {
                    if (compareStrings(aminePrefixes[i].second[j][l], aminePrefixes[i].second[k][l]) < 0) {
                        aminePrefixes[i].second.insert(aminePrefixes[i].second.begin() + k, aminePrefixes[i].second.begin() + j, aminePrefixes[i].second.begin() + j + 1);
                        aminePrefixes[i].second.erase(aminePrefixes[i].second.begin() + j + 1);
                        goto loopjEnd;        
                    }
                }
            }
            loopjEnd:;
        }
    }

    // add amine prefixes to prefixes map
    for (int i = 0; i < maxSize; i++) {
        for (int j = aminePrefixes.size()-1; j >= 0; j--) {
            std::string positionName = aminePrefixes[j].first;
            positionName = positionName[0] + std::string(i, '\'') + positionName.substr(1, positionName.size()-1);

            if (aminePrefixes[j].second.size() > i) {
                for (int k = 0; k < aminePrefixes[j].second[i].size(); k++) {
                    prefixes[aminePrefixes[j].second[i][k]].push_back(positionName);
                }
            }
        }
    }

    return prefixes;
}

std::string getSuffix(std::map<std::string, std::vector<std::string>> prefixes, Chain chain, std::vector<Atom*> principalGroupAtoms, unsigned int principalGroup, bool getPIN = true) {
    std::string name;
    std::string suffix;
    std::vector<std::pair<Atom*, unsigned int>> bonds;
    std::map<std::string, std::vector<std::string>> radicalNames;
    std::vector<std::string> orderedPrefixes;
    std::vector<Atom*> atoms;

    if (chain.carbox || principalGroup >= KETONE) {
        suffix += "-" + std::to_string(chain.principalGroups[chain.principalGroups.size()-1]+1);
        for (int i = (int)chain.principalGroups.size()-2; i >= 0; i--) {
            suffix += "," + std::to_string(chain.principalGroups[i]+1);
        }
        suffix += "-";
    }

    std::vector<std::pair<std::string, std::vector<std::vector<std::string>>>> aminePrefixes;
    switch (principalGroup) {
        case CHAIN:
            suffix = "o";
            break;

        case CARBOXYL:
            if (chain.carbox)
                suffix += greekPrefix[chain.principalGroups.size()] + "carboxílico";
            else
                suffix += greekPrefix[chain.principalGroups.size()] + "óico";

            name = "ácido " + name;
            break;
        
        case CARBOALKOXY:
            if (chain.carbox) {
                suffix += greekPrefix[chain.principalGroups.size()] + "carboxilato de ";
                for (int i = 0; i < chain.chain.size(); i++) {
                    for (int j = 0; j < principalGroupAtoms.size(); j++) {
                        if (principalGroupAtoms[j] == chain.chain[i]) {
                            for (int k = 0; k < chain.chain[i]->bonds.size(); k++) {
                                if (chain.chain[i]->bonds[k]->species == CARBON) {
                                    for (int l = 0; l < chain.chain.size(); l++) {
                                        if (chain.chain[i]->bonds[k] == chain.chain[l]) {
                                            goto loopkEnd;
                                        }
                                    }
                                    radicalNames = getOxyPrefix(radicalNames, chain.chain[i]->bonds[k], std::to_string(chain.chain.size()-i));
                                }
                                loopkEnd:;
                            }
                            break;
                        }
                    }
                }
            }
            else {
                suffix = greekPrefix[chain.principalGroups.size()] + "oato de ";

                for (int i = 0; i < chain.chain.size(); i++) {
                    for (int j = 0; j < principalGroupAtoms.size(); j++) {
                        if (principalGroupAtoms[j] == chain.chain[i]) {
                            radicalNames = getOxyPrefix(radicalNames, chain.chain[i], std::to_string(chain.chain.size()-i));
                            break;
                        }
                    }
                }
            }

            orderedPrefixes = orderPrefixes(radicalNames);
            if (orderedPrefixes.size() == 1) {
                if (orderedPrefixes[0][orderedPrefixes[0].size()-1] == ')') {
                    suffix += greekPrefix[chain.principalGroups.size()] + orderedPrefixes[0].substr(1, orderedPrefixes[0].size()-2) + "o";
                }
                else {
                    suffix += greekPrefix[chain.principalGroups.size()] + orderedPrefixes[0] + "o";    
                }
            }
            else {
                for (int i = 0; i < orderedPrefixes.size(); i++) {
                    suffix += radicalNames[orderedPrefixes[i]][0];
                    for (int j = 1; j < radicalNames[orderedPrefixes[i]].size(); j++) {
                        suffix += "," + radicalNames[orderedPrefixes[i]][j];
                    }
                    
                    if (orderedPrefixes[i][orderedPrefixes[i].size()-1] == ')') {
                        suffix += "-" + greekPrefix[radicalNames[orderedPrefixes[i]].size()] + orderedPrefixes[i].substr(1, orderedPrefixes[i].size()-2);
                    }
                    else {
                        suffix += "-" + greekPrefix[radicalNames[orderedPrefixes[i]].size()] + orderedPrefixes[i];
                    }

                    if (i == orderedPrefixes.size()-1) {
                        suffix += "o";
                    }
                    else {
                        suffix += " ";
                    }
                }
            }
            break;
        
        case HALOFORMYL:
            switch (chain.halogen) {
                case FLUORINE:
                    name = greekPrefix[chain.principalGroups.size()] + "fluoreto de " + name;
                    break;

                case CHLORINE:
                    name = greekPrefix[chain.principalGroups.size()] + "cloreto de " + name;
                    break;

                case BROMINE:
                    name = greekPrefix[chain.principalGroups.size()] + "brometo de " + name;
                    break;

                case IODINE:
                    name = greekPrefix[chain.principalGroups.size()] + "iodeto de " + name;
                    break;
                
                default:
                    name = greekPrefix[chain.principalGroups.size()] + "halogeneto de " + name;
                    break;
            }

            if (chain.carbox)
                suffix += greekPrefix[chain.principalGroups.size()] + "carbonilo";
            else
                suffix += greekPrefix[chain.principalGroups.size()] + "oílo";

            break;

        case CARBOXAMIDE:
            if (chain.carbox) {
                for (int i = 0; i < chain.chain.size(); i++) {
                    for (int j = 0; j < principalGroupAtoms.size(); j++) {
                        if (chain.chain[i] == principalGroupAtoms[j]) {
                            if (principalGroupAtoms.size() > 1 && chain.principalGroups.size() > 1) {
                                aminePrefixes.push_back({"N" + std::to_string(chain.chain.size()-i), {}});
                            }
                            else {
                                aminePrefixes.push_back({"N", {}});
                            }

                            for (int k = 0; k < chain.chain[i]->bonds.size(); k++) {
                                if (getPrincipalGroup({chain.chain[i]->bonds[k]}).principalGroup == CARBOXAMIDE) {
                                    std::vector<std::vector<std::string>> newAminePrefixes;
                                    newAminePrefixes = getAminePrefix(chain.chain[i]->bonds[k]);
                                    aminePrefixes[aminePrefixes.size()-1].second.insert(aminePrefixes[aminePrefixes.size()-1].second.begin(), newAminePrefixes.begin(), newAminePrefixes.end());
                                }
                            }
                            break;
                        }
                    }
                }
            
                suffix += greekPrefix[chain.principalGroups.size()] + "carboxamida";
            }
            else {
                for (int i = 0; i < chain.chain.size(); i++) {
                    for (int j = 0; j < principalGroupAtoms.size(); j++) {
                        if (chain.chain[i] == principalGroupAtoms[j]) {
                            std::vector<std::vector<std::string>> newAminePrefixes;
                            newAminePrefixes = getAminePrefix(chain.chain[i]);
                            if (principalGroupAtoms.size() > 1 && chain.principalGroups.size() > 1) {
                                aminePrefixes.push_back({"N" + std::to_string(chain.chain.size()-i), newAminePrefixes});
                            }
                            else {
                                aminePrefixes.push_back({"N", newAminePrefixes});
                            }
                            break;
                        }
                    }
                    loopiEndAmide:;
                }

                suffix += greekPrefix[chain.principalGroups.size()];
                if (suffix.size() > 0 && suffix[suffix.size()-1] == 'a') {
                    suffix += "mida";
                }
                else {
                    suffix += "amida";
                }
            }

            prefixes = addAminePrefixes(prefixes, aminePrefixes);
            break;
        
        case NITRILE:
            if (chain.carbox) {
                suffix += greekPrefix[chain.principalGroups.size()] + "carbonitrilo";
            }
            else {
                suffix += greekPrefix[chain.principalGroups.size()] + "nitrilo";
            }
            break;
        
        case ALDEHYDE:
            if (chain.carbox) {
                suffix += greekPrefix[chain.principalGroups.size()] + "carbaldeído";
            }
            else {
                suffix += greekPrefix[chain.principalGroups.size()] + "al";
            }
            break;

        case KETONE:
            suffix += greekPrefix[chain.principalGroups.size()] + "ona";
            break;
        
        case HYDROXYL:
            suffix += greekPrefix[chain.principalGroups.size()];
            if (suffix[suffix.size()-1] == 'a') {
                suffix = suffix.substr(0, suffix.size()-1);
            }
            suffix += "ol";
            break;
        
        case SULFHYDRYL:
            suffix += greekPrefix[chain.principalGroups.size()] + "tiol";
            break;
        
        case AMINE:
            for (int i = 0; i < chain.chain.size(); i++) {
                for (int j = 0; j < principalGroupAtoms.size(); j++) {
                    if (chain.chain[i] == principalGroupAtoms[j]) {
                        std::vector<std::vector<std::string>> newAminePrefixes;
                        newAminePrefixes = getAminePrefix(chain.chain[i]);
                        if (principalGroupAtoms.size() > 1 && chain.principalGroups.size() > 1) {
                            aminePrefixes.push_back({"N" + std::to_string(chain.chain.size()-i), newAminePrefixes});
                        }
                        else {
                            aminePrefixes.push_back({"N", newAminePrefixes});
                        }
                        break;
                    }
                }
            }
            prefixes = addAminePrefixes(prefixes, aminePrefixes);

            suffix += greekPrefix[chain.principalGroups.size()];
            if (suffix.size() > 0 && suffix[suffix.size()-1] == 'a') {
                suffix += "mina";
            }
            else {
                suffix += "amina";
            }
            break;
    }

    bool isBenzene = chain.multipleBonds.size() == 3 && chain.doubleBonds.size() == 3 && chain.doubleBonds[0]-chain.doubleBonds[1] == 2 && chain.doubleBonds[1] - chain.doubleBonds[2] == 2;

    if ((chain.isCyclic || chain.chain.size() == 2) && chain.substituents.size() == 1 && (chain.multipleBonds.size() == 0 || isBenzene) && principalGroup >= KETONE || chain.chain.size() == 1 && (principalGroup >= KETONE || chain.carbox) || principalGroup != CHAIN && isBenzene && chain.principalGroups.size() == 1) {
        suffix = suffix.substr(1 + 2*chain.principalGroups.size(), suffix.size());
    }
    
    std::string vowels = "aeiouó";
    for (int i = 0; i < suffix.size(); i++) {
        if (std::isalpha(suffix[i])) {
            for (int j = 0; j < vowels.size(); j++) {
                if (suffix[i] == vowels[j]) {
                    break;
                }
                if (j == vowels.size()-1) {
                    suffix = "o" + suffix;
                }
            }
            break;
        }
    }

    orderedPrefixes = orderPrefixes(prefixes);

    for (int i = 0; i < orderedPrefixes.size(); i++) {
        std::string prefixName = "";
        if (principalGroup == AMINE || principalGroup == CARBOXAMIDE || chain.chain.size() > 1 && !(chain.chain.size() == 2 && principalGroup != CHAIN && principalGroup < KETONE)) {
            prefixName += prefixes[orderedPrefixes[i]][0];
            for (int j = 1; j < prefixes[orderedPrefixes[i]].size(); j++) {
                prefixName += "," + prefixes[orderedPrefixes[i]][j];
            }
            prefixName += "-";
        }
        if (orderedPrefixes[i][0] == 'h' && prefixes[orderedPrefixes[i]].size() > 1) {
            orderedPrefixes[i] = prefixName + greekPrefix[prefixes[orderedPrefixes[i]].size()] + "-" + orderedPrefixes[i];
        }
        else {
            orderedPrefixes[i] = prefixName + greekPrefix[prefixes[orderedPrefixes[i]].size()] + orderedPrefixes[i];
        }
    }

    if (orderedPrefixes.size() > 0) {
        name += orderedPrefixes[0];
        for (int i = 1; i < orderedPrefixes.size(); i++) {
            name += "-" + orderedPrefixes[i];
        }
    }

    if (principalGroup == CHAIN && (chain.isCyclic && chain.multipleBonds.size() == 0 || isBenzene || chain.chain.size() == 2) && prefixes.size() == 1 && prefixes.begin()->second.size() == 1) {
        name = name.substr(2, name.size());
    }

    std::string carbonChain = nameCarbonChain(chain);
    if (carbonChain[0] == 'h' && orderedPrefixes.size() > 0 && name[name.size()-1] != ')') {
        carbonChain = '-' + carbonChain;
    }

    // benzene
    if (carbonChain == "ciclo-hexa-1,3,5-trien" || carbonChain == "ciclo-hexa-2,4,6-trien") {
        carbonChain = "benzen";
    }
    if (getPIN && carbonChain == "etan" && principalGroup < KETONE && chain.principalGroups.size() == 1) {
        if (principalGroup == CARBOXYL) {
            carbonChain = "acético";
            suffix = "";
        }
        else if (principalGroup == CARBOALKOXY) {
            carbonChain = "acetato";
            suffix = suffix.substr(std::string("oato").size(), suffix.size()-std::string("oato").size());
        }
        else {
            carbonChain = "acet";
            if (principalGroup == ALDEHYDE) {
                suffix = "aldeído";
            }
        }
    }

    name += carbonChain;
    name = name + suffix;

    if (name.size() >= std::string("benzen").size() + suffix.size() && name.substr(name.size() - std::string("benzen").size() - suffix.size(), std::string("benzenocarboxilato").size()) == "benzenocarboxilato") {
        name = name.substr(0, name.size() - std::string("benzen").size() - suffix.size()) + "benzoato" + name.substr(name.size() - suffix.size() + std::string("ocarboxilato").size(), suffix.size() - std::string("ocarboxilato").size());
    }

    return name;
}

std::string Molecule::getName(bool getPIN) {
    // find functional groups
    principalGroupInfo info = getPrincipalGroup(atoms);
    std::vector<Atom*> principalGroupAtoms = info.principalGroupAtoms;
    unsigned int principalGroup = info.principalGroup;

    // get possible chains with principal group
    std::vector<Chain> possibleChains;
    for (int i = 0; i < principalGroupAtoms.size(); i++) {
        Chain chain = principalGroupAtoms[i]->getPrincipalChain(std::vector<Atom*>(), principalGroupAtoms, principalGroup);
        chain = getOrderedChain(chain, principalGroupAtoms, principalGroup);
        possibleChains.push_back(chain);
    }
    possibleChains = compareBranches(possibleChains);

    // carbox-
    std::vector<Atom*> ignoredAtoms;
    std::vector<Atom*> carboxGroups;
    if (principalGroup < KETONE && principalGroup != CHAIN && possibleChains[0].chain.size() != 2) {
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

        for (int i = 0; i < principalGroupAtoms.size(); i++) {
            for (int j = 0; j < principalGroupAtoms[i]->bonds.size(); j++) {
                if (principalGroupAtoms[i]->bonds[j]->species == CARBON) {
                    carboxGroups.push_back(principalGroupAtoms[i]->bonds[j]);
                    ignoredAtoms.push_back(principalGroupAtoms[i]);
                    break;
                }
            }
        }
        
        std::vector<Chain> carboxChains;
        for (int i = 0; i < carboxGroups.size(); i++) {
            Chain chain = carboxGroups[i]->getPrincipalChain(std::vector<Atom*>(), carboxGroups, principalGroup, ignoredAtoms);
            chain = getOrderedChain(chain, carboxGroups, principalGroup, ignoredAtoms);
            if (possibleChains.size() == 0 || compareBranchesPropertySize<unsigned int>(std::vector<Chain>{possibleChains[0], chain}, [](Chain c) {return c.principalGroups;}).size() == 1) {
                chain.carbox = true;
                carboxChains.push_back(chain);
            }
        }
        if (carboxChains.size() > 0) {
            principalGroupAtoms = carboxGroups;
            possibleChains = carboxChains;
        }
        else {
            ignoredAtoms.clear();
        }
    }

    // possible chains without principle group
    if (possibleChains.size() == 0) {
        for (int i = 0; i < atoms.size(); i++) {
            if (atoms[i]->species == CARBON) {
                Chain chain = atoms[i]->getPrincipalChain(std::vector<Atom*>(), principalGroupAtoms, principalGroup);
                chain = getOrderedChain(chain, principalGroupAtoms, principalGroup);
                possibleChains.push_back(chain);
            }
        }
    }
    
    // determine principal chain
    principalChain = compareChains(possibleChains, principalGroupAtoms, principalGroup)[0];
    
    if (principalGroup == HALOFORMYL) {
        principalChain.halogen = info.halogen;
    }

    for (int i = 0; i < principalGroupAtoms.size(); i++) {
        for (int j = 0; j < principalChain.chain.size(); j++) {
            if (principalGroupAtoms[i] == principalChain.chain[j]) {
                break;
            }
            if (j == principalChain.chain.size()-1) {
                principalGroupAtoms.erase(principalGroupAtoms.begin() + i);
                i--;
            }
        }
    }
    for (int i = 1; i < principalGroupAtoms.size(); i++) {
        for (int j = 0; j < i; j++) {
            if (principalGroupAtoms[i] == principalGroupAtoms[j]) {
                principalGroupAtoms.erase(principalGroupAtoms.begin() + i);
                i--;
                break;
            }
        }
    }

    // name substituents
    std::map<std::string, std::vector<std::string>> prefixes = getChainPrefixes(principalChain, principalGroupAtoms, principalGroup, ignoredAtoms);

    // get suffix
    std::string name = getSuffix(prefixes, principalChain, principalGroupAtoms, principalGroup, getPIN);
    

    if (name.size() >= std::string("ácido benzenocarboxílico").size() && name.substr(name.size() - std::string("benzenocarboxílico").size(), std::string("benzenocarboxílico").size()) == "benzenocarboxílico") {
        name = name.substr(0, name.size()-std::string("benzenocarboxílico").size()) + "benzóico";
    }
    
    else if (name.size() >= std::string("benzenoílo").size() && name.substr(name.size() - std::string("benzenoílo").size(), std::string("benzenoílo").size()) == "benzenoílo") {
        name = name.substr(0, name.size()-std::string("benzenoílo").size()) + "benzoílo";
    }

    else if (name.size() >= std::string("benzenocarboxamida").size() && name.substr(name.size() - std::string("benzenocarboxamida").size(), std::string("benzenocarboxamida").size()) == "benzenocarboxamida") {
        name = name.substr(0, name.size()-std::string("benzenocarboxamida").size()) + "benzamida";
    }
    
    else if (name.size() >= std::string("benzenocarbaldeído").size() && name.substr(name.size() - std::string("benzenocarbaldeído").size(), std::string("benzenocarbaldeído").size()) == "benzenocarbaldeído") {
        name = name.substr(0, name.size()-std::string("benzenocarbaldeído").size()) + "benzaldeído";
    }

    else if (name.size() >= std::string("benzenol").size() && name.substr(name.size() - std::string("benzenol").size(), std::string("benzenol").size()) == "benzenol") {
        name = name.substr(0, name.size()-std::string("benzenol").size()) + "fenol";
    }

    if (getPIN) {
        if (name == "metoxibenzeno") {
            name = "anisol";
        }
    }


    return name;
}


template <typename T> std::vector<Chain> compareBranchesPropertyPosition(std::vector<Chain> chains, std::vector<T> (*property)(Chain)) {
    for (int i = 1; i < chains.size(); i++) {
        int j = 0;
        while (true) {
            if (j == property(chains[0]).size()) {
                if (chains[0].branches.size() == 0)
                    break;
                
                for (int k = 1; k < chains[0].branches.size(); k++)
                    chains.push_back(addBranch(chains[0], chains[0].branches[k]));
                chains[0] = addBranch(chains[0], chains[0].branches[0]);
                continue;
            }
            if (j == property(chains[i]).size()) {
                if (chains[i].branches.size() == 0)
                    break;

                for (int k = 1; k < chains[i].branches.size(); k++)
                    chains.push_back(addBranch(chains[i], chains[i].branches[k]));
                chains[i] = addBranch(chains[i], chains[i].branches[0]);
                continue;
            }

            if (property(chains[0])[property(chains[0]).size()-j-1] < property(chains[i])[property(chains[i]).size()-j-1]) {
                chains.erase(chains.begin() + i);
                i--;
            }
            else if (property(chains[0])[property(chains[0]).size()-j-1] > property(chains[i])[property(chains[i]).size()-j-1]) {
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

std::vector<Chain> comparePrefixes(std::vector<Chain> chains, std::vector<Atom*> principalGroupAtoms, unsigned int principalGroup, std::vector<Atom*> ignoredAtoms) {
    std::vector<int> index;
    for (int i = 0; i < chains.size(); i++) {
        index.push_back(i);
    }
    
    std::vector<std::map<std::string, std::vector<std::string>>> prefixes;
    for (int i = 0; i < chains.size(); i++) {
        prefixes.push_back(getChainPrefixes(chains[i], principalGroupAtoms, principalGroup, ignoredAtoms));
    }
    std::vector<std::string> prefixOrder = orderPrefixes(prefixes[0]);

    std::string prefix;
    for (int i = 0; i < prefixOrder.size(); i++) {
        prefix = prefixOrder[i];
        for (int j = 1; j < index.size(); j++) {
            if (prefixes[index[0]][prefix].size() > prefixes[index[j]][prefix].size()) {
                index.erase(index.begin() + j);
                j--;
                continue;
            }
            if (prefixes[index[0]][prefix].size() < prefixes[index[j]][prefix].size()) {
                while (j > 0) {
                    index.erase(index.begin());
                    j--;
                }
                continue;
            }

            for (int k = (int)prefixes[index[0]][prefix].size()-1; k >= 0; k--) {
                if (prefixes[index[0]][prefix][k] < prefixes[index[j]][prefix][k]) {
                    index.erase(index.begin() + j);
                    j--;
                    break;
                }
                else {
                    while (j > 0) {
                        index.erase(index.begin());
                        j--;
                    }
                    break;
                }
            }
        }
    }

    std::vector<Chain> bestChains;
    for (int i = 0; i < index.size(); i++) {
        bestChains.push_back(chains[index[i]]);
    }
    return bestChains;
}

std::vector<Chain> compareChains(std::vector<Chain> chains, std::vector<Atom*> principalGroupAtoms, unsigned int principalGroup, std::vector<Atom*> ignoredAtoms) {
    chains = compareBranches(chains);
    
    // lower numbers for suffixes
    chains = compareBranchesPropertyPosition<unsigned int>(chains, [](Chain c) {return c.principalGroups;});

    // not benzene
    if (!(chains[0].isCyclic && chains[0].doubleBonds.size() == 3 && chains[0].chain.size() == 6 && chains[0].doubleBonds[0]-chains[0].doubleBonds[1] == 2 && chains[0].doubleBonds[1]-chains[0].doubleBonds[2] == 2)) {
        // lower numbers for multiple bonds
        chains = compareBranchesPropertyPosition<unsigned int>(chains, [](Chain c) {return c.multipleBonds;});

        // lower numbers for double bonds
        chains = compareBranchesPropertyPosition<unsigned int>(chains, [](Chain c) {return c.doubleBonds;});
    }

    // greatest number of prefixes
    chains = compareBranchesPropertySize<unsigned int>(chains, [](Chain c) {return c.substituents;});

    // lower numbers for prefixes
    chains = compareBranchesPropertyPosition<unsigned int>(chains, [](Chain c) {return c.substituents;});

    // lower numbers for most substituted amine
    if (principalGroup == AMINE || principalGroup == CARBOXAMIDE) {
        int position;
        std::vector<int> bondsSize0;
        std::vector<int> bondsSizei;
        for (int i = 1; i < chains.size(); i++) {
            for (int j = (int)chains[0].principalGroups.size()-1; j >= 0; j--) {
                bondsSize0.clear();
                bondsSizei.clear();
                position = chains[0].chain.size() - chains[0].principalGroups[j] - 1;
                
                if (chains[0].carbox) {
                    for (int k = 0; k < chains[0].chain[position]->bonds.size(); k++) {
                        if (chains[0].chain[position]->bonds[k]->species != CARBON) {
                            continue;
                        }
                        for (int l = 0; l < chains[0].chain.size(); l++) {
                            if (chains[0].chain[l] == chains[0].chain[position]->bonds[k]) {
                                goto skipBonds0;
                            }
                        }
                        for (int l = 0; l < chains[0].chain[position]->bonds[k]->bonds.size(); l++) {
                            if (chains[0].chain[position]->bonds[k]->bonds[l]->species == NITROGEN) {
                                bondsSize0.push_back(chains[0].chain[position]->bonds[k]->bonds[l]->bonds.size());
                            }
                        }
                        skipBonds0:;

                        if (chains[i].chain[position]->bonds[k]->species != CARBON) {
                            continue;
                        }
                        for (int l = 0; l < chains[i].chain.size(); l++) {
                            if (chains[i].chain[l] == chains[i].chain[position]->bonds[k]) {
                                goto skipBondsi;
                            }
                        }
                        for (int l = 0; l < chains[i].chain[position]->bonds[k]->bonds.size(); l++) {
                            if (chains[i].chain[position]->bonds[k]->bonds[l]->species == NITROGEN) {
                                bondsSizei.push_back(chains[i].chain[position]->bonds[k]->bonds[l]->bonds.size());
                            }
                        }
                        skipBondsi:;
                    }
                    while (j > 0 && chains[0].principalGroups[j-1] == chains[0].principalGroups[j]) {
                        j--;
                    }
                }
                else {
                    for (int k = 0; k < chains[0].chain[position]->bonds.size(); k++) {
                        if (chains[0].chain[position]->bonds[k]->species == NITROGEN) {
                            bondsSize0.push_back(chains[0].chain[position]->bonds[k]->bonds.size());
                        }
                        if (chains[i].chain[position]->bonds[k]->species == NITROGEN) {
                            bondsSizei.push_back(chains[i].chain[position]->bonds[k]->bonds.size());
                        }
                    }
                    while (j > 0 && chains[0].principalGroups[j-1] == chains[0].principalGroups[j]) {
                        j--;
                    }
                }

                // cut common values
                for (int k = 0; k < bondsSize0.size(); k++) {
                    for (int l = 0; l < bondsSizei.size(); l++) {
                        if (bondsSize0[k] == bondsSizei[l]) {
                            bondsSize0.erase(bondsSize0.begin() + k);
                            bondsSizei.erase(bondsSizei.begin() + l);
                            k--;
                            l--;
                            break;
                        }
                    }
                }
                // compare maximum
                if (bondsSize0.size() > 0) {
                    int max = bondsSize0[0];
                    for (int k = 1; k < bondsSize0.size(); k++) {
                        if (bondsSize0[k] > max) {
                            max = bondsSize0[k];
                        }
                    }
                    for (int k = 0; k < bondsSizei.size(); k++) {
                        if (bondsSizei[k] > max) {
                            while (i > 0) {
                                chains.erase(chains.begin());
                                i--;
                                goto loopjEnd;
                            }
                        }
                    }
                    chains.erase(chains.begin() + i);
                    i--;
                }
                loopjEnd:;
            }
        }
    }

    // lower numbers for first prefixes
    if (chains.size() > 1 && chains[0].substituents.size() > 0) {
        chains = comparePrefixes(chains, principalGroupAtoms, principalGroup, ignoredAtoms);
    }

    return chains;
}

Chain getOrderedChain(Chain chain, std::vector<Atom*> principalGroupAtoms, unsigned int principalGroup, std::vector<Atom*> ignoredAtoms) {
    std::vector<Chain> principalBranches = {};
    std::vector<Chain> principalChains = {};

    if (chain.branches.size() == 0) {
        principalChains.push_back(chain);
        principalChains.push_back(invertChain(chain));
    }

    else {
        if (chain.chain.size() > 1) {
            // forward
            principalBranches = chain.branches;
            for (int i = 0; i < principalBranches.size(); i++) {
                principalBranches[i] = addBranch(chain, principalBranches[i]);
            }
            principalBranches = compareChains(principalBranches, principalGroupAtoms, principalGroup, ignoredAtoms);
            principalChains.push_back(principalBranches[0]);

            // backward
            principalBranches = chain.branches;
            principalBranches = compareBranches(principalBranches);
            for (int i = 0; i < principalBranches.size(); i++) {
                principalBranches[i] = addBranch(chain, principalBranches[i]);
                principalBranches[i] = invertChain(principalBranches[i]);
            }
            principalBranches = compareChains(principalBranches, principalGroupAtoms, principalGroup, ignoredAtoms);
            principalChains.push_back(principalBranches[0]);
        }
        else {
            std::vector<Chain> firstBranch;
            std::vector<Chain> secondBranch;
            
            for (int i = 0; i < chain.branches.size(); i++) {
                firstBranch = {addBranch(chain, chain.branches[i])};
                firstBranch = compareBranches(firstBranch);
                for (int j = 0; j < firstBranch.size(); j++) {
                    firstBranch[j] = invertChain(firstBranch[j]);
                }
                firstBranch = compareChains(firstBranch, principalGroupAtoms, principalGroup, ignoredAtoms);

                
                secondBranch.clear();
                for (int j = 0; j < chain.branches.size(); j++) {
                    if (chain.branches[i].chain[chain.branches[i].chain.size()-1] == chain.branches[j].chain[chain.branches[j].chain.size()-1]) {
                        continue;
                    }
                    secondBranch.push_back(joinChains(firstBranch[0], chain.branches[j]));
                }

                secondBranch = compareChains(secondBranch, principalGroupAtoms, principalGroup, ignoredAtoms);
                principalChains.push_back(secondBranch[0]);
            }
        }
    }

    for (int i = 0; i < principalChains.size(); i++) {
        principalChains[i].branches.clear();
    }
    principalChains = compareChains(principalChains, principalGroupAtoms, principalGroup, ignoredAtoms);
    return principalChains[0];
}
