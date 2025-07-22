#include <iostream>
#include <vector>
#include <string>
#include <map>

#include "molecule.h"

Molecule* readSmiles(std::string smiles) {
    std::vector<Atom*> atoms;
    unsigned int species;
    std::vector<Atom*> previousAtoms;
    std::map<char, Atom*> labledAtoms;
    bool inNewBranch = true;
    int bondOrder = 1;

    for (int i = 0; i < smiles.size(); i++) {
        if (std::isalpha(smiles[i])) {
            switch (smiles[i]) {
                case 'C':
                    if (smiles.size() > i + 1) {
                        if (smiles[i+1] == 'l') {
                            i++;
                            species = CHLORINE;
                            break;
                        }
                    }
                    species = CARBON;
                    break;
                
                case 'N':
                    species = NITROGEN;
                    break;
                
                case 'O':
                    species = OXYGEN;
                    break;

                case 'F':
                    species = FLUORINE;
                    break;
                
                case 'S':
                    species = SULFUR;
                    break;
                
                case 'B':
                    if (smiles.size() > i + 1) {
                        if (smiles[i+1] == 'r') {
                            i++;
                            species = BROMINE;
                            break;
                        }
                    }
                
                default:
                    std::cout << "Error reading smiles string: unrecognized element \"" << smiles[i] << "\"" << std::endl;
                    break;
            }

            if (atoms.size() == 0) {
                atoms.push_back(new Atom(species, nullptr));
            }
            else {
                atoms.push_back(new Atom(species, previousAtoms[previousAtoms.size()-1]));
                
                for (int j = 1; j < bondOrder; j++) {
                    atoms[atoms.size()-1]->bonds.push_back(previousAtoms[previousAtoms.size()-1]);
                    previousAtoms[previousAtoms.size()-1]->bonds.push_back(atoms[atoms.size()-1]);
                }
                bondOrder = 1;
            }

            if (inNewBranch) {
                previousAtoms.push_back(atoms[atoms.size()-1]);
                inNewBranch = false;
            }
            else {
                previousAtoms[previousAtoms.size()-1] = atoms[atoms.size()-1];
            }
        }

        else if (std::isdigit(smiles[i])) {
            if (labledAtoms.find(smiles[i]) == labledAtoms.end()) {
                labledAtoms[smiles[i]] = atoms[atoms.size()-1];
            }
            else {
                atoms[atoms.size()-1]->bonds.push_back(labledAtoms[smiles[i]]);
                labledAtoms[smiles[i]]->bonds.push_back(atoms[atoms.size()-1]);
                labledAtoms.erase(smiles[i]);
            }
        }

        else {
            switch (smiles[i]) {
                case '(':
                    inNewBranch = true;
                    break;
                
                case ')':
                    previousAtoms.erase(previousAtoms.end()-1);
                    break;
                
                case '=':
                    bondOrder = 2;
                    break;
                
                case '#':
                    bondOrder = 3;
                    break;

                default:
                    std::cout << "Error reading smiles string: unrecognized symbol \"" << smiles[i] << "\"" << std::endl;
                    break;
            }
        }
    }
    
    return new Molecule(atoms);
}

int validateName(std::string name, std::string smiles) {
    std::cout << "Testing \"" << name << "\"..." << std::endl;
    Molecule* molecule = readSmiles(smiles);
    
    std::string generatedName = molecule->getName();
    if (generatedName == name) {
        return 0;
    }
    else {
        std::cout << "Error naming \"" << name << "\" : getName returned \"" << generatedName << "\"" << std::endl;;
        return 1;
    }
}

int main() {
    std::vector<std::pair<std::string, std::string>> tests = {
        // IUPAC Blue Book
        // P-12.1
        {"metano", "C"},
        {"pentano", "CCCCC"},
        {"ciclo-hexano", "C1CCCCC1"},
        {"etenilbenzeno", "C=CC1=CC=CC=C1"},
        {"fenol", "C1=CC=C(C=C1)O"},
        {"ácido acético", "CC(=O)O"},
        {"propan-2-ona", "CC(=O)C"},
        {"5-cloropentan-2-ona", "CC(=O)CCCCl"},
        {"1-etoxipropano", "CCOCCC"},
        {"4-cloropentan-2-ona", "CC(Cl)CC(=O)C"},
        {"ciclo-hexeno", "C1=CCCCC1"},
        // P-13.1
        {"etano", "CC"},
        {"etanotiol", "CCS"},
        {"benzeno", "C1=CC=CC=C1"},
        {"bromobenzeno", "C1=CC=CC=C1Br"},
        // P-13.2.1.1
        {"ciclotetradecano", "C1CCCCCCCCCCCCC1"},
        // P-13.2.2.2
        {"ácido benzóico", "C1=CC=CC=C1C(=O)O"},
        {"ácido hexanóico", "CCCCCC(=O)O"},
        {"ácido 4-formilbenzóico", "C1=CC(C=O)=CC=C1C(=O)O"},
        // P-13.3.3.1
        {"acetonitrilo", "CC#N"}
    };

    for (int i = 0; i < tests.size(); i++) {
        if (validateName(tests[i].first, tests[i].second) == 1) {
            //return 1;
        }
    }

    std::cout << "Testing complete!" << std::endl;
    return 0;
}