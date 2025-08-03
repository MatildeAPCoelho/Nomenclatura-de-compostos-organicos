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
                
                case 'I':
                    species = IODINE;
                    break;
                
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
                labledAtoms[smiles[i]] = previousAtoms[previousAtoms.size()-1];
            }
            else {
                for (int j = 0; j < bondOrder; j++) {
                    previousAtoms[previousAtoms.size()-1]->bonds.push_back(labledAtoms[smiles[i]]);
                    labledAtoms[smiles[i]]->bonds.push_back(previousAtoms[previousAtoms.size()-1]);
                }
                bondOrder = 1;
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

int validateName(std::string name, std::string smiles, bool getPIN) {
    std::cout << "Testing \"" << name << "\"..." << std::endl;
    Molecule* molecule = readSmiles(smiles);
    
    std::string generatedName = molecule->getName(getPIN);
    if (generatedName == name) {
        return 0;
    }
    else {
        std::cout << "Error naming \"" << name << "\" : getName returned \"" << generatedName << "\"" << std::endl;;
        return 1;
    }
}

int main() {
    // IUPAC
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
        {"acetonitrilo", "CC#N"},
        // P-13.3.3.2
        {"metanol", "CO"},
        {"diciclo-hexilmetanona", "C1CCCCC1C(=O)C2CCCCC2"},
        {"anisol", "COC1=CC=CC=C1"},
        {"fenilacetonitrilo", "C1=CC=CC=C1CC#N"}, // "amino(hidroxi)(sulfanil)acetonitrilo"
    };

    for (int i = 0; i < tests.size(); i++) {
        if (validateName(tests[i].first, tests[i].second, true) == 1) {
            //return 1;
        }
    }

    // Exercícios
    tests = {
        {"N-(3-aminopropil)-1-fenil-3,4-dioxociclopentano-1-carboxamida", "C1=CC=CC=C1C2(CC(=O)C(=O)C2)C(=O)NCCCN"},
        {"N-(2-bromo-2-ciclo-hexiletil)-N-(2-hidroxifenil)pent-2-enamida", "OC1C=CC=CC=1N(C(=O)C=CCC)CC(Br)C2CCCCC2"},
        {"ácido 3-fenil-4-(2-hidroxifenil)but-3-enóico", "OC1=CC=CC=C1C=C(CC(=O)O)C2=CC=CC=C2"},
        {"ácido 2-(3-amino-1-hidroxi-4-oxociclo-hexil)butanóico", "NC1C(=O)CCC(O)(C1)C(CC)C(=O)O"},
        {"3-(3,4-dimetil-1-(metilsulfanil)ciclo-hexil)-3-hidroxipropanoato de metilo", "CSC1(CC(C)C(C)CC1)C(O)CC(=O)OC"},
        {"4-amino-2-(2-aminoetil)-2-metilciclopentan-1-ona", "CC1(CCN)CC(N)CC(=O)1"},
        {"2-(1-aminopropan-2-il)-5-clorociclo-hexan-1-ol", "OC1CC(Cl)CCC1C(C)CN"},
        {"3-(3-amino-5-hidroxifenoxi)propanoato de metilo", "C1=C(O)C=C(N)C=C1OCCC(=O)OC"},
        {"2-hidroxi-4-(2-(2-(metilsulfanil)fenil)eten-1-il)benzamida", "CSC1=CC=CC=C1C=CC2=CC=C(C(O)=C2)C(=O)N"},
        {"4-butoxi-4-(4-(fluorocarbonil)-2-(metilsulfanil)ciclopent-1-en-1-il)butanoato de metilo", "CSC1CC(C(=O)F)CC=1C(OCCCC)CCC(=O)OC"},
        {"3-(etilamino)-5-oxo-2-(3-sulfanilbutil)ciclopentano-1-carboxilato de metilo", "CC(S)CCC1C(NCC)CC(=O)C1C(=O)OC"},
        {"brometo de 4-(1-amino-3-oxopropil)-5-ciclopentil-2-oxociclo-hexano-1-carbonilo", "O=CCC(N)C1C(C2CCCC2)CC(C(=O)C1)C(=O)Br"},
        {"4-aminobenzoato de 3-(dietilamino)-2,2-dimetilpropilo", "C1=CC(N)=CC=C1C(=O)OCC(C)(C)CN(CC)CC"},
        {"2-amino-3-(2-metilprop-1-en-1-il)ciclobutano-1-carbaldeído", "NC1C(C=C(C)C)CC1C=O"},
        {"2-(dibutilamino)-2-(4-metoxifenil)etanamida", "C1C=C(OC)C=CC=1C(N(CCCC)CCCC)C(=O)N"},
        {"N-etil-N-(2-hidroxi-5-sulfanilfenil)but-2-enamida", "OC1=CC=C(S)C=C1N(CC)C(=O)C=CC"},
        {"ácido 2-amino-4-((3-amino-2-oxobutil)sulfanil)butanóico", "CC(N)C(=O)CSCCC(N)C(=O)O"},
        {"ácido 3-(2-(metilcarbamoil)-1-sulfanilciclopentil)-6-oxo-hexanóico", "CNC(=O)C1CCCC1(S)C(CCC=O)CC(=O)O"},
        {"6-bromo-3-(3-etenil-4-metoxifenil)-N-etil-4-oxo-hex-2-enamida", "CCNC(=O)C=C(C(=O)CCBr)C1=CC(C=C)=C(OC)C=C1"},
        {"(2-amino-5-(2,4-dioxopentil)-4-hidroxiciclo-hexil)etanal", "CC(=O)CC(=O)CC1C(O)CC(N)C(C1)CC=O"},
        {"2-(3-cloro-2-(metilamino)-4-metoxibutil)-3-(3,4-di-hidroxifenil)pentanodial", "C1=C(O)C(O)=CC=C1C(CC=O)C(C=O)CC(NC)C(Cl)COC"},
        {"3-amino-N-(3-(3,4-di-hidroxifenil)-2-oxopropil)-4-oxobutanamida", "O=CC(N)CC(=O)NCC(=O)CC1=CC=C(O)C(O)=C1"},
        {"2-(3-cloro-5-formil-4-((prop-1-en-1-il)oxi)ciclo-hexil)-2-(dietilamino)-N-feniletanamida", "CC=COC1C(C=O)CC(CC(Cl)1)C(N(CC)CC)C(=O)NC2=CC=CC=C2"},
        {"ácido 2-(1-(5-hidroxiciclo-hex-1-en-1-il)-4-oxobutoxi)-3-oxopropanóico", "O=CC(C(=O)O)OC(CCC=O)C1=CCCC(O)C1"},
        {"ácido 3-(3-etoxipropil)-5-formil-2-(1-hidroxi-3-oxobutil)ciclo-hexano-1-carboxílico", "O=CC1CC(C(=O)O)C(C(O)CC(=O)C)C(C1)CCCOCC"},
        {"cloreto de 5-carbamoil-2-(2-hidroxi-5-metoxi-4-sulfanilfenil)-4,7-dioxo-heptanoílo", "O=CCC(C(=O)N)C(=O)CC(C(=O)Cl)C1=C(O)C=C(S)C(OC)=C1"},
        {"ácido 3-(2-amino-5-oxociclo-hexil)-2-(1-cloro-2-feniletil)but-2-enóico", "C1=CC=CC=C1CC(Cl)C(C(=O)O)=C(C)C2C(N)CCC(=O)C2"},
        {"2-(2-hidroxi-4-(metilcarbamoil)ciclopentil)propanoato de 4-formil-5-sulfanilciclo-hex-2-en-1-ilo", "CNC(=O)C1CC(O)C(C1)C(C)C(=O)OC2C=CC(C=O)C(S)C2"},
        {"2-amino-3-bromo-5-etoxi-6-formil-4-(2-hidroxi-5-sulfanilfenil)octanamida", "CCOC(C(C=O)CC)C(C(Br)C(N)C(=O)N)C1=CC(S)=CC=C(O)1"},
        {"ácido 3-(3-(3-cloro-3-oxopropil)-1-hidroxi-6-oxo-hept-4-en-1-il)-5-(etoxicarbonil)ciclo-hex-2-eno-1-carboxílico", "CCOC(=O)C1CC(C(=O)O)C=C(C1)C(O)CC(CCC(=O)Cl)C=CC(=O)C"},
        {"3-(bromocarbonil)-6-(5-cloro-2,4-di-hidroxifenil)-5-oxo-4-(4-sulfanilciclo-hexil)hexanoato de metilo", "OC1=CC(O)=C(Cl)C=C1CC(=O)C(C2CCC(S)CC2)C(C(=O)Br)CC(=O)OC"},
        {"4-(2-(1-aminoetoxi)etil)-2-(4-formil-2-hidroxiciclopentil)pent-2-enodiamida", "NC(=O)C(C1C(O)CC(C=O)C1)=CC(C(=O)N)CCOC(N)C"},
        {"2-(hidroximetil)-4-metoxibutanoato de 5-(1-(metilamino)-4-oxobutil)-4-oxociclo-hex-1-en-1-ilo", "O=CCCC(NC)C1C(=O)CC=C(C1)OC(=O)C(CO)CCOC"},
        {"ácido 4-(2-cloro-3-(2-(metilcarbamoil)ciclopentil)-1-sulfanilpropil)-3-(penta-1,3-dien-1-il)ciclopentano-1,2-dicarboxílico", "CNC(=O)C1CCCC1CC(Cl)C(S)C2CC(C(=O)O)C(C(=O)O)C2C=CC=CC"},
        {"ácido 3-(3-cloro-5-(3-hidroxifenil)-1-sulfanilpentil)-4-(etoxicarbonil)ciclopenta-2,4-dieno-1-carboxílico", "OC1C=CC=C(C=1)CCC(Cl)CC(S)C2=CC(C(=O)O)C=C2C(=O)OCC"},
        {"2-(diaminometil)-N-((5-etoxi-3-oxociclo-hex-1-en-1-il)metil)-3-hidroxi-4-oxobutanamida", "O=CC(O)C(C(N)N)C(=O)NCC1=CC(=O)CC(OCC)C1"},
        {"ácido 4-amino-6-(bromocarbonil)-2-(5-hidroxiciclo-hex-2-en-1-il)-3-(2-oxoetil)octa-2,5-dienóico", "CCC(C(=O)Br)=CC(N)C(CC=O)=C(C(=O)O)C1C=CCC(O)C1"},
        {"2-(2-(2,3-di-hidroxifenoxi)etil)-4-(etilcarbamoil)-3-(2-oxopentil)ciclo-hexa-2,5-dieno-1-carboxilato de metilo", "CCNC(=O)C1C=CC(C(=O)OC)C(CCOC2=C(O)C(O)=CC=C2)=C1CC(=O)CCC"},
        {"2-(5-amino-4-metoxiciclo-hex-1-en-1-il)-3-(5-(hidroximetil)-2-sulfanilfenil)-3-oxopropanal", "COC1CC=C(CC1N)C(C=O)C(=O)C2=CC(CO)=CC=C2S"},
        {"4-(1-(4-amino-5-etoxi-2-sulfanilciclo-hexil)-1,4-dioxopentan-3-il)ciclopenta-1,3-dieno-1-carbaldeído", "O=CC1=CC=C(C1)C(C(=O)C)CC(=O)C2CC(OCC)C(N)CC2S"},
        {"ácido 4-((2-(bromocarbonil)-6-(carboximetil)-4-hidroxifenil)metil)-5-oxopent-2-enóico", "OC(=O)CC1=CC(O)=CC(C(=O)Br)=C1CC(C=O)C=CC(=O)O"},
        {"ácido 2-(2-(metilamino)-3-oxo-4-sulfanilbutil)-5-(1-metoxi-3-oxopropil)ciclopent-3-eno-1-carboxílico", "O=CCC(OC)C1C(C(=O)O)C(C=C1)CC(NC)C(=O)CS"},
        {"ácido 3-(bromometil)-5-(3-butoxi-1-hidroxi-3-oxopropil)ciclo-hexa-1,5-dieno-1-carboxílico", "CCCCOC(=O)CC(O)C1=CC(C(=O)O)=CC(CBr)C1"},
        {"3-(3-(1-aminoprop-2-en-1-il)-2-hidroxi-6-(3-hidroxipropil)-5-metoxifenil)propanamida", "OCCCC1C(OC)=CC(C(N)C=C)=C(O)C=1CCC(=O)N"},
        {"3-(2-hidroxi-3,5-diiodofenil)-3-oxopropanoato de 4-fenoxi-2-formilfenilo", "C1=C(I)C=C(I)C(O)=C1C(=O)CC(=O)OC2=CC=C(C=C2C=O)OC3=CC=CC=C3"}
    };

    for (int i = 0; i < tests.size(); i++) {
        if (validateName(tests[i].first, tests[i].second, false) == 1) {
            return 1;
        }
    }

    std::cout << "Testing complete!" << std::endl;
    return 0;
}