import pandas as pd
import csv
import pubchempy as pcp
from aimsim.ops.descriptor import Descriptor
from aimsim.chemical_datastructures.molecule import Molecule
from rdkit import Chem
from rdkit.Chem import MolFromSmiles
from aimsim.ops.similarity_measures import SimilarityMeasure

combosFile = open('drugcombs_scored.csv')
combosReader = csv.reader(combosFile)
combosData = list(combosReader)
# ['ID', 'Drug1', 'Drug2', 'Cell line', 'ZIP', 'Bliss', 'Loewe', 'HSA']
# ['1', '5-FU', 'ABT-888', 'A2058', '1.72', '6.26', '-2.75', '5.54']
combosFile2 = open('drug_chemical_info.csv', "r", encoding='latin1')
combosReader2 = csv.reader(combosFile2)
comboChemData = list(combosReader2)
# ['drugName', 'cIds', 'drugNameOfficial', 'molecularWeight', 'smilesString']
# ['Bendamustine', 'CIDs00065628', 'bendamustine', '358.26284', 'CN1C2=C(C=C(C=C2)N(CCCl)CCCl)N=C1CCCC(=O)O']
coconutFile = open('coconut_txt-2.csv')
coconutReader = csv.reader(coconutFile)
cocoData = list(coconutReader)
# ['O=C1OC2C(C(=C)C)CC1C3(O)CC4OC54C(=O)OC[CH]253C', 'CNP0000001']

#PART 4 testDataSet

print("starting part 4")
testDataSet = open('testDataSet.csv', 'w', newline='')
tdsWriter = csv.writer(testDataSet)
tdsWriter.writerow(["ID", "Smiles", "mol_formula", "tpsa","charge","mol_weight","hbdc","hbac", "rbc","hac" ])
#ADD ids then just do the attr
for cnpComp in cocoData:
    print(cnpComp[1])
    myComp = pcp.get_compounds(cnpComp[0], namespace='smiles', as_dataframe=False)
    attrlist1 = []
    attrlist1.append(cnpComp[1])
    attrlist1.append(cnpComp[0])
    attrlist1.append(myComp[0].molecular_formula)
    attrlist1.append(myComp[0].tpsa)
    attrlist1.append(myComp[0].charge)
    attrlist1.append(myComp[0].xlogp)
    attrlist1.append(myComp[0].molecular_weight)
    attrlist1.append(myComp[0].h_bond_donor_count)
    attrlist1.append(myComp[0].h_bond_acceptor_count)
    attrlist1.append(myComp[0].rotatable_bond_count)
    attrlist1.append(myComp[0].heavy_atom_count )
    tdsWriter.writerow(attrlist1)

testDataSet.close()

print("part 4 done")


#PART 3 trainingSet
smileComboFile = open('smileComboSet.csv')
scReader = csv.reader(smileComboFile)
smileCombosData = list(scReader)

print("starting part 3")
trainingSet = open('trainingSet.csv', 'w', newline='')
tsWriter = csv.writer(trainingSet)
titleRow= ["ID", "drug1","smile1", "drug2","smile2", 'ZIP', 'Bliss', 'Loewe', 'HSA']
titleRow = titleRow+ ["mol_formula1", "tpsa1","charge1","mol_weight1","hbdc1","hbac1", "rbc1","hac1"]
titleRow = titleRow+ ["mol_formula2", "tpsa2","charge2","mol_weight2","hbdc2","hbac2", "rbc2","hac2"]
tsWriter.writerow(titleRow)
# ADD ALL THE attr from both compounds + scores
for combo in smileCombosData:
    if(combo[0]=="skip"):
        tsWriter.writerow(["skip"])
        continue
    print(combo[0])
    scores = combosData[int(combo[0])]
    myComp = pcp.get_compounds(combo[2], namespace='smiles', as_dataframe=False)
    attrlist1 = []
    attrlist1.append(myComp[0].molecular_formula)
    attrlist1.append(myComp[0].tpsa)
    attrlist1.append(myComp[0].charge)
    attrlist1.append(myComp[0].xlogp)
    attrlist1.append(myComp[0].molecular_weight)
    attrlist1.append(myComp[0].h_bond_donor_count)
    attrlist1.append(myComp[0].h_bond_acceptor_count)
    attrlist1.append(myComp[0].rotatable_bond_count)
    attrlist1.append(myComp[0].heavy_atom_count )
    myComp2 = pcp.get_compounds(combo[4], namespace='smiles', as_dataframe=False)
    attrlist2 = []
    attrlist2.append(myComp2[0].molecular_formula)
    attrlist2.append(myComp2[0].tpsa)
    attrlist2.append(myComp2[0].charge)
    attrlist2.append(myComp2[0].xlogp)
    attrlist2.append(myComp2[0].molecular_weight)
    attrlist2.append(myComp2[0].h_bond_donor_count)
    attrlist2.append(myComp2[0].h_bond_acceptor_count)
    attrlist2.append(myComp2[0].rotatable_bond_count)
    attrlist2.append(myComp2[0].heavy_atom_count) 
    finalRow = combo+scores[4:8]+attrlist1+attrlist2;
    tsWriter.writerow(finalRow)

trainingSet.close()

print("part 3 done")

#PART 2 validationSet

# smileComboFile = open('smileComboSet.csv')
# scReader = csv.reader(smileComboFile)
# smileCombosData = list(scReader)
#
# print("starting")
# validationSet = open('validationSet.csv', 'w', newline='')
# vsWriter = csv.writer(validationSet)
# for combo in smileCombosData:
#     if(combo[0]=="skip"):
#         continue
#     print(combo[0])
#     check1 = False
#     check2 = False
#     for nut in cocoData:
#         if combo[2]==nut[0]: 
#             check1=True
#             print("matchone "+combo[0])
#         if combo[4]==nut[0]: 
#             check2=True
#             print("matchtwo "+combo[0])
#     if check1==True and check2==True:
#         vsWriter.writerow(comboChemData[int(combo[0])])
#
# validationSet.close()
#
# print("part 2 done")

#PART 1 SmilesComboSet

# smileComboSet = open('smileComboSet.csv', 'w', newline='')
# scWriter = csv.writer(smileComboSet)
# for combo in combosData:
#     smile1 = ""
#     smile2 =""
#     for drug in comboChemData:
#         if combo[1]==drug[0].upper():
#             smile1=drug[4]
#         if combo[2]==drug[0].upper():
#             smile2=drug[4]
#     if smile1!="" and smile2!="":
#         scWriter.writerow([combo[0],combo[1],smile1, combo[2],smile2])
#     else:
#         scWriter.writerow(["skip",combo[1],smile1, combo[2],smile2])
#
# smileComboSet.close()
# print("part 1 done")


