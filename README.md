# ttbbFourFermionEFT
code for EFT 4-fermion interpreation of the ttbb xsec measurement

## production
Instruction based on available installation of MG5_aMC_v2_6_0, Delphes 3.4.1 and Pythia8, using UFO file provide by TOP LHC WG:
```
www.desy.de/~durieux/topbasis/dim6top_LO_UFO.tar.gz.
```
Also includes other models personally created (4Heavy and 2Heavy2light)

Example Delphes_cards and run_cards are in the "cards" directory.

/production/write_restriction_cards_ttbb_fixed.py: creates output directories for each coupling and each scanned Wilson Coefficient, for parallel running on batch system. This creates .txt files that define the output directories and also include the proper restrictions of the model such that only 4-heavy-fermion operators are kept.

/production/generate_restricted_processes.py: generate the process for each of the creates directories (Probably large duplication of work, but convenient for parallel running on batch)

/production/launch_restricted_processes_ttbb_toDelphes.py: launch each of the processes to batch, making sure the parameters of the run_card and the param_card are properly set for each coupling and each Wilson coefficient. Also loads the Delphes card used in this analysis. NOTE: This code is specific to the VUB IIHE batch system!! This by default runs Delphes and outputs Delphes .root files. The final cross section are stored in a .txt file.

## conversion
/conversion/convertDelphesReconstructed_ttbb.py: converts the ROOT output of Delphes to a more flat root tree structure (no user defined classes) based on some specific reconstruction of the ttbb final state.

## crosssection
This directory contains code to derive limits based on the cross section as caluclated by MadGraph (does not rely on the conversion of the Delphes output!!!)

## m4b
Contains code to derive limits based on the M4b cut.

## NeuralNetwork
Code relies on an installation of Keras (Network training), Scikit-Learn (variable preprocessing) and root_numpy (conversion of ROOT tree to numpy arrays)
/NeuralNetwork/trainClassifier_ttbb.py: Neural network training.

## TemplateFits
Code that is used to perform 2D template fits (using ROOFit) on the Neural network outputs. Depends on an installation of Keras (Network training), Scikit-Learn (variable preprocessing) and root_numpy (conversion of ROOT tree to numpy arrays).




