# 2023_PNAS
The code was used in the following publication:

Fully-primed slowly-recovering vesicles mediate presynaptic LTP at neocortical neurons
Weichard I, Taschenberger H, Gsell F, Bornschein G, Ritzau-Jost A, Schmidt H, Kittel RJ, Eilers J, Neher E, Hallermann S, Nerlich J.
PNAS 2023 [doi: 10.1073/pnas.2305460120](https://doi.org/10.1073/pnas.2305460120)

A parallel and a sequential model were used. To execute the code you need to compile the C++ files in the folders /xcode/twopool and run main.cpp in the corresponding folders (parallel and sequential). You have to change the import and export folder names in line 10 and 11 of main.cpp to the appropriate paths on your computer. The xcode project file is also provided (twopool.xcodeproj). The code was run with the xcode application (version 15.4) on macOS 14.5.

First, the text file named apTimesSamp300.txt in folder /in is loaded. It contains three columns:  the time of stimuli (in ms), the weighting for the chi2 calculation, and a here not used factor for postsynaptic depression. Next, the text file 1_300Hz.txt in folder /in is loaded. The file contains six columns representing the EPSC amplitudes before and after LTP induction of three example experiments (cf. howManyExperiments = 6, line 20, main.cpp).

During execution, several text files are written in the folder /out of the corresponding folders (parallel and sequential). In the /out folders, the Mathematica files plot.nb (version 12, also provided as a PDF file plot.pdf) allows importing the text files and ploting the best-fit simulation results superimposed on the experimental data.

For more information contact hallermann@medizin.uni-leipzig.de
