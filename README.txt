This code contains a Matlab implementation of the Mirrortree-IPA as described in "Phylogenetic correlations can suffice to infer protein partners from sequences", by Guillaume Marmier, Martin Weigt and Anne-Florence Bitbol, DOI: 10.1101/670877.
Archived version: 10.5281/zenodo.3377592

Briefly, this is an algorithm iterative pairing algorithm (IPA) based on the Mirrortree method [1,2] that aims to predict interaction partners among paralogs from two protein families, just from their sequences. Here it is applied to our "standard dataset" of 5064 sequences of cognate histidine kinases and response regulators from the P2CS database (http://www.p2cs.org/).

This algorithm is a Mirrortree-based variant of the DCA-IPA introduced in "Inferring interaction partners from protein sequences", by Anne-Florence Bitbol, Robert S. Dwyer, Lucy J. Colwell, and Ned S. Wingreen, Proc. Natl. Acad. Sci. U.S.A 113 (43) 12180-12185 (2016), DOI: 10.1073/pnas.1606762113, and of the MI-IPA introduced in "Inferring interaction partners from protein sequences using mutual information", by Anne-Florence Bitbol, PLoS Comput Biol 14(11):e1006401 (2018), DOI: 10.1371/journal.pcbi.1006401.

In order to use the code, please run "Mirrortree_IPA_main" under Matlab.

The source code is freely available under the GNU GPLv3 license (unless otherwise indicated in specific files).

If you find this code useful for your research, please cite the associated reference, "Phylogenetic correlations can suffice to infer protein partners from sequences", by Guillaume Marmier, Martin Weigt and Anne-Florence Bitbol, DOI: 10.1101/670877.

References about the Mirrortree method:
[1] Pazos F, Valencia A. Similarity of phylogenetic trees as indicator of protein-protein interaction. Protein Eng Des Sel. 2001;14(9):609–614.
[2] Ochoa D, Pazos F. Studying the co-evolution of protein families with the Mirrortree web server. Bioinformatics. 2010;26(10):1370–1371, http://csbg.cnb.csic.es/mtserver.
