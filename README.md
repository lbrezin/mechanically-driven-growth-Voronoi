This code was used for the simulation of the paper "Mechanically-driven growth and competition in a Voronoi model of tissues", by Louis Brezin and Kirill Korolev.

This code is almost identical to the cellGPU code developped by Daniel Sussman available at https://dmsussman.gitlab.io/cellGPUdocumentation/index.html 
We thank him for making the code available to the scientific community!

Some modifications were made to consider mechanically-driven growth in the main file cellDivision.cpp. 
Some of the src files were modified compared to the original code in order to:
- Export specific data from the simulation
- Allow different mechanical moduli for different population types
- Incorporate pressure computations within the code

The analysis folder contains files in python to analyse the data
