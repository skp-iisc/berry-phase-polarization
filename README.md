# Berry Phase Polarization

Calculation of Berry phase polarization using Quantum MASALA, Quantum ESPRESSO, and tight-binding model. 

## Using Quantum ESPRESSO (QE)

In the QE package, the example 10 is dedicated for this purpose. The results of [Ref. 4] is reproduced by using QE.

## Using Quantum MASALA (QTM)

The inclusion of electric field is currently not present in this package. This work intends to do that. For benchmarking the results, the results from QE is used.

## Using tight-binding model

This is based on a 3 band tight binding Hamiltonian [Ref. 2] from which dynamical Berry phase polarization is is calculated. The results matches with the results of [Ref. 5].


Main references:

1. R. D. King-Smith and David Vanderbilt; PRB 47, 1651 (1993).
2. R. W. Nunes and David Vanderbilt; PRL 73, 712 (1994).
3. R. W. Nunes and Xavier Gonze; PRB 63, 155107 (2001).
4. Ivo Souza, Jorge Iniguez, and David Vanderbilt; PRL 89, 117602 (2002).
5. Ivo Souza, Jorge Iniguez, and David Vanderbilt; PRB 69, 085106 (2004).
