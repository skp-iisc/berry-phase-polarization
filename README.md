# Berry Phase Polarization

Calculation of Berry phase polarization using Quantum MASALA, Quantum ESPRESSO, and tight-binding model. 

## Using Quantum ESPRESSO (QE)

In the QE package, the example 10 is dedicated for this purpose.

Benchmarking: [Ref. 4].

## Using Quantum MASALA (QTM)

The inclusion of electric field is currently not present in this package. This work does that. 

Benchmarking: QE results, [Ref. 4] is used.

## Using tight-binding model

This deals with time dependent electric fields, but using a 3 band tight binding Hamiltonian [Ref. 2] from which dynamical Berry phase polarization is calculated.

Benchmarking: [Ref. 5].

Main references:

1. King-Smith and Vanderbilt; PRB 47, 1651 (1993).
2. Nunes and Vanderbilt; PRL 73, 712 (1994).
3. Nunes and Gonze; PRB 63, 155107 (2001).
4. Souza, Iniguez, and Vanderbilt; PRL 89, 117602 (2002).
5. Souza, Iniguez, and Vanderbilt; PRB 69, 085106 (2004).
