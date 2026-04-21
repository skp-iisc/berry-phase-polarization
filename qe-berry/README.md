# Calculations on Quantum ESPRESSO

Calculations can be done for 3 different types of pseudopotentials, by giving input in the `.sh` file. The lattice optimization results is shown as `.png` files inside the sub-folders.

- Run `python calc_alas.py`. Output:
```
Check alas/alat.png.
Z* = 2.1105086808753
eps_inf = 8.424098614478504.
eps_stat = 10.007195153894429.
```

- Run `python calc_alp.py`. Output:
```
Check alp/alat.png.
Z* = 2.220569851366986
eps_inf = 7.489919283596497.
eps_stat = 9.405485152764415.
```

Benchmarking: *Souza, Iniguez, and Vanderbilt; PRL 89, 117602 (2002)*.
