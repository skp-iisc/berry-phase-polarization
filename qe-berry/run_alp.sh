#!/bin/bash

# INPUT pseudopotential, check from ../pseudo folder.
Al_PP="Al_ONCV_PBE-1.2.upf"
P_PP="P_ONCV_PBE-1.2.upf"

if [[ -d alp ]];
then
   echo "CALCULATIONS FOR AlP"
else
   echo "CALCULATIONS FOR AlP"
   mkdir alp
fi
cd alp

data_alat="alat.txt"
if [[ -f $data_alat ]]; then
    rm $data_alat
fi
echo "Lattice(Bohr)    Total Energy(Ry)" > $data_alat

for a in $(seq 10.21 0.01 10.30); do

    mkdir -p folder-$a
    cd folder-$a
    
    cat > alp_scf.in << EOF
&control
    calculation='scf',
    restart_mode='from_scratch',
    prefix='alp',
    pseudo_dir='../../../pseudo/',
    outdir='./'
 /
 &system
    ibrav= 1, celldm(1)=$a, nat= 8, ntyp= 2,
    ecutwfc = 25.0
 /
 &electrons
    diagonalization='david',
    conv_thr =  1.0d-8,
    mixing_beta = 0.5,
    startingwfc='random'
 /
ATOMIC_SPECIES
 Al  26.982 $Al_PP
 P  30.974 $P_PP
ATOMIC_POSITIONS alat
 Al -0.125 -0.125 -0.125
 Al  0.375  0.375 -0.125
 Al  0.375 -0.125  0.375
 Al -0.125  0.375  0.375
 P  0.125  0.125  0.125
 P  0.625  0.625  0.125
 P  0.625  0.125  0.625
 P  0.125  0.625  0.625
K_POINTS {automatic}
4 4 4 0 0 0
EOF
    echo "Running scf for celldm(1) = $a Bohr..."
    mpirun -np 4 pw.x < alp_scf.in > alp_scf.out

    tot_en=$(grep "!    total energy" alp_scf.out | tail -n 1 | awk '{print $5}')

    cd ../

    if [ -z "$tot_en" ]; then
        echo "Error: pw.x failed at a=$a. Skipping..."
    else
        echo "$a  $tot_en" >> $data_alat
    fi

done

echo "Check $data_alat for the results."
rm -r folder*

# Clamped-Ion SCF

cat > efield0.in << EOF
&control
    calculation='scf'
    restart_mode='from_scratch',
    prefix='alp',
    lelfield=.true.,
    nberrycyc=1
    pseudo_dir='../../pseudo/',
    outdir='../tmp/'
    tprnfor=.true.
 /
 &system
    ibrav= 1, celldm(1)=10.28, nat=  8, ntyp= 2,
    ecutwfc = 40.0
 /
 &electrons
    diagonalization='david',
    conv_thr =  1.0d-8,
    mixing_beta = 0.5,
    startingwfc='random',
    efield_cart(1)=0.0d0,efield_cart(2)=0.0d0,efield_cart(3)=0.0d0
 /
ATOMIC_SPECIES
 Al  26.982 $Al_PP
 P  30.974 $P_PP
ATOMIC_POSITIONS alat
 Al -0.125 -0.125 -0.125
 Al  0.375  0.375 -0.125
 Al  0.375 -0.125  0.375
 Al -0.125  0.375  0.375
 P  0.125  0.125  0.125
 P  0.625  0.625  0.125
 P  0.625  0.125  0.625
 P  0.125  0.625  0.625
K_POINTS {automatic}
4 4 4 0 0 0
EOF

echo "Running alp/efield0.in..."
mpirun -np 4 pw.x < efield0.in > efield0.out
echo "Done."

cat > efield1.in << EOF
 &control
    calculation='scf'
    restart_mode='from_scratch',
    prefix='alp',
    lelfield=.true.,
    nberrycyc=3
    pseudo_dir='../../pseudo/',
    outdir='../tmp/'
    tprnfor=.true.
 /
 &system
    ibrav= 1, celldm(1)=10.28, nat=  8, ntyp= 2
    ecutwfc = 40.0
 /
 &electrons
    diagonalization='david',
    conv_thr =  1.0d-8,
    mixing_beta = 0.5,
    startingwfc='random',
    efield_cart(1)=0.001d0,efield_cart(2)=0.0d0,efield_cart(3)=0.0d0
 /
ATOMIC_SPECIES
 Al  26.982 $Al_PP
 P  30.974 $P_PP
ATOMIC_POSITIONS alat
 Al -0.125 -0.125 -0.125
 Al  0.375  0.375 -0.125
 Al  0.375 -0.125  0.375
 Al -0.125  0.375  0.375
 P  0.125  0.125  0.125
 P  0.625  0.625  0.125
 P  0.625  0.125  0.625
 P  0.125  0.625  0.625
K_POINTS {automatic}
4 4 4 0 0 0
EOF

echo "Running alp/efield1.in..."
mpirun -np 4 pw.x < efield1.in > efield1.out
echo "Done."

# Relaxed-Ion Calculations

cat > relax_efield0.in << EOF
 &control
    calculation='relax'
    restart_mode='from_scratch',
    prefix='alp',
    lelfield=.true.,
    nberrycyc=1
    pseudo_dir='../../pseudo/',
    outdir='../tmp/'
    tprnfor=.true.
 /
 &system
    ibrav= 1, celldm(1)=10.28, nat=  8, ntyp= 2
    ecutwfc = 40.0
 /
 &electrons
    diagonalization='david',
    conv_thr =  1.0d-8,
    mixing_beta = 0.5,
    startingwfc='random',
    efield_cart(1)=0.0d0,efield_cart(2)=0.0d0,efield_cart(3)=0.0d0
 /
 &ions
    ion_dynamics='bfgs'
 /
ATOMIC_SPECIES
 Al  26.982 $Al_PP
 P  30.974 $P_PP
ATOMIC_POSITIONS alat
 Al -0.125 -0.125 -0.125
 Al  0.375  0.375 -0.125
 Al  0.375 -0.125  0.375
 Al -0.125  0.375  0.375
 P  0.125  0.125  0.125
 P  0.625  0.625  0.125
 P  0.625  0.125  0.625
 P  0.125  0.625  0.625
K_POINTS {automatic}
4 4 4 0 0 0
EOF

echo "Running alp/relax_efield0.in..."
mpirun -np 4 pw.x < relax_efield0.in > relax_efield0.out
echo "Done."

cat > relax_efield1.in << EOF
 &control
    calculation='relax'
    restart_mode='from_scratch',
    prefix='alp',
    lelfield=.true.,
    nberrycyc=3
    pseudo_dir='../../pseudo/',
    outdir='../tmp/'
    tprnfor=.true.
 /
 &system
    ibrav= 1, celldm(1)=10.28, nat=  8, ntyp= 2
    ecutwfc = 40.0
 /
 &electrons
    diagonalization='david',
    conv_thr =  1.0d-8,
    mixing_beta = 0.5,
    startingwfc='random',
    efield_cart(1)=0.001d0,efield_cart(2)=0.0d0,efield_cart(3)=0.0d0
 /
 &ions
    ion_dynamics='bfgs'
 /
ATOMIC_SPECIES
 Al  26.982 $Al_PP
 P  30.974 $P_PP
ATOMIC_POSITIONS alat
 Al -0.125 -0.125 -0.125
 Al  0.375  0.375 -0.125
 Al  0.375 -0.125  0.375
 Al -0.125  0.375  0.375
 P  0.125  0.125  0.125
 P  0.625  0.625  0.125
 P  0.625  0.125  0.625
 P  0.125  0.625  0.625
K_POINTS {automatic}
4 4 4 0 0 0
EOF

echo "Running alp/relax_efield1.in..."
mpirun -np 4 pw.x < relax_efield1.in > relax_efield1.out
echo "Done."

cd ../

rm -r tmp/
