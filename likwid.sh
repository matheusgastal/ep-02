#!/bin/bash

# Arquivo de entrada do resolveEDO
INPUT="teste.dat"  # substitua pelo seu arquivo ou remova se for digitar manualmente

# Executa o resolveEDO normalmente e mostra a saída
./resolveEDO < "$INPUT"

# Captura apenas o contador FP_ARITH_INST_RETIRED_SCALAR_DOUBLE usando LIKWID
sudo likwid-perfctr -g FLOPS_DP --marker GaussSeidel ./resolveEDO < "$INPUT" | \
awk '/FP_ARITH_INST_RETIRED_SCALAR_DOUBLE/ {print $1","$2}'
