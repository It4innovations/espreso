mpirun -n 128 --bind-to core    --map-by core    espreso -c espreso.ecf SQUARE4  8 16   8  4  100 100
mpirun -n  32 --bind-to L3cache --map-by L3cache espreso -c espreso.ecf SQUARE4  4  8  16  8  100 100
mpirun -n   8 --bind-to numa    --map-by numa    espreso -c espreso.ecf SQUARE4  2  4  32 16  100 100
mpirun -n   2 --bind-to socket  --map-by socket  espreso -c espreso.ecf SQUARE4  1  2  64 32  100 100
mpirun -n   1 --bind-to none    --map-by node    espreso -c espreso.ecf SQUARE4  1  1  64 64  100 100
