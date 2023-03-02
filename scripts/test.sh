
res=`find  tests/heatTransfer/ -name espreso.ecf | grep 2D`
res="${res}
`find  tests/linearElasticity/ -name espreso.ecf | grep plane`"
readarray -td$'\n' ecfs2D <<<"$res";

res=`find  tests/heatTransfer/ -name espreso.ecf | grep 3D`
res="${res}
`find  tests/linearElasticity/ -name espreso.ecf | grep volume`"
readarray -td$'\n' ecfs3D <<<"$res";

for ecf in "${ecfs3D[@]}"
do
   echo "${ecf}"
done