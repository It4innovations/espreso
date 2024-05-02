
./waf configure
./waf $1
. tests/scripts/groups/feti.sh

./waf configure --intwidth=64
./waf $1
. tests/scripts/groups/feti.sh
