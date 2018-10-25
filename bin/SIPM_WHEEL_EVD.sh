

BIN_DIR=$( cd $(dirname ${BASH_SOURCE[0]}) && pwd )
WHEEL_DIR=$( dirname $BIN_DIR )
EVENT_DISPLAY_DIR="$WHEEL_DIR/eventDisplay"

#source $EVENT_DISPLAY_DIR/sipmwheel/config/setup.sh

echo
echo "Running event display..."
echo 

root -l $EVENT_DISPLAY_DIR/EventDisplay.C

