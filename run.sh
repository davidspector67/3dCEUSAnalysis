SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
progPath="$SCRIPT_DIR/main.py"
source "$SCRIPT_DIR/venv/bin/activate"
python $progPath