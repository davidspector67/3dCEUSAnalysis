SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
virtualenv --python==python3.9 venv
source "$SCRIPT_DIR/venv/bin/activate"
pip install -r "$SCRIPT_DIR/pyPackages.txt"