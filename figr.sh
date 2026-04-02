#!/bin/bash

# Check whether this script was called from figr's root directory.
if [ ! -f "$(pwd)/toolchain/util.sh" ]; then
    echo "figr: ERROR > You must call this script from within figr's root folder."

    exit 1
fi

# Load utility script
. "$(pwd)/toolchain/util.sh"

# Print startup message immediately for user feedback
log "Starting..."

# Handle upgrading from older figr build systems
if [ -d "$(pwd)/bootstrap" ] || [ -d "$(pwd)/dependencies" ] || [ -f "$(pwd)/build/figr.lock.yaml" ]; then
    error "Please remove, if applicable, the following directories:"
    error "* dependencies/"
    error "* bootstrap/"
    error "* build/"
    error "These are incompatible with modern versions of figr."

    exit 1
fi

# If the user wishes to run the "load" script
if [ "$1" '==' 'load' ] && [ "$2" != "--help" ] && [ "$2" != "-h" ]; then
    # Check if the script is being sourced (required for load to work)
    if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
        echo ""
        echo "figr: ERROR > The 'load' command must be run with 'source' to set environment variables."
        echo ""
        echo "  Instead of:  ./figr.sh load ..."
        echo "  Please use:  source ./figr.sh load ..."
        echo "          or:  . ./figr.sh load ..."
        echo ""
        echo "  Example:     source ./figr.sh load -c p -m g"
        echo ""
        exit 1
    fi
    shift; . "$(pwd)/toolchain/bootstrap/modules.sh" "$@"; return
elif [ "$1" '==' "format" ] && [ "$2" != "--help" ] && [ "$2" != "-h" ]; then
    . "$(pwd)/toolchain/bootstrap/python.sh" "$@"

    shift; . "$(pwd)/toolchain/bootstrap/format.sh" "$@"; exit 0
elif [ "$1" '==' "venv" ]; then
    shift; . "$(pwd)/toolchain/bootstrap/python.sh" "$@"; return
elif [ "$1" '==' "clean" ] && [ "$2" != "--help" ] && [ "$2" != "-h" ]; then
    rm -rf "$(pwd)/build"; exit 0
fi

mkdir -p "$(pwd)/build"

. "$(pwd)/toolchain/bootstrap/cmake.sh"
. "$(pwd)/toolchain/bootstrap/python.sh" "$@"

# init command: just bootstrap the environment and exit (no Python command)
if [ "$1" '==' 'init' ]; then
    exit 0
fi

echo

# Run the main.py bootstrap script
# If only flags given (no command), show help without passing flags
has_command=false
for arg in "$@"; do
    case "$arg" in -*) ;; *) has_command=true; break ;; esac
done
if [ "$has_command" = true ]; then
    python3 "$(pwd)/toolchain/main.py" "$@"
else
    python3 "$(pwd)/toolchain/main.py"
fi
code=$?

echo

if [ $code -ne 0 ]; then
    error "main.py finished with a $code exit code."
fi

# Deactivate the Python virtualenv in case the user "source"'d this script
log "(venv) Exiting the$MAGENTA Python$COLOR_RESET virtual environment."
deactivate

# Exit with proper exit code
exit $code
