#!/bin/bash

# Base directory of the Intratumoral folder
INTRATUMORAL_DIR="/home/rungger/.conda/envs/Intratumoral"

# Loop through all subdirectories in Intratumoral
for env_dir in "$INTRATUMORAL_DIR"/env-*; do
    if [ -d "$env_dir" ]; then
        echo "Deleting environment: $env_dir"
        conda env remove -p "$env_dir" -y
    fi
done

echo "All unwanted environments deleted."
