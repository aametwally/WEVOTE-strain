#!/bin/bash

# Copyright 2013-2015, Derrick Wood <dwood@cs.jhu.edu>
#
# This file is part of the Kraken taxonomic classification system.
#
# Kraken is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Kraken is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Kraken.  If not, see <http://www.gnu.org/licenses/>.

set -e

VERSION="1.0"

if [ -z "$1" ] || [ -n "$2" ]
then
  echo "Usage: $(basename $0) KRAKEN_DIR"
  exit 64
fi

if [ "$1" = "KRAKEN_DIR" ]
then
  echo "Please replace \"KRAKEN_DIR\" with the name of the directory"
  echo "that you want to install Kraken in."
  exit 1
fi

# Perl cmd used to canonicalize dirname - "readlink -f" doesn't work
# on OS X.
export KRAKEN_DIR=$(perl -MCwd=abs_path -le 'print abs_path(shift)' "$1")

mkdir -p "$KRAKEN_DIR"
#Mladen's comment - here is where u make all the executables
make -C src install
#Mladen - This goes in each file, and substitutes KRAKEN_DIR for the actual directory
for file in scripts/*
do
  perl -pl -e 'BEGIN { while (@ARGV) { $_ = shift; ($k,$v) = split /=/, $_, 2; $H{$k} = $v } }'\
           -e 's/#####=(\w+)=#####/$H{$1}/g' \ #This is the part that gets deleted and replaced in each filen
           "KRAKEN_DIR=$KRAKEN_DIR" "VERSION=$VERSION" \
           < "$file" > "$KRAKEN_DIR/$(basename $file)"
  if [ -x "$file" ]
  then
    chmod +x "$KRAKEN_DIR/$(basename $file)"
  fi
done

echo
echo "Kraken installation complete."
echo
echo "To make things easier for you, you may want to copy/symlink the following"
echo "files into a directory in your PATH:"
for file in $KRAKEN_DIR/kraken*
do
  if [ -x "$file" ]
  then
    echo "  $file"
  fi
done