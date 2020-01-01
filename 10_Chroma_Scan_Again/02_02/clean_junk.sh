#!/bin/bash

# clean ghost files
rm PTC/*~
rm Input/*~
rm lib/*~
rm lib/spacecharge_tunespread/*.pyc

# clean python files
rm lib/*.pyc
rm *.pyc

# ptc cleanup
rm junk.txt
rm Maxwellian_bend_for_ptc.txt
rm Negative_node.OUT
rm internal_mag_pot.txt
rm fort.18
