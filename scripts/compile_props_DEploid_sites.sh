#!/usr/bin/env bash
tail -n 1 run_DEploid/*.classic.prop |
sed -e 's/.*\///' -e 's/.classic.prop.*//' |
tr '\n' '\t' |
sed 's/\t\t/\n/g' > run_DEploid/filtered_non-swga_AF.DEploid.props
