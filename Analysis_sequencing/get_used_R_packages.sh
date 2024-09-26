#!/bin/bash
rm used_R_packages
grep library -h ./figures/*/*.r | tr '; ' '\n' | sed 's/\r//g' | grep "\S" | sort | uniq | grep -o '([^)]*)' | tr -d '()' > used_R_packages_pre
while read p; do echo -n $p; Rscript -e "as.character(packageVersion('$p'))" | sed s/\\[1\\]//g | sed s/\"//g; done < used_R_packages_pre > used_R_packages
rm used_R_packages_pre
