#!/bin/bash
cd "${0%/*}"
cd ../fwsw_results/dadi_results/
python ../../scripts/253_plot_best_dadi.py FLLG FLCC IM 0.1054 3.7014 0.5469 0.2949 0.1553
python ../../scripts/253_plot_best_dadi.py ALFW ALST SC2mG 2.5116 12.0516 19.5199 14.2786 0.5088 1.6496 1.8945 2.4924 1.4216 0.3354 0.6241
python ../../scripts/253_plot_best_dadi.py LAFW ALST SC2mG 1.7085 1.8422 11.3515 7.6879 0.0338 4.6925 0.0681 5.3442 0.5737 0.1248 0.9376
python ../../scripts/253_plot_best_dadi.py TXFW TXCC IM2m 0.3205 6.3369 9.0838 10.8524 1.3404 0.1808 0.6585 0.4870
