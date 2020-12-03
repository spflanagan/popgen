#!/bin/bash
cd "${0%/*}"
cd ../fwsw_results/dadi_results/
python ../../scripts/253_plot_best_dadi.py FLLG FLCC IM 0.1054 3.7014 0.5469 0.2949 0.1553
python ../../scripts/253_plot_best_dadi.py ALFW ALST SC2mG 1.9278 26.0526 2.7965 1.1743 28.4259 0.1592 9.9064 8.0659 0.8581 0.1573 0.5531
python ../../scripts/253_plot_best_dadi.py LAFW ALST SC2mG 0.8340 2.3917 11.3725 2.7406 0.8783 5.8332 3.8272 0.1672 0.4249 0.3549 0.9390
python ../../scripts/253_plot_best_dadi.py TXFW TXCC IM2m 0.3205 6.3369 9.0838 10.8524 1.3404 0.1808 0.6585 0.4870
