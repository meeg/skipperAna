#!/bin/bash
#./spurious_charge.py -e 60 plots_12 skp_1hr_*_12.root
#./spurious_charge.py -e 60 plots_13 skp_1hr_*_13.root
##./spurious_charge.py plots_14 skp_1hr_*_14.root
#./spurious_charge.py -e 60 plots_15 skp_1hr_*_15.root
#
#./spurious_charge.py -e 120 plots_2hr_12 skp_2hr_*_12.root
#./spurious_charge.py -e 120 plots_2hr_13 skp_2hr_*_13.root
#./spurious_charge.py -e 120 plots_2hr_15 skp_2hr_*_15.root

#./spurious_charge.py -e 120 plots_noHclks_12 skp_noHclks_*_12.root
#./spurious_charge.py -e 120 plots_noHclks_13 skp_noHclks_*_13.root
#./spurious_charge.py -e 120 plots_noHclks_15 skp_noHclks_*_15.root

#./spurious_charge.py -e 120 plots_vddOFF_2hr_12 20MAY2019_DarkCurrent_vddOFF/skp_2hr_*_12.root
#./spurious_charge.py -e 120 -g 250 plots_vddOFF_2hr_13 20MAY2019_DarkCurrent_vddOFF/skp_2hr_*_13.root
#./spurious_charge.py -e 120 plots_vddOFF_2hr_15 20MAY2019_DarkCurrent_vddOFF/skp_2hr_*_15.root

./spurious_charge.py -e 120 plots_vddtestOFF_2hr_12 23MAY2019_DarkCurrent_VDDtest/skp_2hr_vddOFF_*_12.root
./spurious_charge.py -e 120 -g 250 plots_vddtestOFF_2hr_13 23MAY2019_DarkCurrent_VDDtest/skp_2hr_vddOFF_*_13.root
./spurious_charge.py -e 120 plots_vddtestOFF_2hr_15 23MAY2019_DarkCurrent_VDDtest/skp_2hr_vddOFF_*_15.root

./spurious_charge.py -e 120 plots_vddtestON_2hr_12 23MAY2019_DarkCurrent_VDDtest/skp_2hr_vddON_*_12.root
./spurious_charge.py -e 120 -g 250 plots_vddtestON_2hr_13 23MAY2019_DarkCurrent_VDDtest/skp_2hr_vddON_*_13.root
./spurious_charge.py -e 120 plots_vddtestON_2hr_15 23MAY2019_DarkCurrent_VDDtest/skp_2hr_vddON_*_15.root
