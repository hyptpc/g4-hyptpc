--------------------------
| 2023.01.15             |
| Modified by Haein Lee  |
--------------------------


In directory "LED_modified" : LED data at the KEK 2nd building and at ELPH
In directory "simul_ELPH" : Geant4 simulation data of ELPH setup


gain_check.cc : Code for getting one photon gain (not perfect fitting result,,,because of the limited time.)


analysis_compare.cc
 - Pedestal fitting for Trigger counters (T4 - T7) 
 - Pedestal fitting for Individual channel and SUM channel of BAC   (With no beam data, using clock generator)
 - Conversion from ADC to Npe (I just brought one photon gain values from the "gain_check.cc". If you want to change this value, you have to modify by yourself... )
 - Cut condition using ADC & TDC of Trigger counters and TDC of BAC (no cut condition for ADC of BAC)
 - Npe distribution of BAC SUM Channel (There is no result of Individual channels.)
 - Npe distribution of Simulation  (SUM of four MPPCs)
 - Efficiency at each position

* You have to change a value "run_directory!! >>> modify to the right directory where the ELPH data exists.
* You can change the Y position. The detailed explanation about the position is in "position_information.png".


