Short SEDRs from the various Voyager encounters.
These files are documented in Voyager Memo 200, Revised June 20, 2022

The original SSEDR for Neptune that came with the PLS data is the "ssedr.v2.nep.txt". While the PLS data cover a period from 1989-233T19:55 to 1989-241T20:00, the Neptune SSEDR contained only 2 days of data, from 1989-236T00:00 to 1989-238T04:00. In order to fill in the gaps, the PLS SSEDR file was combined with the cruise SSEDR, "v2_shortsedr_1989.156.00_1989.275.00.Rob_(cruise_SEDR).csv". This resulted to a new file "ssedr.v2.nep.modified".

In the in the new SSEDR file we kept all the data of the original Neptune SSEDR file, and we expanded the rest as follow:
- We used the time for each measurement
- For the state of V2 on the "Neptune Centered, Neptune True Prime Meridian and Equator of Date" we calculated the state based on the best fit frame called "V2_NEPTUNE_PLS"
- For the V2 state on ECL50 we extrapolated the numbers from the cruise SSEDR using a linear approximation
- The "Cartesian Position of Triton, Neptune Centered, Neptune True Prime Meridian and Equator of Date" outside the PLS SSEDR has not been calculated
- The other three columns before the rotation matrix are there for the new file columns to agree with the old one - they are also not used in PLS SSEDR
- The SC->ECL50 rotattion matrix has been calculated based on the cruise SSEDR: the closest times were getting the specific value. E.g. if the cruise SSEDR had two rotation matrices at DOY238T00:00:00 and DOY238T02:00:00, all the times until 01:00:00 would get the first rotation matrix, and from 01:00 they would get the values of the second rotation matrix. (That's the best we can do, given that it is very complicated to extrapolate exact values for the rotation matrices. It is adequate though, as outside the PLS SSEDR there are no rotations while V2 is in the magnetosphere/magnetosheath)


Additionally there are two auxilliary files:
1) constructed.nep.rotation.sc.ecl50
It containes all the rotation matrices for the spacecraft frame to ECL50. While there are seemingly large datagaps (e.g. between DOY234 and 236, or 238 and 240), cruise SSEDR tends not to have multiple datapoints if Voyager does not rotate, i.e. Voyager did not rotate between the last datapoint in DOY234 to the first in DOY236. Therefore it can somewhat be safely assumed that the rotation matrix is the same for the time in-between.

2) constructed.nep.v2.state.ecl50
It contains Voyager 2 state in ECL50. The state in the period covered by the PLS SSEDR was kept as it is, but outside that period was linearly extrapolated based on the cruise SSEDR file. While this can potentially introduce an uncertainty in the location but also speed of Voyager (as the gravitational forces do not allow a perfectly linear movement, but also we are lacking information of boosters burning), the differences are adequately small that do not contribute almost any uncertainties during the analysis. (e.g. the speed and location between the gaps were fitted almost perfectly using a linear fitting).